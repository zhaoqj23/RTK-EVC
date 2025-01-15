#include "artest.h"

using namespace std;
using namespace Eigen;


/** 
 * @brief Test satellite system (m=0:GPS/SBS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN).
 * 
 * @param sys The satellite system identifier.
 * @param m The expected system identifier (0: GPS/SBS, 1: GLO, 2: GAL, 3: BDS, 4: QZS, 5: IRN).
 * 
 * @return Returns 1 if the system matches, otherwise returns 0.
 */
static int test_sys(int sys, int m)
{
	switch (sys) {
	case SYS_GPS: return m == 0;
	case SYS_SBS: return m == 0;
	case SYS_GLO: return m == 1;
	case SYS_GAL: return m == 2;
	case SYS_CMP: return m == 3;
	case SYS_QZS: return m == 4;
	case SYS_IRN: return m == 5;
	}
	return 0;
}

/** 
 * @brief Calculates the variance.
 * 
 * @param sat The satellite identifier.
 * @param sys The satellite system identifier.
 * @param el The elevation angle of the satellite.
 * @param bl The baseline length.
 * @param dt The differential age.
 * @param f The carrier frequency identifier.
 * @param opt Pointer to the options.
 * 
 * @return The calculated variance.
 */
static double varerr(int sat, int sys, double el, double bl, double dt, int f,
                     const prcopt_t *opt)
{
    double a,b,c=opt->err[3]*bl/1E4,d=CLIGHT*opt->sclkstab*dt,fact=1.0;
    double sinel=sin(el);
    int nf=(opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf;
    
    if (f>=nf) fact=opt->eratio[f-nf];
    if (fact<=0.0) fact=opt->eratio[0];
    fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    a=fact*opt->err[1];
    b=fact*opt->err[2];
    
    return 2.0*(opt->ionoopt==IONOOPT_IFLC?3.0:1.0)*(a*a+b*b/sinel/sinel+c*c)+d*d;
}


/** 
 * @brief Partial ambiguity resolution for the node.
 * 
 * @param graph Pointer to the Ceres problem object.
 * @param node Pointer to the node information.
 * 
 * @return The number of ambiguities fixed.
 */
int resamb_partial(ceres::Problem* graph, node_info *node) {
	int nfix, nsdambi, col, row;
	double minweight = 9999.9;
	int minweight_index = -1;
	set<int> ambi_index_tofix;
	vector<int> sdambi_index_inorder, sdambi_index_inorder_origin;
	MatrixXd Qac, Qc, Qc_origin, Qac_origin;
	double *bias = mat(node->nb, 1);
	int *ix = imat(node->nb, 2);
	fillup_amb_weight_flag(node); // Fill the weights

	find_ambiindex_tofix(ambi_index_tofix, node);
	nfix = select_ddidx(node, ambi_index_tofix, ix);
	if (nfix <= 0) return 0;

	// Sort the ambiguities
	nsdambi = sdambi_inorder(ix, nfix, sdambi_index_inorder);
	sdambi_index_inorder_origin = sdambi_index_inorder;

	// Initialize the covariance matrix
	Qc = MatrixXd::Zero(nsdambi, nsdambi); Qac = MatrixXd::Zero(3, nsdambi);
	Qc_origin = MatrixXd::Zero(nsdambi, nsdambi); Qac_origin = MatrixXd::Zero(3, nsdambi);

	// Estimate the covariance matrix
	if (!covariance_estimation(graph, Qc, Qac, sdambi_index_inorder, node)) {
		nfix = 0;
		return nfix;
	}
	Qc_origin = Qc; Qac_origin = Qac;

	while ((nfix = ar_lambda(graph, node, bias, ix, nfix, Qc, Qac, sdambi_index_inorder)) <= 0) {

		minweight = 9999.9;
		minweight_index = -1;
		// Exclude the ambiguity with the smallest weight
		for (auto index : ambi_index_tofix) {
			if (sdambi_info_list[index]->weight < minweight) {
				minweight = sdambi_info_list[index]->weight;
				minweight_index = index;
			}
		}

		if (minweight_index >= 0) {
			sdambi_info_list[minweight_index]->flag = 1;
			find_ambiindex_tofix(ambi_index_tofix, node);
			nfix = select_ddidx(node, ambi_index_tofix, ix);
			if (nfix <= MINNB) return 0;
			nsdambi = sdambi_inorder(ix, nfix, sdambi_index_inorder);
			Qc.resize(nsdambi, nsdambi); Qac.resize(3, nsdambi);
			Qc.setZero(); Qac.setZero();
			for (int i = 0; i < nsdambi; i++) {
				row = find_index(sdambi_index_inorder_origin, sdambi_index_inorder[i]);
				Qac.col(i) = Qac_origin.col(row);
				for (int j = 0; j < nsdambi; j++) {
					col = find_index(sdambi_index_inorder_origin, sdambi_index_inorder[j]);
					Qc(i, j) = Qc_origin(row, col);
				}
			}
			trace(4, "Qc = \n"); tracemat(4, Qc.data(), sdambi_index_inorder.size(), sdambi_index_inorder.size(), 3, 10);
			trace(4, "Qac = \n"); tracemat(4, Qac.data(), 3, sdambi_index_inorder.size(), 3, 10);
		}
	}

	free(bias); free(ix);
	return nfix;
}


/** 
 * @brief Finds the index of a given ambiguity in the original ambiguity list.
 * 
 * @param sdambi_index_inorder_origin The original list of ambiguity indices.
 * @param sdambi_index The ambiguity index to find.
 * 
 * @return The index of the ambiguity in the original list, or -1 if not found.
 */
int find_index(vector<int>& sdambi_index_inorder_origin, int sdambi_index) {
	auto it = find(sdambi_index_inorder_origin.begin(), sdambi_index_inorder_origin.end(), sdambi_index);
	if (it != sdambi_index_inorder_origin.end()) {
		return distance(sdambi_index_inorder_origin.begin(), it);
	}
	return -1;
}


/** 
 * @brief Orders the ambiguities in the node based on the given ambiguity indexes.
 * 
 * @param ix The array of ambiguity indexes.
 * @param nfix The number of ambiguities to be fixed.
 * @param sdambi_index_inorder The ordered ambiguity index list to be filled.
 * 
 * @return The number of ordered ambiguities.
 */
int sdambi_inorder(int *ix, int nfix, vector<int> &sdambi_index_inorder) {
	int nsdambi, i;
	sdambi_index_inorder.clear();
	for (i = 0; i < nfix; i++) {
		if (find(sdambi_index_inorder.begin(), sdambi_index_inorder.end(), ix[2 * i]) == sdambi_index_inorder.end()) {
			sdambi_index_inorder.push_back(ix[2 * i]);
		}
		if (find(sdambi_index_inorder.begin(), sdambi_index_inorder.end(), ix[2 * i + 1]) == sdambi_index_inorder.end()) {
			sdambi_index_inorder.push_back(ix[2 * i + 1]);
		}
	}
	nsdambi = sdambi_index_inorder.size();
	return nsdambi;
}


/** 
 * @brief Finds the ambiguity indexes to be fixed from the related ambiguity list.
 * 
 * @param ambi_index_tofix The set to hold the indexes of the ambiguities to be fixed.
 * @param node Pointer to the node information.
 */
void find_ambiindex_tofix(set<int> &ambi_index_tofix, node_info *node) {
	int i;
	ambi_index_tofix.clear();
	for (i = 0; i < node->nb; i++) {
		if (sdambi_info_list[node->ddfix[i].ref_sdambi_index]->flag > 1) {
			ambi_index_tofix.insert(node->ddfix[i].ref_sdambi_index);
		}
		if (sdambi_info_list[node->ddfix[i].uref_sdambi_index]->flag > 1) {
			ambi_index_tofix.insert(node->ddfix[i].uref_sdambi_index);
		}
	}
}



/** 
 * @brief Performs ambiguity resolution using the lambda method.
 * 
 * 
 * @param graph Pointer to the Ceres problem object.
 * @param node Pointer to the node information.
 * @param bias Array to hold the integer ambiguities.
 * @param ix Array containing the indexes of the ambiguities to be fixed.
 * @param nfix The number of ambiguities to be fixed.
 * @param Qc The covariance matrix of the ambiguities.
 * @param Qac The cross-covariance matrix between the ambiguities and the position.
 * @param sdambi_index_inorder List of ordered ambiguity indexes.
 * 
 * @return The number of ambiguities that were successfully fixed. Returns 0 if the resolution fails.
 */
int ar_lambda(ceres::Problem *graph, node_info *node, double *bias, int *ix, int nfix, MatrixXd &Qc, MatrixXd &Qac, vector<int> sdambi_index_inorder) {
    int i, j, nsdambi, index, info;
    double s[2];
    MatrixXd y, b, db, Qab, QQ, D, Qb, xaa, Paa, xf, Pf;

    nsdambi = sdambi_index_inorder.size();

    // Initialize matrices and vectors
    y = MatrixXd::Zero(nfix, 1);
    Qb = MatrixXd::Zero(nfix, nfix); 
    Qab = MatrixXd::Zero(3, nfix);
    D = MatrixXd::Zero(nfix, nsdambi); 
    b = MatrixXd::Zero(nfix, 2);
    xaa = MatrixXd::Zero(3, 1); 
    Paa = MatrixXd::Zero(3, 3);
    xf = MatrixXd::Zero(3, 1); 
    Pf = MatrixXd::Zero(3, 3);
    
    for (i = 0; i < 3; i++) {
        xf(i, 0) = node->xf[i];
        for (j = 0; j < 3; j++) {
            Pf(i, j) = node->Pf[i + node->na * j];
        }
    }

    // Construct the single-differencial matrix D
    for (i = 0; i < nfix; i++) {
        index = find(sdambi_index_inorder.begin(), sdambi_index_inorder.end(), ix[2 * i]) - sdambi_index_inorder.begin();
        D(i, index) = 1.0;
        index = find(sdambi_index_inorder.begin(), sdambi_index_inorder.end(), ix[2 * i + 1]) - sdambi_index_inorder.begin();
        D(i, index) = -1.0;
    }

    trace(4, "D = \n"); tracemat(4, D.data(), nfix, nsdambi, 3, 5);

    // Calculate the float DD ambiguities vector
    for (i = 0; i < nfix; i++) {
        y(i, 0) = sdambi_info_list[ix[2 * i]]->sdbias - sdambi_info_list[ix[2 * i + 1]]->sdbias;
    }
    trace(4, "y = \n"); tracemat(4, y.data(), 1, nfix, 5, 5);

    // Calculate Qb, which is the covariance matrix of the DD ambiguities
    Qb = D * Qc * (D.transpose());
    trace(4, "Qb = \n"); tracemat(4, Qb.data(), nfix, nfix, 3, 10);

    // Calculate the cross-covariance matrix Qab
    Qab = Qac * (D.transpose());
    trace(4, "Qab = \n"); tracemat(4, Qab.data(), 3, nfix, 3, 10);

    // Resolve the ambiguities through the lambda method
    if (!(info = lambda(nfix, 2, y.data(), Qb.data(), b.data(), s))) {
        
        // Log the results of the lambda method
        trace(4, "N(1)="); tracemat(4, b.data(), 1, nfix, 5, 3);
        trace(4, "N(2)="); tracemat(4, b.data() + nfix, 1, nfix, 5, 3);

        node->ratio = s[0] > 0 ? (float)(s[1] / s[0]) : 0.0f;
        if (node->ratio > 999.9) node->ratio = 999.9f;

        /* Validation by the popular ratio-test */
        if (s[0] <= 0.0 || s[1] / s[0] >= option.thresar[0]) {

            // Restore the fixed ambiguities
            for (i = 0; i < nfix; i++) {
                bias[i] = b(i, 0);
                y(i, 0) -= b(i, 0);
            }

            // Calculate the position and covariance
            Eigen::FullPivLU<Eigen::MatrixXd> lu(Qb);
            if (lu.isInvertible()) {
                Eigen::MatrixXd Qbinv = Qb.inverse();
                xaa = xf - Qab * Qbinv * y;
                Paa = Pf - Qab * Qbinv * (Qab.transpose());

                // DR consistency check
                if (consistency_check(xaa)) {
                    trace(3, "resamb : validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                          nfix, s[0] == 0.0 ? 0.0 : s[1] / s[0], s[0], s[1]);

                    // Restore SD ambiguity and update node's state
                    sdambi_restore(nfix, ix, bias, node);
                    memcpy(node->xa, xaa.data(), sizeof(double) * 3);
                    memcpy(node->xf, xaa.data(), sizeof(double) * 3);
                    for (i = 0; i < 3; i++) {
                        for (j = 0; j < 3; j++) {
                            node->Pf[i + node->na * j] = Paa(i, j);
                            node->Pa[i + node->na * j] = Paa(i, j);
                        }
                    }
                } else {
                    nfix = 0;
                }
                
            } else {
                nfix = 0;
            } 
        } else { /* Validation failed */
            trace(3, "ambiguity validation failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                  nfix, s[1] / s[0], s[0], s[1]);
            nfix = 0;
        }
    } else {
        trace(3, "lambda error (info=%d)\n", info);
        nfix = 0;
    }

    return nfix;
}


/** 
 * @brief Performs consistency check on the fixed solution using Dead Reckoning.
 * 
 * 
 * @param xaa The fixed solution to be checked.
 * 
 * @return Returns `true` if the consistency check passes, i.e., the error is below the threshold. 
 *         Returns `false` if the consistency check fails.
 */
bool consistency_check(MatrixXd &xaa) {
    size_t i, j;
    double weight = 0.0, total_weight = 0.0;
    std::vector<size_t> fix_index;
    Vector3d xa = Vector3d::Zero();

    // Find all fixed epochs in a certain period.
    for (i = node_list.size() - CONP; i < node_list.size() - 1; i++) {
        if (node_list[i]->stat == SOLQ_FIX) {
            bool iscon = true;
            // Check if the epochs have interruptions.
            for (j = i; j < node_list.size() - 1; j++) {
                if (node_list[j + 1]->nepoch - node_list[j]->nepoch > 1) {
                    iscon = false;
                    break;
                }
            }
            if(iscon) {
                fix_index.push_back(i);
            } 
        }
    }

    // If no fixed nodes found, abandon the consistency check.
    if (fix_index.size() == 0) {
        return true;
    }

    // Check the consistency by Dead Reckoning
    for (auto index : fix_index) {
        Vector3d xf = Vector3d::Zero();
        memcpy(xf.data(), node_list[index]->xf, sizeof(double) * 3);
        Dead_Reckoning(index, xf.data());
        weight = 1.0 / (node_list.size() - index);
        xa += xf * weight;
        total_weight += weight;
    }

    // Calculate the weighted average
    if (total_weight > 0) {
        xa /= total_weight;
    }
    
    // Calculate the consistency error (distance between fixed and calculated solution)
    double conerr = (xa - xaa).norm();
    if (conerr < CON_THRE) {
        return true;  // Consistency check passed
    } else {
        trace(3, "Consistency check failed\n");
        return false;  // Consistency check failed
    }
}


/** 
 * @brief Performs Dead Reckoning to estimate the position from the given index.
 * 
 * 
 * @param index The starting index of the node from which Dead Reckoning begins.
 * @param xf The DR result
 */
void Dead_Reckoning(int index, double *xf) {
    size_t i;
    int j;
    double va[3] = { 0.0 }; 

    for (i = index; i < node_list.size() - 1; i++) {
        for (j = 0; j < 3; j++) {
            va[j] = (node_list[i]->xf[j + 3] + node_list[i + 1]->xf[j + 3]) / 2;  
            xf[j] += va[j] * node_list[i + 1]->tt;  // Update position estimate by integrating velocity over time
        }
    }
}

/** 
 * @brief Restores the SD ambiguities.
 * 
 * @param nfix The number of fixed ambiguities.
 * @param ix The index array that contains pairs of ambiguity indexes to be restored.
 * @param bias The array containing the integer ambiguities.
 * @param node Pointer to the node information.
 */
void sdambi_restore(int nfix, int *ix, double *bias, node_info *node) {
    int i;
    
    // Loop through each ambiguity pair and restore the SD ambiguity bias and flag
    for (i = 0; i < nfix; i++) {
        sdambi_info_list[ix[2 * i + 1]]->sdbias = sdambi_info_list[ix[2 * i]]->sdbias - bias[i];
        sdambi_info_list[ix[2 * i]]->flag = 3;
        sdambi_info_list[ix[2 * i + 1]]->flag = 3;
    }

    // Reset the node's double-differenced (ddfix) information
    for (i = 0; i < node->nb; i++) {
        // Reset all fields for each ddfix entry
        node->ddfix[i].refsat = 0;
        node->ddfix[i].ref_sdambi_index = -1;
        node->ddfix[i].urefsat = 0;
        node->ddfix[i].uref_sdambi_index = -1;
        node->ddfix[i].f = -1;
        node->ddfix[i].code = -1;
        node->ddfix[i].ref_index = -1;
        node->ddfix[i].uref_index = -1;
        node->ddfix[i].isValid = true;
        node->ddfix[i].isPrior = true;
        node->ddfix[i].var = 0.0;
        node->ddfix[i].ddbias = 0.0;
    }

    node->nb = 0;  // Reset the number of integer ambiguities

    for (i = 0; i < nfix; i++) {
        node->ddfix[i].code = 0;
        node->ddfix[i].f = sdambi_info_list[ix[2 * i]]->f;
        node->ddfix[i].refsat = sdambi_info_list[ix[2 * i]]->sat;
        node->ddfix[i].ref_sdambi_index = ix[2 * i];

        node->ddfix[i].f = sdambi_info_list[ix[2 * i + 1]]->f;
        node->ddfix[i].urefsat = sdambi_info_list[ix[2 * i + 1]]->sat;
        node->ddfix[i].uref_sdambi_index = ix[2 * i + 1];

        node->ddfix[i].ddbias = bias[i];
    }

    node->nb = nfix; 
}


/** 
 * @brief Estimates the covariance matrix.
 * 
 * @param graph Pointer to the Ceres problem object.
 * @param Qc Output covariance matrix for SD ambiguities.
 * @param Qac Output covariance matrix between the position and SD ambiguity.
 * @param sdambi_index_inorder Ordered list of SD ambiguity indexes.
 * @param node Pointer to the node information.
 * @return true if covariance computation is successful, false otherwise.
 */
bool covariance_estimation(ceres::Problem* graph, MatrixXd &Qc, MatrixXd &Qac, std::vector<int> &sdambi_index_inorder, node_info* node) {
	size_t i, j;
	double cov[1];
	ceres::Covariance::Options cov_options;
	cov_options.sparse_linear_algebra_library_type = ceres::EIGEN_SPARSE;
	cov_options.algorithm_type = ceres::DENSE_SVD;
	cov_options.null_space_rank = -1;
	cov_options.num_threads = 8;
	cov_options.apply_loss_function = false;
	ceres::Covariance covariance(cov_options);

	std::vector<std::pair<const double*, const double*>> covariance_blocks;

	// Add SD ambiguity covariance blocks
	for (i = 0; i < sdambi_index_inorder.size(); i++) {
		for (j = i; j < sdambi_index_inorder.size(); j++) {
			covariance_blocks.push_back(make_pair(&(sdambi_info_list[sdambi_index_inorder[i]]->sdbias), &(sdambi_info_list[sdambi_index_inorder[j]]->sdbias)));
		}
	}

	// Add covariance blocks between position and SD ambiguity
	for (i = 0; i < sdambi_index_inorder.size(); i++) {
		covariance_blocks.push_back(make_pair(&(node->xf[0]), &(sdambi_info_list[sdambi_index_inorder[i]]->sdbias)));
	}

	// Add position covariance block
	covariance_blocks.push_back(make_pair(&(node->xf[0]), &(node->xf[0])));

	// Compute covariance
	if (!covariance.Compute(covariance_blocks, graph)) {
		trace(3, "Error: Covariance computation failed.\n");
		return false;
	}

	// Fill the Qc matrix (SD ambiguity covariance)
	for (i = 0; i < sdambi_index_inorder.size(); i++) {
		for (j = i; j < sdambi_index_inorder.size(); j++) {
			covariance.GetCovarianceBlock(&(sdambi_info_list[sdambi_index_inorder[i]]->sdbias), &(sdambi_info_list[sdambi_index_inorder[j]]->sdbias), cov);
			Qc(i, j) = cov[0];
			Qc(j, i) = Qc(i, j);
		}
	}
	trace(4, "Qc = \n"); tracemat(4, Qc.data(), sdambi_index_inorder.size(), sdambi_index_inorder.size(), 3, 10);

	// Fill the Qac matrix (covariance between position and SD ambiguity)
	Eigen::VectorXd covariance_xambi(3 * 1);
	for (i = 0; i < sdambi_index_inorder.size(); i++) {
		covariance.GetCovarianceBlock(&(node->xf[0]), &(sdambi_info_list[sdambi_index_inorder[i]]->sdbias), covariance_xambi.data());
		Qac.col(i) = covariance_xambi;
	}
	trace(4, "Qac = \n"); tracemat(4, Qac.data(), 3, sdambi_index_inorder.size(), 3, 10);

	// Get position covariance (Pf)
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Pf(3, 3);
	Pf.setIdentity();
	covariance.GetCovarianceBlock(&(node->xf[0]), &(node->xf[0]), Pf.data());
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			node->Pf[i + node->na * j] = Pf(i, j);
		}
	}
	trace(4, "Pf = \n"); tracemat(4, Pf.data(), 3, 3, 3, 10);

	return true;
}


/** 
 * @brief Selects SD ambiguity pairs to fix based on weight and system type.
 * 
 * @param node Pointer to the node information.
 * @param ambi_index_tofix Set of SD ambiguity indices to fix.
 * @param ix Array to store selected ambiguity pairs.
 * @return The number of selected ambiguity pairs.
 */
int select_ddidx(node_info *node, set<int> &ambi_index_tofix, int *ix) {
    int nfix = 0, m, f, sat, maxweight = -1, maxindex = -1, i;

    // Initialize ix to -1
    for (i = 0; i < node->nb; i++) {
        ix[2 * i] = -1;
        ix[2 * i + 1] = -1;
    }

    // Loop through system types and frequencies to select ambiguity pairs
    for (m = 0; m < 6; m++) {
        for (f = 0; f < option.nf; f++) {

            vector<int> fixindex;
            for (auto index : ambi_index_tofix) {
                sat = sdambi_info_list[index]->sat;
                if (test_sys(node->ssat[sat - 1].sys, m) && f == sdambi_info_list[index]->f) {
                    fixindex.push_back(index);
                }
            }

            // Select the reference satellite with the highest weight
            if (fixindex.size() > 1) {
                maxweight = -1; maxindex = -1;
                for (auto index : fixindex) {
                    if (sdambi_info_list[index]->weight > maxweight) {
                        maxweight = sdambi_info_list[index]->weight;
                        maxindex = index;
                    }
                }

                for (auto index : fixindex) {
                    if (index == maxindex) continue;
                    ix[nfix * 2] = maxindex;
                    ix[nfix * 2 + 1] = index;
                    nfix++;
                }
            }
        }
    }

    return nfix;
}

/** 
 * @brief Fills up the weight for each ambiguity based on satellite system and frequency.
 * 
 * @param node Pointer to the node information.
 */
void fillup_amb_weight_flag(node_info *node) {
    int i;
    set<int> fixsatindex;

    // Collect ambiguity indices to be fixed
    for (i = 0; i < node->nb; i++) {
        fixsatindex.insert(node->ddfix[i].ref_sdambi_index);
        fixsatindex.insert(node->ddfix[i].uref_sdambi_index);
    }

    // Calculate and assign weight for each ambiguity
    for (auto index : fixsatindex) {
        int sat = sdambi_info_list[index]->sat;
        int f = sdambi_info_list[index]->f;
        sdambi_info_list[index]->weight = 0.001 / varerr(sat, node->ssat[sat - 1].sys, node->ssat[sat - 1].azel[1], 
        node->bl, node->dt, f, &option);
    }
}


