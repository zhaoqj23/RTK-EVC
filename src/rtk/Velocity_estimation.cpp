#include "../global/global.h"
#include "../factors/dopp_factor.h"
#include "../global/marginalization_info.h"
#include "../global/marginalization_factor.h"
#include "../factors/cvh_factor.h"
#include "../factors/cah_factor.h"
#include "../rtklib/rtklib.h"
#include "../global/MinSet_Generate.h"
#include <random>
#include <set>

#define MAXITR 10

/** 
 * @brief Performs the VFG optimization process.
 * 
 */
void VFG_optimization() {

    // Solver configuration
    ceres::Problem::Options problem_options;
    problem_options.enable_fast_removal = true; // Enable fast removal of residuals during optimization
    ceres::Problem graph(problem_options);
    ceres::Solver solver;
    ceres::Solver::Summary summary;
    ceres::Solver::Options options;
    
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.use_nonmonotonic_steps = true;
    options.initial_trust_region_radius = 1e6;
    options.dynamic_sparsity = true;
    options.num_threads = 8;
    options.max_num_iterations = 50;

    // Add velocity and clock frequency deviation parameters
    for (size_t k = 0; k < node_list.size(); k++) {
        graph.AddParameterBlock(&(node_list[k]->xf[3]), 3); 
        graph.AddParameterBlock(&(node_list[k]->clockshift), 1); 
    }

    // Add Doppler factors with Huber loss for Doppler-shift measurements
    double robust_control_dop = 0.10;
    for (size_t k = 0; k < node_list.size(); k++) {
        Add_doppler_factor(&graph, node_list[k], robust_control_dop);
    }

    // Add CVH factors for consecutive nodes
    if (node_list.size() >= 2) {
        for (size_t k = 0; k < node_list.size() - 1; k++) {
            if (node_list[k + 1]->nepoch - node_list[k]->nepoch > 1) continue;
            Add_cvh_factor(&graph, node_list[k], node_list[k + 1]);
        }
    }

    // Add CAH factors for triplets of consecutive nodes
    if (node_list.size() >= 3) {
        for (size_t k = 0; k < node_list.size() - 2; k++) {
            if ((node_list[k + 2]->nepoch - node_list[k + 1]->nepoch > 1) || (node_list[k + 1]->nepoch - node_list[k]->nepoch > 1))
                continue;
            Add_cah_factor(&graph, node_list[k], node_list[k + 1], node_list[k + 2]);
        }
    }

    // Add VFG marginalization prior factor
    if (vfg_last_marginalization_info && vfg_last_marginalization_info->isValid()) {
        auto factor = new MarginalizationFactor(vfg_last_marginalization_info);
        graph.AddResidualBlock(factor, nullptr, vfg_last_marginalization_parameter_blocks);
    }

    // Solve the optimization problem
    solver.Solve(options, &graph, &summary);
    std::cout << "VFG Optimization: " << summary.BriefReport() << "\n";

    return;
}


/** 
 * @brief Adds Doppler factors to the optimization problem.
 * 
 * @param graph The optimization problem.
 * @param node The node information.
 * @param robust_control_dop A robust control parameter for Doppler-shift measurements. If greater than zero, 
 *                            a Huber loss function will be applied.
 */
void Add_doppler_factor(ceres::Problem* graph, node_info* node, double robust_control_dop) {
    int i, f, sat;
    
    // Check if Doppler data should be used for this node
    if (!node->USE_Dopp) return;

    for (i = 0; i < node->nu; i++) {
        for (f = 0; f < NFREQ; f++) {
            sat = node->obs[i].sat;
            if (!node->Dopp_Valid[sat - 1][f]) continue;

            // If robust control is enabled, use a Huber loss function for the residual
            if (robust_control_dop > 0.0) {
                ceres::LossFunction* loss_function = new ceres::HuberLoss(robust_control_dop);
                auto doppfactor = new Dopp_Factor(node, i, f);
                graph->AddResidualBlock(doppfactor, loss_function, &(node->xf[3]), &(node->clockshift));
            } else {
                auto doppfactor = new Dopp_Factor(node, i, f);
                graph->AddResidualBlock(doppfactor, nullptr, &(node->xf[3]), &(node->clockshift));
            }
        }
    }
    return;
}

/** 
 * @brief Generates random combinations of prior combinations for RANSAC.
 * 
 * @param combinations The vector to store the selected random combinations.
 * @param combinations_prior The vector of prior combinations to sample from.
 */
void generate_random_combinations(vector<vector<int>> &combinations, vector<vector<int>> &combinations_prior) {
    // Set random number seed for reproducibility
    unsigned seed = 42;
    std::mt19937 generator(seed);

    int N = combinations_prior.size();
    std::uniform_int_distribution<int> distribution(0, N - 1);
    std::set<int> indexs;

    // Randomly select unique indices from combinations_prior
    while (indexs.size() < RANSAC_ITERATION) {
        int index = distribution(generator);
        indexs.insert(index);
    }

    // Add selected combinations to the output vector
    for (auto index : indexs) {
        combinations.push_back(combinations_prior.at(index));
    }

    return;
}


/** 
 * @brief Performs RANSAC-based Doppler shift quality control for a node.
 * 
 * @param node The node information.
 * @return Returns true if a sufficient set of valid Doppler-shift measurements is found, false otherwise.
 */
bool RANSAC_Doppler_shift_quality_control(node_info *node) {
    size_t j, max_value = 0;
    int i, f, dopflg, max_index = 0;
    vector<int> doppler_set;
    vector<vector<int>> Inner_Points_Sets;

    // Reset Doppler validity flags for all satellites and frequencies
    for (i = 0; i < MAXSAT; i++) {
        for (f = 0; f < NFREQ; f++) {
            node->Dopp_Valid[i][f] = false;
        }
    }

    // Collect valid Doppler-shift measurements
    for (i = 0; i < node->nu; i++) {
        for (f = 0; f < NFREQ; f++) {
            if (Is_Dop_Valid(node, i, f)) {
                dopflg = (i << 8) | f;
                doppler_set.push_back(dopflg);
            }
        }
    }

    node->nd = doppler_set.size();
    node->ndopp = node->nd;

    // If there are fewer valid Doppler-shift measurements than the minimum required, skip RANSAC
    if (doppler_set.size() <= MIN_SET) {
        node->USE_Dopp = false;
        node->ndopp = doppler_set.size();
        return false;
    }
    // If there are enough measurements for validation, mark them as valid without RANSAC
    else if (doppler_set.size() <= MIN_RANSAC) {
        node->USE_Dopp = true;
        node->ndopp = doppler_set.size();
        // Mark all Doppler-shift measurements as valid
        for (j = 0; j < doppler_set.size(); j++) {
            dopflg = doppler_set.at(j);
            node->Dopp_Valid[node->obs[(dopflg >> 8) & 0xFF].sat - 1][dopflg & 0xFF] = true;
        }
        return false;
    }

    // Generate random combinations of Doppler-shift measurements for RANSAC iterations
    vector<vector<int>> combinations_prior;
    CombinationGenerator generator(node->ndopp, MIN_SET);
    combinations_prior = generator.getCombinations();

    vector<vector<int>> combinations;
    generate_random_combinations(combinations, combinations_prior);

    // Perform RANSAC iterations to find valid inlier sets
    for (auto combination : combinations) {
        vector<int> Inner_Points;
        Dopp_RANSAC_Single_Iterition(node, doppler_set, combination, Inner_Points);
        Inner_Points_Sets.push_back(Inner_Points);
    }

    // Identify the largest inlier set from all RANSAC iterations
    for (j = 0; j < Inner_Points_Sets.size(); j++) {
        if (Inner_Points_Sets.at(j).size() > max_value) {
            max_value = Inner_Points_Sets.at(j).size();
            max_index = j;
        }
    }

    // If RANSAC doesn't find a sufficient number of inliers, mark them as valid anyway
    node->ndopp = Inner_Points_Sets.at(max_index).size();

    if (node->ndopp <= MIN_SET) {
        node->USE_Dopp = true;
        for (j = 0; j < Inner_Points_Sets.at(max_index).size(); j++) {
            dopflg = doppler_set.at(Inner_Points_Sets.at(max_index).at(j));
            node->Dopp_Valid[node->obs[(dopflg >> 8) & 0xFF].sat - 1][dopflg & 0xFF] = true;
        }
        return false;
    }

    // Mark Doppler-shift measurements in the largest inlier set as valid
    for (j = 0; j < Inner_Points_Sets.at(max_index).size(); j++) {
        dopflg = doppler_set.at(Inner_Points_Sets.at(max_index).at(j));
        node->Dopp_Valid[node->obs[(dopflg >> 8) & 0xFF].sat - 1][dopflg & 0xFF] = true;
    }

    return true;
}


/** 
 * @brief Performs a single RANSAC iteration to find inliers based on Doppler-shift measurements.
 * 
 * @param node The node information.
 * @param doppler_set A list of all valid Doppler-shift measurement indices.
 * @param Minset_index The indices of the smallest subset of Doppler-shift measurements used for estimation.
 * @param Inner_Points A vector that will store the indices of the inliers.
 */
void Dopp_RANSAC_Single_Iterition(node_info *node, vector<int> &doppler_set, vector<int> &Minset_index, vector<int> &Inner_Points) {
    double x[4] = { 0 }, Q[16] = { 0 };
    double *v, *H;
    int dopn = doppler_set.size(), i;
    vector<int> all(dopn);
    iota(all.begin(), all.end(), 0);  // Initialize all indices from 0 to dopn-1

    v = mat(dopn, 1);  // Allocate memory for residual vector
    H = mat(4, dopn);  // Allocate memory for measurement matrix

    // Estimate model parameters using the given minimal set of measurements
    Minset_Estimate(node, doppler_set, Minset_index, x, Q);

    // Compute residuals for all Doppler-shift measurements
    resdop_multifreq(node, doppler_set, all, v, H, x);

    // Evaluate which measurements are inliers by checking if their residuals are below the threshold
    for (i = 0; i < dopn; i++) {
        if (fabs(v[i]) < DOP_THRE) {
            Inner_Points.push_back(i);  // Store inlier indices
        }
    }

    free(v);  
    free(H); 

    return;
}


/** 
 * @brief Performs least squares estimation for the minimal set of Doppler-shift measurements.
 * 
 * @param node The node information.
 * @param doppler_set A list of all valid Doppler-shift measurement indices.
 * @param Minset_index The indices of the smallest subset of Doppler-shift measurements used for estimation.
 * @param x The estimated model parameters.
 * @param Q The covariance matrix of the estimated parameters.
 */
void Minset_Estimate(node_info *node, vector<int> &doppler_set, vector<int> &Minset_index, double *x, double *Q) {
    int i, j;
    double *v, *H, dx[4] = { 0 };

    // Allocate memory for residual vector and measurement matrix
    v = mat(MIN_SET, 1);  
    H = mat(4, MIN_SET);

    // Initialize the model parameters
    for (j = 0; j < 3; j++) x[j] = node->xf[j + 3];

    for (i = 0; i < MAXITR; i++) {
        // Compute the residuals and the measurement matrix for the current parameters
        resdop_multifreq(node, doppler_set, Minset_index, v, H, x);

        // Perform least squares estimation to obtain the parameters
        if (lsq(H, v, 4, MIN_SET, dx, Q)) break;

        for (j = 0; j < 4; j++) x[j] += dx[j];

        if (norm(dx, 4) < 1E-6) {
            break;
        }
    }

    free(v);  
    free(H);
    
    return;
}


/** 
 * @brief Computes residuals and the Jacobian for the Doppler-shift measurements.
 * 
 * @param node The node information.
 * @param doppler_set A list of all valid Doppler-shift measurement indices.
 * @param Index_set The indices of Doppler-shift measurements used in the current iteration.
 * @param v The vector to store the computed residuals.
 * @param H The matrix to store the computed Jacobian.
 * @param x The current estimate of the model parameters.
 */
void resdop_multifreq(node_info *node, vector<int> &doppler_set, vector<int> &Index_set, double *v, double *H, double *x) {
    size_t i;
    int j, satindex, f, dopflg;
    double vs[3], e[3], freq, rate;

    for (i = 0; i < Index_set.size(); i++) {
        dopflg = doppler_set.at(Index_set.at(i));  // Get the Doppler-shift flag
        satindex = (dopflg >> 8) & 0xFF;
        f = dopflg & 0xFF; 
        
        // Get the frequency
        freq = sat2freq(node->obs[satindex].sat, node->obs[satindex].code[f], &navs);


		// Compute the measurement residual and Jacobian 
        for (j = 0; j < 3; j++) {
            vs[j] = node->rs[j + 3 + satindex * 6] - x[j];
            e[j] = node->e[j + satindex * 3];
        }

        rate = dot(vs, e, 3) + OMGE / CLIGHT * (node->rs[4 + satindex * 6] * node->xf[0] + node->rs[1 + satindex * 6] * x[0] -
            node->rs[3 + satindex * 6] * node->xf[1] - node->rs[satindex * 6] * x[1]);

        v[i] = -node->obs[satindex].D[f] * CLIGHT / freq - (rate + x[3] - CLIGHT * node->dts[1 + satindex * 2]);

        // Fill in the Jacobian matrix (H)
        for (j = 0; j < 4; j++) {
            H[j + i * 4] = (j < 3) ? -e[j] : 1.0;
        }
    }

    return;
}


/** 
 * @brief Checks if a Doppler-shift measurement is valid.
 * 
 * @param node The node information.
 * @param index Index of the Doppler-shift measurement in the node's observations.
 * @param f Frequency index.
 * 
 * @return bool True if the Doppler-shift measurement is valid, false otherwise.
 */
bool Is_Dop_Valid(node_info *node, int index, int f) {
    double freq;
    
    freq = sat2freq(node->obs[index].sat, node->obs[index].code[f], &navs);
    
    if (node->obs[index].D[f] == 0.0 || freq == 0.0 || 
        !node->ssat[node->obs[index].sat - 1].vsat[f] || 
        satexclude(node->obs[index].sat, node->var[index], node->svh[index], &(option)) || 
        norm(node->rs + 3 + index * 6, 3) <= 0.0) {
        return false; // Invalid Doppler-shift
    }
    
    return true; // Valid Doppler-shift
}




/** 
 * @brief Adds a CVH factor to the optimization problem.
 * 
 * @param graph The optimization problem.
 * @param node0 The first node in the sequence.
 * @param node1 The second node in the sequence.
 */
void Add_cvh_factor(ceres::Problem* graph, node_info* node0, node_info* node1) {
    auto cv_factor = new CVH_Factor(VEL_CON_VAR_CV);  // Create CVH factor with specified variance
    graph->AddResidualBlock(cv_factor, nullptr, &(node0->xf[3]), &(node1->xf[3]));  // Add residual block
}


/** 
 * @brief Adds a CAH factor to the optimization problem.
 * 
 * @param graph The optimization problem.
 * @param node0 The first node in the sequence.
 * @param node1 The second node in the sequence.
 * @param node2 The third node in the sequence.
 */
void Add_cah_factor(ceres::Problem* graph, node_info* node0, node_info* node1, node_info* node2) {
    auto ca_factor = new CAH_Factors(VEL_CON_VAR_CA);  // Create CAH factor with specified variance
    graph->AddResidualBlock(ca_factor, nullptr, &(node0->xf[3]), &(node1->xf[3]), &(node2->xf[3]));  // Add residual block
}


/** 
 * @brief Performs the VFG marginalization process.
 * 
 */
void VFG_marginalization() {
    
	// Check for interruption conditions
    bool uncon01, uncon12;
    uncon01 = (node_list[1]->nepoch - node_list[0]->nepoch) > 1;
    uncon12 = (node_list[2]->nepoch - node_list[1]->nepoch) > 1;
    
    if (uncon01 || vfg_last_marginalization_parameter_blocks.empty()) {
        vfg_last_marginalization_parameter_blocks.clear();
        vfg_last_marginalization_info = NULL;
        return;
    }

    std::shared_ptr<MarginalizationInfo> marginalization_info = std::make_shared<MarginalizationInfo>();

    if (vfg_last_marginalization_info && vfg_last_marginalization_info->isValid()) {
        std::vector<int> marginilized_index;
        for (size_t k = 0; k < vfg_last_marginalization_parameter_blocks.size(); k++) {
            if (vfg_last_marginalization_parameter_blocks[k] == &(node_list[0]->xf[3])) {
                marginilized_index.push_back(static_cast<int>(k));
            }
        }
        auto factor = std::make_shared<MarginalizationFactor>(vfg_last_marginalization_info);
        auto residual = std::make_shared<ResidualBlockInfo>(
            factor, nullptr, vfg_last_marginalization_parameter_blocks, marginilized_index);
        marginalization_info->addResidualBlockInfo(residual);
    }

    // Add factors
    if (node_list[0]->USE_Dopp) {
        int i, f, sat;
        for (i = 0; i < node_list[0]->nu; i++) {
            for (f = 0; f < NFREQ; f++) {
                sat = node_list[0]->obs[i].sat;
                if (!node_list[0]->Dopp_Valid[sat - 1][f]) continue;
                auto factor = std::make_shared<Dopp_Factor>(node_list[0], i, f);
                auto residual = std::make_shared<ResidualBlockInfo>(
                    factor, nullptr,
                    std::vector<double*>{&(node_list[0]->xf[3]), &(node_list[0]->clockshift)}, std::vector<int>{0, 1});
                marginalization_info->addResidualBlockInfo(residual);
            }
        }
    }

    if (!uncon01) {
        auto factor = std::make_shared<CVH_Factor>(VEL_CON_VAR_CV);
        auto residual = std::make_shared<ResidualBlockInfo>(
            factor, nullptr,
            std::vector<double*>{&(node_list[0]->xf[3]), &(node_list[1]->xf[3])}, std::vector<int>{0});
        marginalization_info->addResidualBlockInfo(residual);
    }

    if ((!uncon01) && (!uncon12)) {
        auto factor = std::make_shared<CAH_Factors>(VEL_CON_VAR_CA);
        auto residual = std::make_shared<ResidualBlockInfo>(
            factor, nullptr,
            std::vector<double*>{&(node_list[0]->xf[3]), &(node_list[1]->xf[3]), &(node_list[2]->xf[3])}, std::vector<int>{0});
        marginalization_info->addResidualBlockInfo(residual);
    }

    // Do marginalization
    marginalization_info->marginalization();

    // Get the new parameter block pointers
    std::unordered_map<long, double*> address;
    for (size_t k = 1; k < node_list.size(); k++) {
        address[reinterpret_cast<long>(&(node_list[k]->xf[3]))] = &(node_list[k]->xf[3]);
    }

	// Update the marginalization information
    vfg_last_marginalization_parameter_blocks = marginalization_info->getParamterBlocks(address);
    vfg_last_marginalization_info = std::move(marginalization_info);

    return;
}
