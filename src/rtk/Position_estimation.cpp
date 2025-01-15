#include "../global/global.h"
#include "../global/marginalization_info.h"
#include "../global/marginalization_factor.h"
#include "../factors/dd_l_factor.h"
#include "../factors/dd_p_factor.h"
#include "../factors/ambi_prior_factor.h"
#include "../factors/fix_pos_prior_factor.h"
#include "../factors/cvh_factor.h"
#include "../factors/vic_factor.h"
#include "../ar/artest.h"
#include "../rtklib/rtklib.h"
#include <ceres/normal_prior.h>

/** 
 * @brief Performs the PFG optimization.
 * 
 * @return Returns true if the optimization converges and returns false if the optimization falied.
 */
bool PFG_optimization() {

    // Solver configuration
    ceres::Problem::Options problem_options;
    problem_options.enable_fast_removal = true;
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

    // Add parameters and residuals
    std::vector<ceres::ResidualBlockId> BackendId;
    Add_Parameter_and_Residual(&graph, BackendId);

    // Solve the problem
    solver.Solve(options, &graph, &summary);
    std::cout << "PFG Optimization: " << summary.BriefReport() << "\n";

    // Check solution validity
    if (!summary.IsSolutionUsable()) {
        trace(3, "Check optimization failure!\n");
		return false;
    }

    // Perform ambiguity resolution
    Ambiguity_resolution(&graph, BackendId);

    return true;
}

/** 
 * @brief Performs ambiguity resolution for the current node.
 * 
 * @param graph The ceres problem object.
 * @param BackendId A list of residual block IDs that are restored during constructing the PFG.
 */
void Ambiguity_resolution(ceres::Problem* graph, std::vector<ceres::ResidualBlockId> &BackendId) {
	// Initialize node status and reset the fix count
	node_list.back()->stat = SOLQ_FLOAT;
	node_list.back()->nfix = 0;

	if (node_list.back()->nb >= MINNB) {
		node_info *node = node_list.back();

		// Clear previous parameters and residuals, then add new ones
		Remove_All_Parameters_and_Residuals(graph);
		BackendId.clear();
		Add_dd_factor(graph, node, 0.0, 0.0, BackendId);

		// Attempt ambiguity resolution and update node status
		int nb = resamb_partial(graph, node);
		if (nb >= MINNB) {
			node->stat = SOLQ_FIX;
			fixcnt++;
		} else {
			node->stat = SOLQ_FLOAT;
			fixcnt = 0;
		}

		// Update node's fixed count and number of integer ambiguities
		node->nfix = fixcnt;
		node->nb = nb;
	}
}

/** 
 * @brief Adds parameters and residuals to the ceres optimization problem.
 * 
 * @param graph The ceres problem object.
 * @param BackendId List of residual block IDs to be updated.
 */
void Add_Parameter_and_Residual(ceres::Problem* graph, std::vector<ceres::ResidualBlockId> &BackendId) {
	
	// Add parameter blocks
	for (size_t k = 0; k < node_list.size(); k++) {
		graph->AddParameterBlock(&(node_list[k]->xf[0]), 3);
	}
	for (size_t k = 0; k < sdambi_info_list.size(); k++) {
		if (!sdambi_info_list[k]->Ismarg && sdambi_info_list[k]->num > 0) {
			graph->AddParameterBlock(&(sdambi_info_list[k]->sdbias), 1);
		}
	}

	for (size_t k = 0; k < node_list.size(); k++) {
		graph->AddParameterBlock(&(node_list[k]->xf[3]), 3);
		graph->AddParameterBlock(&(node_list[k]->clockshift), 1);
	}

	// Add residual blocks for DD pseudorange and DD carrier phase factors
	double robust_control_ddp = 1.0;
	double robust_control_ddl = 0.5;
	for (size_t k = 0; k < node_list.size(); k++) {
		Add_dd_factor(graph, node_list[k], robust_control_ddp, robust_control_ddl, BackendId);
	}

	// Add residual blocks for VIC factors between consecutive nodes
	double robust_control_vic = 0.0;
	if (node_list.size() >= 10) {
		for (size_t k = 0; k < node_list.size() - 1; k++) {
			if (node_list[k + 1]->nepoch - node_list[k]->nepoch > 1) continue;
			if (!have_first_fix) continue;
			Add_vic_factor(graph, node_list[k], node_list[k + 1], robust_control_vic);
		}
	}

	// Add prior factors for fixed ambiguities, position and float ambiguities
	set<int> priorset;
	for (size_t k = 0; k < node_list.size(); k++) {
		if (node_list[k]->stat == SOLQ_FIX) {
			Add_fix_ambi_prior(graph, node_list[k], priorset);
			Add_fix_pos_prior(graph, node_list[k]);
		}
		Add_float_ambi_prior(graph, node_list[k], priorset);
	}

	// Add PFG marginalization prior factor
	Add_pfg_marginalization_prior_factor(graph);
}


/** 
 * @brief Removes all parameters and residuals from the ceres optimization problem.
 * 
 * @param graph The ceres problem object.
 */
void Remove_All_Parameters_and_Residuals(ceres::Problem* graph) {
	// Remove all residual blocks from the problem
	std::vector<ceres::ResidualBlockId> residual_block_ids;
	graph->GetResidualBlocks(&residual_block_ids);
	for (ceres::ResidualBlockId residual_block_id : residual_block_ids) {
		graph->RemoveResidualBlock(residual_block_id);
	}

	// Remove all parameter blocks from the problem
	std::vector<double*> parameter_blocks;
	graph->GetParameterBlocks(&parameter_blocks);
	for (double* parameter_block : parameter_blocks) {
		graph->RemoveParameterBlock(parameter_block);
	}
}


/** 
 * @brief Adds a fix position prior factor to the ceres optimization problem.
 * 
 * @param graph The ceres problem object.
 * @param node The node information.
 */
void Add_fix_pos_prior(ceres::Problem* graph, node_info* node) {
	auto pos_prior_factor = new POS_Fix_Prior(node, FIX_POS_PRIOR_VAR);
	graph->AddResidualBlock(pos_prior_factor, nullptr, &(node->xf[0]));
}


/** 
 * @brief Adds a VIC factor to the ceres optimization problem.
 * 
 * @param graph The ceres problem instance.
 * @param node_pre The previous node in the sequence.
 * @param node_now The current node in the sequence.
 * @param robust_control_vic The control parameter for applying a robust loss function. 
 * If greater than zero, Huber loss is applied.
 */
void Add_vic_factor(ceres::Problem* graph, node_info* node_pre, node_info* node_now, double robust_control_vic) {
    // If robust control is enabled, use Huber loss to handle outliers
	if (robust_control_vic > 0) {
        ceres::LossFunction* loss_function = new ceres::HuberLoss(robust_control_vic);
        auto tufactor = new VIC_Factor(node_pre, node_now);
        // Add the residual block with the loss function
        graph->AddResidualBlock(tufactor, loss_function, &(node_pre->xf[0]), &(node_now->xf[0]));
    }
    else {
        auto tufactor = new VIC_Factor(node_pre, node_now);
        graph->AddResidualBlock(tufactor, nullptr, &(node_pre->xf[0]), &(node_now->xf[0]));
    }
}


/** 
 * @brief Adds a marginalization prior factor to the ceres optimization problem.
 * 
 * @param graph The ceres problem instance.
 */
void Add_pfg_marginalization_prior_factor(ceres::Problem* graph) {
    // Check if the last marginalization info is valid before adding the factor
    if (pfg_last_marginalization_info && pfg_last_marginalization_info->isValid()) {
        // Create the marginalization factor
        auto factor = new MarginalizationFactor(pfg_last_marginalization_info);
        graph->AddResidualBlock(factor, nullptr, pfg_last_marginalization_parameter_blocks);
    }
}


/** 
 * @brief Adds float ambiguity prior factors to the ceres optimization problem.
 * 
 * @param graph The ceres problem instance.
 * @param node The node information.
 * @param priorset A set to track added ambiguity pairs and prevent duplication.
 */
void Add_float_ambi_prior(ceres::Problem* graph, node_info *node, set<int> &priorset) {
    int i = 0;
    int refindex, urefindex, vflg;

    for (i = 0; i < node->nv; i++) {

        // Skip invalid or pseudorange measurements
        if (!node->ddmeas[i].isValid || node->ddmeas[i].code) continue;

        // Only add priors for ambiguity measurements marked as "prior" (have little residual)
        if (!node->ddmeas[i].isPrior) continue;

        refindex = node->ddmeas[i].ref_sdambi_index;
        urefindex = node->ddmeas[i].uref_sdambi_index;

        // Create a unique flag for this reference-uref pair to avoid duplication
        vflg = (refindex << 16) | (urefindex << 8);
        if (priorset.find(vflg) != priorset.end()) {
            continue; // Skip if this pair has already been added
        }

        // Create the prior factor and add it to the optimization problem
        double var = FLOAT_PRIOR_VAR;
        auto prior = new AMBI_Prior(node->ddmeas[i].ddbias, var);
        graph->AddResidualBlock(prior, nullptr, &(sdambi_info_list.at(node->ddmeas[i].ref_sdambi_index)->sdbias),
            &(sdambi_info_list.at(node->ddmeas[i].uref_sdambi_index)->sdbias));

        // Add the ambiguity pair to the set to avoid duplicates
        priorset.insert(vflg);
    }
}


/** 
 * @brief Adds fixed ambiguity prior factors to the ceres optimization problem.
 * 
 * @param graph The ceres problem instance.
 * @param node The node information.
 * @param priorset A set to track added ambiguity pairs and prevent duplication.
 */
void Add_fix_ambi_prior(ceres::Problem* graph, node_info *node, set<int> &priorset) {
    int i = 0;
    int refindex, urefindex, vflg;

    // Only add prior factors if the node has fixed ambiguities
    if (node->nb <= 0 || node->stat != SOLQ_FIX) {
        return;
    }

    for (i = 0; i < node->nb; i++) {
        refindex = node->ddfix[i].ref_sdambi_index;
        urefindex = node->ddfix[i].uref_sdambi_index;

        // Create a unique flag for this reference-uref pair to avoid duplication
        vflg = (refindex << 16) | (urefindex << 8);
        if (priorset.find(vflg) != priorset.end()) {
            continue; // Skip if this pair has already been added
        }

        // Create the prior factor and add it to the optimization problem
        double var = FIX_AMBI_PRIOR_VAR;
        auto prior = new AMBI_Prior(node->ddfix[i].ddbias, var);
        graph->AddResidualBlock(prior, nullptr, &(sdambi_info_list.at(node->ddfix[i].ref_sdambi_index)->sdbias),
            &(sdambi_info_list.at(node->ddfix[i].uref_sdambi_index)->sdbias));

        // Add the ambiguity pair to the set to avoid duplicates
        priorset.insert(vflg);
    }
}


/** 
 * @brief Adds DD pseudorange and DD carrier phase factors to the ceres optimization problem.
 * 
 * @param graph The ceres problem instance.
 * @param node The node information.
 * @param robust_control_ddp The robust control parameter for DD_P factors (Huber loss function).
 * @param robust_control_ddl The robust control parameter for DD_L factors (Huber loss function).
 * @param BackendId A vector to store the added residual block IDs.
 */
void Add_dd_factor(ceres::Problem* graph, node_info* node, double robust_control_ddp, double robust_control_ddl, std::vector<ceres::ResidualBlockId>& BackendId) {
    
    int i = 0;
    for (i = 0; i < node->nv; i++) {
        if (!node->ddmeas[i].isValid) continue; // Skip invalid measurements
        if (node->ddmeas[i].code) {
            if (robust_control_ddp > 0.0) {
                // Apply robust control (Huber loss) for DD_P factor
                ceres::LossFunction* loss_function = new ceres::HuberLoss(robust_control_ddp);
                auto ddpfactor = new DD_P_Factor(node, i);
                ceres::ResidualBlockId residual_block_id = graph->AddResidualBlock(ddpfactor, loss_function, &(node->xf[0]));
                BackendId.push_back(residual_block_id); // Store the residual block ID
            }
            else {
                auto ddpfactor = new DD_P_Factor(node, i);
                ceres::ResidualBlockId residual_block_id = graph->AddResidualBlock(ddpfactor, nullptr, &(node->xf[0]));
                BackendId.push_back(residual_block_id); // Store the residual block ID
            }
        }
        else {
            if (robust_control_ddl > 0.0) {
                // Apply robust control (Huber loss) for DD_L factor
                ceres::LossFunction* loss_function = new ceres::HuberLoss(robust_control_ddl);
                auto ddlfactor = new DD_L_Factor(node, i);
                ceres::ResidualBlockId residual_block_id = graph->AddResidualBlock(ddlfactor, loss_function, &(node->xf[0]), &(sdambi_info_list.at(node->ddmeas[i].ref_sdambi_index)->sdbias),
                    &(sdambi_info_list.at(node->ddmeas[i].uref_sdambi_index)->sdbias));
                BackendId.push_back(residual_block_id); // Store the residual block ID
            }
            else {
                auto ddlfactor = new DD_L_Factor(node, i);
                ceres::ResidualBlockId residual_block_id = graph->AddResidualBlock(ddlfactor, nullptr, &(node->xf[0]), &(sdambi_info_list.at(node->ddmeas[i].ref_sdambi_index)->sdbias),
                    &(sdambi_info_list.at(node->ddmeas[i].uref_sdambi_index)->sdbias));
                BackendId.push_back(residual_block_id); // Store the residual block ID
            }
        }
    }
    
    return;
}


/** 
 * @brief Performs the PFG marginalization process.
 * 
 */
void PFG_marginalization() {

    bool uncon01;
    uncon01 = (node_list[1]->nepoch - node_list[0]->nepoch) > 1; // Check if the epochs are uncontinuous

    if (uncon01) {
        pfg_last_marginalization_parameter_blocks.clear();
        pfg_last_marginalization_info = NULL;
        return;
    }
    
    std::shared_ptr<MarginalizationInfo> marginalization_info = std::make_shared<MarginalizationInfo>();

    if (pfg_last_marginalization_info && pfg_last_marginalization_info->isValid()) {
        std::vector<int> marginilized_index;
        
        for (size_t k = 0; k < pfg_last_marginalization_parameter_blocks.size(); k++) {
            if (pfg_last_marginalization_parameter_blocks[k] == &(node_list[0]->xf[0])) {
                marginilized_index.push_back(static_cast<int>(k));
            }
        }

        auto factor = std::make_shared<MarginalizationFactor>(pfg_last_marginalization_info);
        auto residual = std::make_shared<ResidualBlockInfo>(
            factor, nullptr, pfg_last_marginalization_parameter_blocks, marginilized_index);
        marginalization_info->addResidualBlockInfo(residual); 
    }

	// Add factors
    if (!uncon01) {
        auto factor = std::make_shared<VIC_Factor>(node_list[0], node_list[1]);
        auto residual = std::make_shared<ResidualBlockInfo>(
            factor, nullptr,
            std::vector<double*>{&(node_list[0]->xf[0]), &(node_list[1]->xf[0])}, std::vector<int>{0});
        marginalization_info->addResidualBlockInfo(residual); // Add to marginalization info
    }

    if (node_list[0]->stat == SOLQ_FIX) {
        auto factor = std::make_shared<POS_Fix_Prior>(node_list[0], FIX_POS_PRIOR_VAR);
        auto residual = std::make_shared<ResidualBlockInfo>(
            factor, nullptr,
            std::vector<double*>{&(node_list[0]->xf[0])}, std::vector<int>{0});
        marginalization_info->addResidualBlockInfo(residual); // Add to marginalization info
    }

    // Do marginalization
    marginalization_info->marginalization();

    // Get the new parameter block pointers
    std::unordered_map<long, double*> address;
    for (size_t k = 1; k < node_list.size(); k++) {
        address[reinterpret_cast<long>(&(node_list[k]->xf[0]))] = &(node_list[k]->xf[0]);
    }

    // Update the marginalization information
    pfg_last_marginalization_parameter_blocks = marginalization_info->getParamterBlocks(address);
    pfg_last_marginalization_info = std::move(marginalization_info);

    return;
}

