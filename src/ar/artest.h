#pragma once
#include "../global/global.h"

#define MINNB 8
#ifdef HIGH_SAMPLE
    #define CONP 50
    #define CON_THRE 0.5
#else
    #define CONP 5
    #define CON_THRE 1.5
#endif


extern int resamb_partial(ceres::Problem* graph, node_info *node);
extern void fillup_amb_weight_flag(node_info *node);
extern int ar_lambda(ceres::Problem *graph, node_info *node, double *bias, int *ix, int nfix, MatrixXd &Qc, MatrixXd &Qac, vector<int> sdambi_index_inorder);
extern int select_ddidx(node_info *node, set<int> &ambi_index_tofix, int *ix);
extern bool covariance_estimation(ceres::Problem* graph, MatrixXd &Qc, MatrixXd &Qac, std::vector<int> &sdambi_index_inorder, node_info* node);
extern void sdambi_restore(int nfix, int *ix, double *bias, node_info *node);
extern void find_ambiindex_tofix(set<int> &ambi_index_tofix, node_info *node);
extern int sdambi_inorder(int *ix, int nfix, vector<int> &sdambi_index_inorder);
extern int find_index(vector<int> &sdambi_index_inorder_origin, int sdambi_index);
extern bool consistency_check(MatrixXd &xaa);
extern void Dead_Reckoning(int index, double *xf);
extern void Ambiguity_resolution(ceres::Problem* graph, std::vector<ceres::ResidualBlockId> &BackendId);

