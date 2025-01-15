#pragma once
#include "../rtklib/rtklib.h"
#include "residual_block_info.h"
#include "marginalization_info.h"
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <ctime>
#include <algorithm>
#include <deque>
#include <eigen3/Eigen/Core>
#include <ceres/ceres.h>
#include <fstream>

using namespace std;
using namespace Eigen;


// Variance of prior factors
#define FIX_AMBI_PRIOR_VAR 0.001
#define FIX_POS_PRIOR_VAR 0.00001
#define FLOAT_PRIOR_VAR 1

// RANSAC parameters
#define DOP_THRE 0.5
#define RANSAC_ITERATION 108
#define MIN_SET 4
#define MIN_RANSAC 10

// Different parameters for sample rates 10 Hz/1 Hz
// #define HIGH_SAMPLE
#ifdef HIGH_SAMPLE
	#define VEL_CON_VAR_CV 1
	#define VEL_CON_VAR_CA 1
#else
	#define VEL_CON_VAR_CV 5
	#define VEL_CON_VAR_CA 5
#endif

// Global variables
extern string rover_path, base_path, nav_path, rtklib_path, rtkevc_path, config_path;

extern deque<node_info*> node_list;
extern deque<sdambi_info*> sdambi_info_list;

extern prcopt_t option;
extern nav_t navs;

extern size_t window_size;
extern int epoch_now;
extern int fixcnt;
extern int head_flg;
extern bool have_first_fix;

extern std::shared_ptr<MarginalizationInfo> pfg_last_marginalization_info;
extern std::vector<double*> pfg_last_marginalization_parameter_blocks;

extern std::shared_ptr<MarginalizationInfo> vfg_last_marginalization_info;
extern std::vector<double*> vfg_last_marginalization_parameter_blocks;


// Global
extern void outputdata(node_info* node, rtk_t *rtk);
extern void output_xfdata(node_info* node);
extern void outputhead();


// Velocity estimation: RANSAC Doppler-shift quality control
extern bool RANSAC_Doppler_shift_quality_control(node_info *node);
extern void Dopp_RANSAC_Single_Iterition(node_info *node, vector<int> &doppler_set, vector<int> &Minset_index, vector<int> &Inner_Points);
extern bool Is_Dop_Valid(node_info *node, int index, int f);
extern void Minset_Estimate(node_info *node, vector<int> &doppler_set, vector<int> &Minset_index, double *x, double *Q);
extern void resdop_multifreq(node_info *node, vector<int> &doppler_set, vector<int> &Index_set, double *v, double *H, double *x);

// Velocity estimation: construct and solve VFG
extern void Add_doppler_factor(ceres::Problem* graph, node_info* node, double robust_control_dop);
extern void Add_cvh_factor(ceres::Problem* graph, node_info* node0, node_info* node1);
extern void Add_cah_factor(ceres::Problem* graph, node_info* node0, node_info* node1, node_info* node2);
extern void VFG_optimization();
extern void VFG_marginalization();


// Position estimation: GSSM modeling
extern void ddmeas_init(rtk_t *rtk);
extern void ddfix_init(rtk_t *rtk);
extern int save_node_data(rtk_t *rtk, const int stat, const prcopt_t *popt, const nav_t *nav, const obsd_t *obs, const int nobs,
							double *rs, double *dts, double *var, int *svh, double *y, double *R, double* azel, double* e);
extern void rtkinformation(node_info* node, rtk_t* rtk, std::set<int> &fixflag);
extern void init_sdambi_infor(rtk_t *rtk, double sdbias, double var, int epoch_s, int sdflg, int stat);
extern void attach_new_sdambi(node_info* node, rtk_t* rtk, int sdflg, int stat);
extern void sdambi_merge(rtk_t *rtk, sdambi_info *sdambi, int stat);
extern void sdambi_gssm(map<int, int> *snmap, node_info *node, rtk_t *rtk, int sdflg, int stat);


// Position estimation: construct and solve PFG
extern void Add_float_ambi_prior(ceres::Problem* graph, node_info* node, set<int> &priorset);
extern void Add_fix_ambi_prior(ceres::Problem* graph, node_info *node, set<int> &priorset);
extern void Add_fix_pos_prior(ceres::Problem* graph, node_info* node);
extern void Add_vic_factor(ceres::Problem* graph, node_info* node_pre, node_info* node_now, double robust_control_vic);
extern void Add_pfg_marginalization_prior_factor(ceres::Problem* graph);
extern void Add_dd_factor(ceres::Problem* graph, node_info* node, double robust_control_ddp, double robust_control_ddl, std::vector<ceres::ResidualBlockId>& BackendId);
extern void Add_Parameter_and_Residual(ceres::Problem* graph, std::vector<ceres::ResidualBlockId> &BackendId);
extern bool PFG_optimization();
extern void Remove_All_Parameters_and_Residuals(ceres::Problem* graph);
extern int zdres(int base, const obsd_t* obs, int n, const double* rs, const double* dts, const double* var, const int* svh, const nav_t* nav, const double* rr, const prcopt_t* opt, int index, double* y, double* e, double* azel, double* freq);

// Position estimation: Sliding window and marginalization
extern void index_adjust(int index);
extern void PFG_marginalization();
extern void sliding_window();
extern void node_free(node_info *node);







