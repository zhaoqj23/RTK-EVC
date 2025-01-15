#include "global.h"

deque<node_info*> node_list;
deque<sdambi_info*> sdambi_info_list;

size_t window_size = 100;
int epoch_now = 0;
int head_flg = 1;
int fixcnt = 0;
bool have_first_fix = false;

// TODO: Paths need to be modified!
string rover_path = "/home/zhaoqj23/Documents/RTK-EVC/dataset/0928/Rover/cpt.obs";
string base_path = "/home/zhaoqj23/Documents/RTK-EVC/dataset/0928/Base/0928basertcm.obs";
string nav_path = "/home/zhaoqj23/Documents/RTK-EVC/dataset/0928/Base/0928basertcm.nav";
string config_path = "/home/zhaoqj23/Documents/RTK-EVC/dataset/0928/rtkevc.conf";
string rtklib_path = "/home/zhaoqj23/Documents/RTK-EVC/dataset/0928/Rover/RTKLIB.pos";
string rtkevc_path = "/home/zhaoqj23/Documents/RTK-EVC/dataset/0928/Rover/RTK-EVC.pos";

// Shared pointers for marginalization information and parameter blocks
std::shared_ptr<MarginalizationInfo> pfg_last_marginalization_info;
std::vector<double*> pfg_last_marginalization_parameter_blocks;

std::shared_ptr<MarginalizationInfo> vfg_last_marginalization_info;
std::vector<double*> vfg_last_marginalization_parameter_blocks;

/**
 * @brief Outputs data to the RTK-EVC file, calling header output if it's the first data.
 * 
 * @param node Pointer to the node information.
 * @param rtk Pointer to the RTK data structure.
 */
void outputdata(node_info *node, rtk_t *rtk) {
    if (head_flg) {
        outputhead();
        head_flg = 0;
    }
    output_xfdata(node);
}

/**
 * @brief Outputs the node's position and velocity data to the RTK-EVC file.
 * 
 * @param node Pointer to the node information.
 */
void output_xfdata(node_info *node) {
    int week;
    double second, blh[3], vec[3];

    // Convert from ECEF to latitude, longitude, height (blh) and ENU (East-North-Up) coordinates
    ecef2pos(node->xf, blh);
    ecef2enu(blh, &(node->xf[3]), vec);
    
    second = time2gpst(node->time, &week);

    // Open RTK-EVC file for appending
    FILE* f_ceres = fopen(rtkevc_path.c_str(), "a+");

    // Define the radian to degree conversion factor
    double rad2deg = 57.29577951308232;

    // Write the data into the file in a tab-delimited format
    fprintf(f_ceres, "%d\t%f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\t%.5f\t%d\n", 
            week, second, 
            blh[0] * rad2deg, blh[1] * rad2deg, blh[2], 
            vec[0], vec[1], vec[2], 
            node->stat, node->ratio, node->nfix);

    // Close the file
    fclose(f_ceres);
}

/**
 * @brief Outputs the header for the RTK-EVC data file.
 */
void outputhead() {
    // Open RTK-EVC file for writing the header
    FILE* f_ceres = fopen(rtkevc_path.c_str(), "w");
    
    // Write the header line
    fprintf(f_ceres, "\tGPST\tlatitude(deg)\tlongitude(deg)\theight(m)\tVe(m/s)\tVn(m/s)\tVu(m/s)\tQ\tRatio\tnfix\n");
    
    // Close the file
    fclose(f_ceres);
}
