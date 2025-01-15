#include <stdarg.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "./rtklib/rtklib.h"
#include "./global/global.h"

prcopt_t option;

int main(int argc, char **argv)
{
    // Initialize Google Logging system
    google::InitGoogleLogging("RTK-EVC");
    
    // Open trace file for logging (Path needs to be modified!)
    traceopen("/home/zhaoqj23/Documents/RTK-EVC/Debug/rtkevc.trace");
    tracelevel(1);

    // Set up default options
    prcopt_t prcopt = prcopt_default;
    solopt_t solopt = solopt_default;
    filopt_t filopt = { "" };
    gtime_t ts = { 0 }, te = { 0 };
    double tint = 0.0;

    // Input and output file paths
    char *infile[3], *outfile = (char*)"";
	infile[0] = (char*)rover_path.c_str();
    infile[1] = (char*)base_path.c_str();
    infile[2] = (char*)nav_path.c_str();
    outfile = (char*)rtklib_path.c_str();

    // Reset system options and load custom configuration
    resetsysopts();
    loadopts((char*)config_path.c_str(), sysopts);
    getsysopts(&prcopt, &solopt, &filopt);
    
    // Store the configuration into the global option variable
    option = prcopt;

    // Perform processing
    postpos(ts, te, tint, 0.0, &prcopt, &solopt, &filopt, infile, 3, outfile, "", "");

    // Close the trace log file
    traceclose();

    return 0;
}
