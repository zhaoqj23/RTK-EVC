# RTK-EVC: Real-Time Kinematic Positioning with Enhanced Velocity Constraints via Factor Graph Optimization in GNSS-challenged Environments
We open-source the source code and experiment dataset of a novel Real-time Kinematic positioning architecture with Enhanced Velocity Constraints (RTK-EVC) via Factor Graph Optimization in GNSS-challenged Environments.
## How to Use This Library
The library was developed based on  **[RTKLIB](https://github.com/tomojitakasu/RTKLIB)**  2.4.3 b34 and tested in Ubuntu 20.04. The results may be different for different OS. If any problem was found by you, please propose an issue or report it to zhaoqj23@mails.tsinghua.edu.cn.
### Requirements
1) Ubuntu 20.04 with the newest compiler is recommended.
2) Eigen3
3) Ceres Solver

Eigen 3.3.7 and Ceres 2.1.0 have been tested and work properly. The dataset path and config path need to be modified in ```src/global/global.cpp```. The trace log path need to be modified in ```src/main.cpp```.
### Clone the repository
```git clone https://github.com/zhaoqj23/RTK-EVC.git```
### Download the dataset
The dataset of RTK-EVC can be obtained via **[Onedrive](https://1drv.ms/f/c/f185ecb4d7c54710/EmEiek9-Cz1PmxUmrju4Xy0B2WJMHYBJSbt1KmWuKx3Ndg?e=hn8mwY)**. The experiment was carried out in Beijing, China, around Tsinghua University and Zhongguancun, with multiple GNSS-challenged scenarios. IMU-ISA-100C and Inertial Explorer were used to get the ground truth.
### Build the library
```
cd ~/RTK-EVC
mkdir build && cd build
cmake ..
make -j8
```
### Run demo
```
cd ~/RTK-EVC
./bin/RTK-EVC
```
