#pragma once
#include "../global/global.h"

using namespace std;
using namespace Eigen;

class CVH_Factor : public ceres::CostFunction {

public:
	explicit CVH_Factor(double std) {
		_std = std;

		*mutable_parameter_block_sizes() = std::vector<int>{3, 3};
		set_num_residuals(3);
	};

	bool Evaluate(const double* const* parameters, double* residuals, double** jacobians) const override {

		// velocity vectors
		MatrixXd vf0(3, 1);
		for (int i = 0; i < 3; i++) {
			vf0(i, 0) = parameters[0][i];
		}

		MatrixXd vf1(3, 1);
		for (int i = 0; i < 3; i++) {
			vf1(i, 0) = parameters[1][i];
		}

		// Compute residual
		Map<Matrix<double, 3, 1>> residual(residuals);
		residual = vf0 - vf1; // residual = v0 - v1
		residual /= sqrt(_std); // Normalize by standard deviation

		// Compute Jacobians
		if (jacobians) {
			if (jacobians[0]) {
				Map<Matrix<double, 3, 3, Eigen::RowMajor>> jacobian0(jacobians[0]);
				jacobian0.setZero();
				jacobian0 = Matrix3d::Identity();
				jacobian0 /= sqrt(_std); // Jacobian for vf0
			}
			if (jacobians[1]) {
				Map<Matrix<double, 3, 3, Eigen::RowMajor>> jacobian1(jacobians[1]);
				jacobian1.setZero();
				jacobian1 = -Matrix3d::Identity();
				jacobian1 /= sqrt(_std); // Jacobian for vf1
			}
		}
		return true;
	}

public:
	double _std;
};