#pragma once
#include "../global/global.h"

using namespace std;
using namespace Eigen;

class POS_Fix_Prior : public ceres::CostFunction {

public:
	explicit POS_Fix_Prior(node_info *node, double var_prior)
		: _node(move(node)) {

		_var_prior = var_prior; // Standard deviation

		*mutable_parameter_block_sizes() = std::vector<int>{3}; // One parameters: 3D position

		set_num_residuals(3);
	};

	bool Evaluate(const double* const* parameters, double* residuals, double** jacobians) const override {

		Matrix<double, 3, 1> vf = Matrix<double, 3, 1>::Zero();
		for (int i = 0; i < 3; i++) {
			vf(i, 0) = parameters[0][i];
		}

		// Compute the residuals
		Map<Matrix<double, 3, 1>> residual(residuals);
		for (int i = 0; i < 3; i++) {
			residual(i, 0) = vf(i, 0) - _node->xa[i];
		}
		residual /= sqrt(_var_prior); // Standardize the residual
		
		// Compute Jacobians
		if (jacobians) {
			if (jacobians[0]) {
				Map<Matrix<double, 3, 3, Eigen::RowMajor>> jacobian(jacobians[0]);
				jacobian = Matrix3d::Identity();
				jacobian /= sqrt(_var_prior);
			}
		}
		return true;
	}

public:
	double _var_prior;
	node_info *_node;
};