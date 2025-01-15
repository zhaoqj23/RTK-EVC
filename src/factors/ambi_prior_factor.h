#pragma once
#include "../global/global.h"

using namespace std;
using namespace Eigen;

class AMBI_Prior : public ceres::CostFunction {

public:

	explicit AMBI_Prior(double ddbias, double var_prior) {
		_ddbias = ddbias;
		_var_prior = var_prior;

		// Two parameters: SD ambiguity of the reference satellite and unreference satellite
		*mutable_parameter_block_sizes() = std::vector<int>{1, 1};  
		set_num_residuals(1);
	};


	bool Evaluate(const double* const* parameters, double* residuals, double** jacobians) const override {

		// Reference and unreference ambiguity
		Matrix<double, 1, 1> refbias = Matrix<double, 1, 1>::Zero();
		refbias(0, 0) = parameters[0][0];

		Matrix<double, 1, 1> urefbias = Matrix<double, 1, 1>::Zero();
		urefbias(0, 0) = parameters[1][0];

		// Calculate residual
		Map<Matrix<double, 1, 1>> residual(residuals);
		residual(0, 0) = refbias(0, 0) - urefbias(0, 0) - _ddbias;
		residual(0, 0) /= sqrt(_var_prior);  // Normalize by prior variance

		// Compute Jacobians
		if (jacobians) {
			if (jacobians[0]) {
				Map<Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_ref(jacobians[0]);
				jacobian_ref.setIdentity();
				jacobian_ref /= sqrt(_var_prior);  // Jacobian for reference ambiguity
			}
			if (jacobians[1]) {
				Map<Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_uref(jacobians[1]);
				jacobian_uref.setIdentity();
				jacobian_uref /= (-sqrt(_var_prior));  // Jacobian for unreference ambiguity
			}
		}

		return true;
	}

public:
	double _ddbias, _var_prior;
};
