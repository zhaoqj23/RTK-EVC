#pragma once
#include "../global/global.h"

using namespace std;
using namespace Eigen;
class VIC_Factor : public ceres::CostFunction {

public:
	explicit VIC_Factor(node_info* node, node_info* noden)
		: _node(move(node)), _noden(move(noden)) {

		*mutable_parameter_block_sizes() = std::vector<int>{3, 3};

		set_num_residuals(3);

		_sqrt_information.resize(3, 3);
		_sqrt_information.setIdentity();
		_Q.resize(3, 3);
		_Sj = 0.00003; // Power Spectral Density
		_ddt = timediff(noden->time, _node->time);
		_Q.setIdentity();
		_Q.block(0, 0, 3, 3) = _Sj * _ddt * _ddt * Eigen::Matrix3d::Identity();
		_sqrt_information = (_Q.inverse()).llt().matrixL().transpose();
	};

	bool Evaluate(const double* const* parameters, double* residuals, double** jacobians) const override {

		// parameter vectors
		Map<Matrix<double, 3, 1>> residual(residuals);

		MatrixXd xf0(3, 1);
		for (int i = 0; i < 3; i++) {
			xf0(i, 0) = parameters[0][i];
		}

		MatrixXd xf1(3, 1);
		for (int i = 0; i < 3; i++) {
			xf1(i, 0) = parameters[1][i];
		}

		// Compute the residuals
		for (int i = 0; i < 3; i++) {
			residual(i, 0) = xf1(i, 0) - xf0(i, 0) - (_node->xf[i + 3] + _noden->xf[i + 3]) * _ddt / 2;
		}
		trace(4, "epoch%d: VIC Factor: Raw Residual:", _noden->nepoch); tracemat(3, residual.data(), 1, 3, 5, 5);
		trace(4, "epoch%d: VIC Factor: Residual Norm: %.3lf\n", _noden->nepoch, (_sqrt_information * residual).norm());
		residual = _sqrt_information * residual;

		// Compute Jacobians
		if (jacobians) {
			if (jacobians[0]) {
				Map<Matrix<double, 3, 3, Eigen::RowMajor>> jacobian_xf0(jacobians[0]);
				jacobian_xf0.setZero();
				jacobian_xf0 = -Matrix3d::Identity();
				jacobian_xf0 = _sqrt_information * jacobian_xf0;
			}

			if (jacobians[1]) {
				Map<Matrix<double, 3, 3, Eigen::RowMajor>> jacobian_xf1(jacobians[1]);
				jacobian_xf1.setZero();
				jacobian_xf1 = Matrix3d::Identity();
				jacobian_xf1 = _sqrt_information * jacobian_xf1;
			}
		}

		return true;
	}

public:
	node_info* _node, *_noden;
	MatrixXd _Q, _sqrt_information;
	double _ddt, _Sj;
};
