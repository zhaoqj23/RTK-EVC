#pragma once
#include "../global/global.h"

using namespace std;
using namespace Eigen;
class Dopp_Factor : public ceres::CostFunction {

public:
	explicit Dopp_Factor(node_info* node, int satindex, int f)
		: _node(move(node)) {

		_satindex = satindex; _f = f;

		*mutable_parameter_block_sizes() = std::vector<int>{3, 1}; // Two parameters: 3D velocity and clock frequency deviation

		set_num_residuals(1);
	};

	bool Evaluate(const double* const* parameters, double* residuals, double** jacobians) const override {
		
		int i, j;
		Vector3d vf;
		double vs[3], e[3], freq, rate, clockshift;

		// parameter vectors
		for (i = 0; i < 3; i++) {
			vf(i) = parameters[0][i];
		}
		clockshift = parameters[1][0];
		
		Map<Matrix<double, Dynamic, Dynamic>> residual(residuals, 1, 1);
		
		// Restore data from node information
		freq = sat2freq(_node->obs[_satindex].sat, _node->obs[_satindex].code[_f], &navs);
		for (j = 0; j < 3; j++) {
			vs[j] = _node->rs[j + 3 + _satindex * 6] - vf(j);
			e[j] = _node->e[j + _satindex * 3];
		}

		rate = dot(vs, e, 3) + OMGE / CLIGHT * (_node->rs[4 + _satindex * 6] * _node->xf[0] + _node->rs[1 + _satindex * 6] * vf(0) -
			_node->rs[3 + _satindex * 6] * _node->xf[1] - _node->rs[_satindex * 6] * vf(1));

		// Compute the residuals
		residual(0, 0) = -_node->obs[_satindex].D[_f] * CLIGHT / freq - (rate + clockshift - CLIGHT * _node->dts[1 + _satindex * 2]);

		// Compute Jacobians
		if (jacobians) {
			if (jacobians[0]) {
				Map<Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_vf(jacobians[0]);
				jacobian_vf.setZero();
				for (j = 0; j < 3; j++) {
					jacobian_vf(0, j) = e[j];
				}
				jacobian_vf = jacobian_vf;
			}
			if (jacobians[1]) {
				Map<Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_clockshift(jacobians[1]);
				jacobian_clockshift.setZero();
				jacobian_clockshift(0, 0) = -1;
				jacobian_clockshift = jacobian_clockshift;
			}
		}
		return true;
	}

private:
	node_info *_node;
	int _satindex, _f;
};