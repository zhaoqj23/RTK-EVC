#pragma once
#include "../global/global.h"

using namespace std;
using namespace Eigen;

class DD_P_Factor : public ceres::CostFunction {

public:

	explicit DD_P_Factor(node_info* node, int index)
		: _node(move(node)),
		_ddmeas(move(node->ddmeas[index])) {

		_index = index;
		_std = sqrt(_ddmeas.var); // Standard deviation from the measurement variance

		*mutable_parameter_block_sizes() = std::vector<int>{3}; // One parameter: 3D position

		set_num_residuals(1);
	};

	bool Evaluate(const double* const* parameters, double* residuals, double** jacobians) const override {

		// parameter vectors
		MatrixXd xf(3, 1);
		for (int i = 0; i < 3; i++) {
			xf(i, 0) = parameters[0][i];
		}

		// Reference and unreference satellite information
		double rs_ref[6], rs_uref[6], dts_ref[2], dts_uref[2];
		double var_ref, var_uref;
		int svh_ref, svh_uref;
		int nf = option.ionoopt == IONOOPT_IFLC ? 1 : option.nf;

		// Restore data from node information
		memcpy(rs_ref, _node->rs + 6 * _node->iu[_ddmeas.ref_index], sizeof(double) * 6);
		memcpy(dts_ref, _node->dts + 2 * _node->iu[_ddmeas.ref_index], sizeof(double) * 2);
		var_ref = _node->var[_node->iu[_ddmeas.ref_index]];
		svh_ref = _node->svh[_node->iu[_ddmeas.ref_index]];

		memcpy(rs_uref, _node->rs + 6 * _node->iu[_ddmeas.uref_index], sizeof(double) * 6);
		memcpy(dts_uref, _node->dts + 2 * _node->iu[_ddmeas.uref_index], sizeof(double) * 2);
		var_uref = _node->var[_node->iu[_ddmeas.uref_index]];
		svh_uref = _node->svh[_node->iu[_ddmeas.uref_index]];

		// Initialize matrices
		MatrixXd y_ref = MatrixXd::Zero(nf * 2, 1);
		MatrixXd e_ref = MatrixXd::Zero(3, 1);
		MatrixXd azel_ref = MatrixXd::Zero(2, 1);
		MatrixXd freq_ref = MatrixXd::Zero(nf, 1);

		MatrixXd y_uref = MatrixXd::Zero(nf * 2, 1);
		MatrixXd e_uref = MatrixXd::Zero(3, 1);
		MatrixXd azel_uref = MatrixXd::Zero(2, 1);
		MatrixXd freq_uref = MatrixXd::Zero(nf, 1);

		// Compute the zd residuals
		zdres(0, &(_node->obs[_node->iu[_ddmeas.ref_index]]), 1, rs_ref, dts_ref, &(var_ref), &(svh_ref),
			&(navs), xf.data(), &(option), 0, y_ref.data(), e_ref.data(), azel_ref.data(), freq_ref.data());

		zdres(0, &(_node->obs[_node->iu[_ddmeas.uref_index]]), 1, rs_uref, dts_uref, &(var_uref), &(svh_uref),
			&(navs), xf.data(), &(option), 0, y_uref.data(), e_uref.data(), azel_uref.data(), freq_uref.data());


		// Calculate the dd residuals
		Map<Matrix<double, 1, 1>> residual(residuals);
		int f = _ddmeas.f;
		
		residual(0, 0) = y_ref(f + nf, 0) - y_uref(f + nf, 0) - _node->d_base[_index];

		trace(4, "epoch%d: dd_p_factor P%d %d - %d: raw residual: %.10lf standardized residual: %.10lf\n", _node->nepoch, f + 1, _node->nsat[_ddmeas.ref_index], _node->nsat[_ddmeas.uref_index], residual(0, 0), residual(0, 0) / _std);

		residual(0, 0) /= _std; // Standardize the residual

		// Check for invalid residuals
		if (residual.hasNaN() || !residual.allFinite()) {
			trace(3, "dd_l_factor: residual has nan!\n");
		}

		// Compute Jacobians
		if (jacobians) {
			if (jacobians[0]) {
				Map<Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_xf(jacobians[0]);
				jacobian_xf.setZero();
				jacobian_xf.block(0, 0, 1, 3) = e_ref.transpose() - e_uref.transpose();
				jacobian_xf = jacobian_xf / _std;
			}
		}
		return true;
	}

public:
	node_info* _node;
	ddinfo _ddmeas;
	int _index;
	double _std;
};