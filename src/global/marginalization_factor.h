#pragma once
#include "../global/global.h"
#include <memory>

#include "marginalization_info.h"

using namespace std;
using namespace Eigen;


class MarginalizationFactor : public ceres::CostFunction {

public:
    MarginalizationFactor() = delete;
    explicit MarginalizationFactor(std::shared_ptr<MarginalizationInfo> marg_info)
        : marg_info_(std::move(marg_info)) {

        for (auto size : marg_info_->remainedBlockSize()) {
            mutable_parameter_block_sizes()->push_back(size);
        }

        set_num_residuals(marg_info_->remainedSize());
    }



    bool Evaluate(const double *const *parameters, double *residuals, double **jacobians) const override {
        int marginalizaed_size = marg_info_->marginalizedSize();
        int remained_size      = marg_info_->remainedSize();

        const std::vector<int> &remained_block_index     = marg_info_->remainedBlockIndex();
        const std::vector<int> &remained_block_size      = marg_info_->remainedBlockSize();
        const std::vector<double *> &remained_block_data = marg_info_->remainedBlockData();

        Eigen::VectorXd dx(remained_size);
        for (size_t i = 0; i < remained_block_size.size(); i++) {
            int size  = remained_block_size[i];
            int index = remained_block_index[i] - marginalizaed_size;

            Eigen::VectorXd x  = Eigen::Map<const Eigen::VectorXd>(parameters[i], size);
            Eigen::VectorXd x0 = Eigen::Map<const Eigen::VectorXd>(remained_block_data[i], size);

            // dx = x - x0
            if (size == POSE_GLOBAL_SIZE) {
                Eigen::Quaterniond dq(Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() *
                                      Eigen::Quaterniond(x(6), x(3), x(4), x(5)));

                dx.segment(index, 3)     = x.head<3>() - x0.head<3>();
                dx.segment(index + 3, 3) = 2.0 * dq.vec();
                if (dq.w() < 0) {
                    dx.segment<3>(index + 3) = -2.0 * dq.vec();
                }
            } else {
                dx.segment(index, size) = x - x0;
            }
        }

        // e = e0 + J0 * dx
        Eigen::Map<Eigen::VectorXd>(residuals, remained_size) =
            marg_info_->linearizedResiduals() + marg_info_->linearizedJacobians() * dx;

        if (jacobians) {

            for (size_t i = 0; i < remained_block_size.size(); i++) {
                if (jacobians[i]) {
                    int size       = remained_block_size[i];
                    int index      = remained_block_index[i] - marginalizaed_size;
                    int local_size = marg_info_->localSize(size);

                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobian(
                        jacobians[i], remained_size, size);

                    // J = J0
                    jacobian.setZero();
                    jacobian.leftCols(local_size) = marg_info_->linearizedJacobians().middleCols(index, local_size);
                }
            }
        }

        return true;
    }

public:
    std::shared_ptr<MarginalizationInfo> marg_info_;
};
