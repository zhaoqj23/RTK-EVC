#ifndef MARGINILAZATION_INFO_H
#define MARGINILAZATION_INFO_H

#include <memory>
#include <unordered_map>

#include "residual_block_info.h"

class MarginalizationInfo {

public:
    MarginalizationInfo() = default;

    ~MarginalizationInfo() {
        for (auto &block : parameter_block_data_)
            delete[] block.second;
    }

    bool isValid() const {
        return isvalid_;
    }

    static int localSize(int size) {
        return size == POSE_GLOBAL_SIZE ? POSE_LOCAL_SIZE : size;
    }

    static int globalSize(int size) {
        return size == POSE_LOCAL_SIZE ? POSE_GLOBAL_SIZE : size;
    }

    void addResidualBlockInfo(const std::shared_ptr<ResidualBlockInfo> &blockinfo) {
        factors_.push_back(blockinfo);

        const auto &parameter_blocks = blockinfo->parameterBlocks();
        const auto &block_sizes      = blockinfo->parameterBlockSizes();

        for (size_t k = 0; k < parameter_blocks.size(); k++) {
            parameter_block_size_[reinterpret_cast<long>(parameter_blocks[k])] = block_sizes[k];
        }

        for (int index : blockinfo->marginalizationParametersIndex()) {
            parameter_block_index_[reinterpret_cast<long>(parameter_blocks[index])] = 0;
        }
    }

    bool marginalization() {

        if (!updateParameterBlocksIndex()) {
            isvalid_ = false;

            releaseMemory();

            return false;
        }

        preMarginalization();

        constructEquation();

        schurElimination();

        linearization();

        releaseMemory();

        return true;
    }

    std::vector<double *> getParamterBlocks(std::unordered_map<long, double *> &address) {
        std::vector<double *> remained_block_addr;

        remained_block_data_.clear();
        remained_block_index_.clear();
        remained_block_size_.clear();

        for (const auto &block : parameter_block_index_) {

            if (block.second >= marginalized_size_) {
                remained_block_data_.push_back(parameter_block_data_[block.first]);
                remained_block_size_.push_back(parameter_block_size_[block.first]);
                remained_block_index_.push_back(parameter_block_index_[block.first]);
                remained_block_addr.push_back(address[block.first]);
            }
        }
        return remained_block_addr;
    }

    const Eigen::MatrixXd &linearizedJacobians() {
        return linearized_jacobians_;
    }

    const Eigen::VectorXd &linearizedResiduals() {
        return linearized_residuals_;
    }

    int marginalizedSize() const {
        return marginalized_size_;
    }

    int remainedSize() const {
        return remained_size_;
    }

    const std::vector<int> &remainedBlockSize() {
        return remained_block_size_;
    }

    const std::vector<int> &remainedBlockIndex() {
        return remained_block_index_;
    }

    const std::vector<double *> &remainedBlockData() {
        return remained_block_data_;
    }

private:
    void linearization() {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(Hp_);
        Eigen::VectorXd S = Eigen::VectorXd((saes2.eigenvalues().array() > EPS).select(saes2.eigenvalues().array(), 0));
        Eigen::VectorXd S_inv =
            Eigen::VectorXd((saes2.eigenvalues().array() > EPS).select(saes2.eigenvalues().array().inverse(), 0));

        Eigen::VectorXd S_sqrt     = S.cwiseSqrt();
        Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();

        // J0 = S^{1/2} * V^T
        linearized_jacobians_ = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();
        // e0 = -{J0^T}^{-1} * bp = - S^{-1/2} * V^T * bp
        linearized_residuals_ = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose() * -bp_;
    }

    // Hp * dx_r = bp
    void schurElimination() {
        // H0 * dx = b0
        Eigen::MatrixXd Hmm = 0.5 * (H0_.block(0, 0, marginalized_size_, marginalized_size_) +
                                     H0_.block(0, 0, marginalized_size_, marginalized_size_).transpose());
        Eigen::MatrixXd Hmr = H0_.block(0, marginalized_size_, marginalized_size_, remained_size_);
        Eigen::MatrixXd Hrm = H0_.block(marginalized_size_, 0, remained_size_, marginalized_size_);
        Eigen::MatrixXd Hrr = H0_.block(marginalized_size_, marginalized_size_, remained_size_, remained_size_);
        Eigen::VectorXd bmm = b0_.segment(0, marginalized_size_);
        Eigen::VectorXd brr = b0_.segment(marginalized_size_, remained_size_);

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(Hmm);
        Eigen::MatrixXd Hmm_inv =
            saes.eigenvectors() *
            Eigen::VectorXd((saes.eigenvalues().array() > EPS).select(saes.eigenvalues().array().inverse(), 0))
                .asDiagonal() *
            saes.eigenvectors().transpose();

        // Hp = Hrr - Hrm * Hmm^-1 * Hmr
        Hp_ = Hrr - Hrm * Hmm_inv * Hmr;
        // bp = br - Hrm * Hmm^-1 * bm
        bp_ = brr - Hrm * Hmm_inv * bmm;
    }

    void constructEquation() {
        H0_ = Eigen::MatrixXd::Zero(local_size_, local_size_);
        b0_ = Eigen::VectorXd::Zero(local_size_);

        for (const auto &factor : factors_) {
            for (size_t i = 0; i < factor->parameterBlocks().size(); i++) {
                int row0 = parameter_block_index_[reinterpret_cast<long>(factor->parameterBlocks()[i])];
                int rows = parameter_block_size_[reinterpret_cast<long>(factor->parameterBlocks()[i])];
                rows     = (rows == POSE_GLOBAL_SIZE) ? POSE_LOCAL_SIZE : rows;

                Eigen::MatrixXd jacobian_i = factor->jacobians()[i].leftCols(rows);
                for (size_t j = i; j < factor->parameterBlocks().size(); ++j) {
                    int col0 = parameter_block_index_[reinterpret_cast<long>(factor->parameterBlocks()[j])];
                    int cols = parameter_block_size_[reinterpret_cast<long>(factor->parameterBlocks()[j])];
                    cols     = (cols == POSE_GLOBAL_SIZE) ? POSE_LOCAL_SIZE : cols;

                    Eigen::MatrixXd jacobian_j = factor->jacobians()[j].leftCols(cols);

                    // H = J^T * J
                    if (i == j) {
                        // Hmm, Hrr
                        H0_.block(row0, col0, rows, cols) += jacobian_i.transpose() * jacobian_j;
                    } else {
                        // Hmr, Hrm = Hmr^T
                        H0_.block(row0, col0, rows, cols) += jacobian_i.transpose() * jacobian_j;
                        H0_.block(col0, row0, cols, rows) = H0_.block(row0, col0, rows, cols).transpose();
                    }
                }
                // b = - J^T * e
                b0_.segment(row0, rows) -= jacobian_i.transpose() * factor->residuals();
            }
        }
    }

    bool updateParameterBlocksIndex() {
        int index = 0;
        for (auto &block : parameter_block_index_) {
            block.second = index;
            index += localSize(parameter_block_size_[block.first]);
        }
        marginalized_size_ = index;

        for (const auto &block : parameter_block_size_) {
            if (parameter_block_index_.find(block.first) == parameter_block_index_.end()) {
                parameter_block_index_[block.first] = index;
                index += localSize(block.second);
            }
        }
        remained_size_ = index - marginalized_size_;

        local_size_ = index;

        return marginalized_size_ > 0;
    }

    void preMarginalization() {
        for (const auto &factor : factors_) {
            factor->Evaluate();

            std::vector<int> block_sizes = factor->parameterBlockSizes();
            for (size_t k = 0; k < block_sizes.size(); k++) {
                long addr = reinterpret_cast<long>(factor->parameterBlocks()[k]);
                int size  = block_sizes[k];

                if (parameter_block_data_.find(addr) == parameter_block_data_.end()) {
                    auto *data = new double[size];
                    memcpy(data, factor->parameterBlocks()[k], sizeof(double) * size);
                    parameter_block_data_[addr] = data;
                }
            }
        }
    }

    void releaseMemory() {
        factors_.clear();
    }

private:

    Eigen::MatrixXd H0_, Hp_;
    Eigen::VectorXd b0_, bp_;

    std::unordered_map<long, int> parameter_block_size_;
    std::unordered_map<long, int> parameter_block_index_;
    std::unordered_map<long, double *> parameter_block_data_;

    std::vector<int> remained_block_size_; 
    std::vector<int> remained_block_index_; 
    std::vector<double *> remained_block_data_;

    int marginalized_size_{0};
    int remained_size_{0};
    int local_size_{0};

    std::vector<std::shared_ptr<ResidualBlockInfo>> factors_;

    const double EPS = 1e-8;

    Eigen::MatrixXd linearized_jacobians_;
    Eigen::VectorXd linearized_residuals_;

    bool isvalid_{true};
};

#endif // MARGINILAZATION_INFO_H
