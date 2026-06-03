#include "map_x_rotation.h"
#include "../../../free_peaks.h"
#include <random>
#include <cmath>

namespace ofec::free_peaks {

    void MapXRotation::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
        MapXBase::initialize(pro, subspace_name, param);

        size_t dim = m_subspace_var_ranges.size();
        if (dim <= 1) return;

        unsigned seed = 0;
        if (param.has("seed")) {
            seed = static_cast<unsigned>(std::get<Real>(param.at("seed")));
        }
        else {
            seed = static_cast<unsigned>(std::hash<std::string>{}(subspace_name));
        }

        generateRandomRotationMatrix(dim, seed);
    }

    // Generate a random orthogonal matrix via Gram-Schmidt on a
    // Gaussian matrix.  Ensures det(R) = +1 (proper rotation).
    void MapXRotation::generateRandomRotationMatrix(size_t dim, unsigned seed) {
        m_rotation_matrix.assign(dim, std::vector<Real>(dim, Real(0.0)));

        std::mt19937 rng(seed);
        std::normal_distribution<Real> gauss(0.0, 1.0);

        // 1. Random Gaussian matrix
        for (size_t i = 0; i < dim; ++i)
            for (size_t j = 0; j < dim; ++j)
                m_rotation_matrix[i][j] = gauss(rng);

        // 2. Orthogonalize columns
        gramSchmidt(m_rotation_matrix);

        // 3. Force determinant +1 (if det == -1, flip last column)
        Real det = computeDeterminant(m_rotation_matrix);
        if (det < Real(0.0)) {
            for (size_t i = 0; i < dim; ++i)
                m_rotation_matrix[i][dim - 1] = -m_rotation_matrix[i][dim - 1];
        }
    }

    // Gram-Schmidt orthonormalization of columns
    void MapXRotation::gramSchmidt(std::vector<std::vector<Real>>& A) {
        size_t n = A.size();
        if (n == 0) return;

        for (size_t k = 0; k < n; ++k) {
            Real norm = 0.0;
            for (size_t i = 0; i < n; ++i)
                norm += A[i][k] * A[i][k];
            norm = std::sqrt(norm);

            if (norm < Real(1e-12)) {
                norm = Real(1e-12);
                A[k][k] = Real(1.0); // fallback to identity column
            }

            for (size_t i = 0; i < n; ++i)
                A[i][k] /= norm;

            for (size_t j = k + 1; j < n; ++j) {
                Real dot = 0.0;
                for (size_t i = 0; i < n; ++i)
                    dot += A[i][k] * A[i][j];
                for (size_t i = 0; i < n; ++i)
                    A[i][j] -= dot * A[i][k];
            }
        }
    }

    //LU-style determinant for low dimensions (2,3,4,...)
    Real MapXRotation::computeDeterminant(const std::vector<std::vector<Real>>& mat) {
        size_t n = mat.size();
        if (n == 0) return Real(1.0);
        if (n == 1) return mat[0][0];
        if (n == 2)
            return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
        if (n == 3)
            return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
            - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
            + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

        std::vector<std::vector<Real>> a = mat;
        Real det = Real(1.0);
        for (size_t i = 0; i < n; ++i) {
            size_t pivot = i;
            Real max_val = std::abs(a[i][i]);
            for (size_t r = i + 1; r < n; ++r) {
                if (std::abs(a[r][i]) > max_val) {
                    max_val = std::abs(a[r][i]);
                    pivot = r;
                }
            }
            if (max_val < Real(1e-12)) return Real(0.0);
            if (pivot != i) {
                std::swap(a[i], a[pivot]);
                det = -det;
            }
            det *= a[i][i];
            for (size_t r = i + 1; r < n; ++r) {
                Real factor = a[r][i] / a[i][i];
                for (size_t c = i; c < n; ++c)
                    a[r][c] -= factor * a[i][c];
            }
        }
        return det;
    }

    // Apply rotation: x_ = R * x_
    void MapXRotation::transfer(std::vector<Real>& x_, const std::vector<Real>& var) {
        if (m_rotation_matrix.empty() || x_.size() <= 1) return;

        size_t dim = x_.size();
        std::vector<Real> rotated(dim, Real(0.0));
        for (size_t i = 0; i < dim; ++i)
            for (size_t j = 0; j < dim; ++j)
                rotated[i] += m_rotation_matrix[i][j] * x_[j];

        x_ = std::move(rotated);
    }

} // namespace ofec::free_peaks