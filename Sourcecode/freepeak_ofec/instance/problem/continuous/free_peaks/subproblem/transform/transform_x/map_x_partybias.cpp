#include "map_x_partybias.h"
#include "instance/problem/continuous/free_peaks/free_peaks.h"
#include <cmath>

namespace ofec {
    namespace free_peaks {

        void MapXPartyBias::initialize(Problem* prob, const std::string& subspace_name, const ParameterMap& param) {
            X_TransformBase::initialize(prob, subspace_name, param);

            m_party_id = param.at("party_id").getValue<int>();
            m_magnitude = param.at("magnitude").getValue<Real>();
            m_rotation_angle = param.at("rotation_angle").getValue<Real>();
            m_condition_number = param.at("condition_number").getValue<Real>();

            int dim = prob->numberVariables();
            m_shift.resize(dim, 0.0);

            // Create party-specific shift
            for (int i = 0; i < dim; ++i) {
                m_shift[i] = m_magnitude * std::sin(m_party_id * 2.1 + i * 1.3);
            }

            // Initialize rotation matrix
            m_rotation_matrix.resize(dim, dim);
            m_rotation_matrix.identify();

            if (dim >= 2 && m_rotation_angle != 0.0) {
                // Create 2D rotation in first two dimensions
                double c = std::cos(m_rotation_angle);
                double s = std::sin(m_rotation_angle);
                m_rotation_matrix[0][0] = static_cast<Real>(c);
                m_rotation_matrix[0][1] = static_cast<Real>(-s);
                m_rotation_matrix[1][0] = static_cast<Real>(s);
                m_rotation_matrix[1][1] = static_cast<Real>(c);
            }

            // Initialize scaling matrix for ill-conditioning
            m_scaling_matrix.resize(dim, dim);
            m_scaling_matrix.identify();

            if (m_condition_number > 1.0 && dim >= 2) {
                // Apply different scaling to create ill-conditioning
                double scale_factor = std::pow(m_condition_number, 1.0 / (dim - 1));
                for (int i = 0; i < dim; ++i) {
                    m_scaling_matrix[i][i] = static_cast<Real>(std::pow(scale_factor, i));
                }
            }

            m_initialized = true;
        }

        void MapXPartyBias::transfer(std::vector<Real>& x) const {
            if (!m_initialized || x.empty()) return;

            std::vector<Real> transformed = x;

            // Apply transformations in sequence
            applyBias(transformed);
            applyRotation(transformed);
            applyIllConditioning(transformed);

            x = transformed;
        }

        void MapXPartyBias::applyBias(std::vector<Real>& x) const {
            for (size_t i = 0; i < x.size(); ++i) {
                x[i] += m_shift[i];
                wrapToUnitInterval(x[i]);
            }
        }

        void MapXPartyBias::applyRotation(std::vector<Real>& x) const {
            if (m_rotation_matrix.row() == 0) return;

            int n = static_cast<int>(x.size());
            std::vector<Real> rotated(n, 0.0);

            // Shift to center
            for (int i = 0; i < n; ++i) {
                x[i] -= m_shift[i];
            }

            // Apply rotation: y = R * x
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    rotated[i] += m_rotation_matrix[i][j] * x[j];
                }
            }

            // Shift back and wrap
            for (int i = 0; i < n; ++i) {
                rotated[i] += m_shift[i];
                wrapToUnitInterval(rotated[i]);
            }

            swap(rotated, x);
        }

        void MapXPartyBias::applyIllConditioning(std::vector<Real>& x) const {
            if (m_scaling_matrix.row() == 0) return;

            int n = static_cast<int>(x.size());
            std::vector<Real> scaled(n, 0.0);

            // Shift to center
            for (int i = 0; i < n; ++i) {
                x[i] -= m_shift[i];
            }

            // Apply scaling: y = S * x
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    scaled[i] += m_scaling_matrix[i][j] * x[j];
                }
                // Wrap after scaling
                scaled[i] += m_shift[i];
                wrapToUnitInterval(scaled[i]);
            }

            swap(scaled, x);
        }

        void MapXPartyBias::wrapToUnitInterval(Real& val) const {
            while (val < 0.0) val += 1.0;
            while (val > 1.0) val -= 1.0;
        }

    } // namespace free_peaks
} // namespace ofec