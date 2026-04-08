#include "map_x_partybias.h"
#include "instance/problem/continuous/free_peaks/free_peaks.h"
#include <cmath>

namespace ofec {
    namespace free_peaks {

        void MapXPartyBias::addInputParameters() {
            m_input_parameters.add("party_id", new InputParameterValueTypeInt(m_party_id));
            m_input_parameters.add("magnitude", new EnumeratedReal(m_magnitude, m_recommended_magnitudes, 0.02));
            m_input_parameters.add("rotation_angle", new EnumeratedReal(m_rotation_angle, m_recommended_angles, 0.0));
            m_input_parameters.add("condition_number", new EnumeratedReal(m_condition_number, m_recommended_conditions, 1.0));
        }

        void MapXPartyBias::initialize(Problem* prob, const std::string& subspace_name, const ParameterMap& param) {
            // Call base initialize
            X_TransformBase::initialize(prob, subspace_name, param);

            // Extract party-specific parameters using has() and get()
            if (param.has("party_id")) {
                m_party_id = param.get<int>("party_id");
            }

            if (param.has("magnitude")) {
                m_magnitude = param.get<Real>("magnitude");
            }

            if (param.has("rotation_angle")) {
                m_rotation_angle = param.get<Real>("rotation_angle");
            }

            if (param.has("condition_number")) {
                m_condition_number = param.get<Real>("condition_number");
            }

            // Initialize rotation matrix if rotation is specified
            int numVar = prob->numberVariables();
            if (m_rotation_angle != 0.0 && numVar >= 2) {
                m_rotation_matrix.resize(numVar, numVar);
                m_rotation_matrix.identify();
                
                // Create 2D rotation in the first two dimensions
                // For higher dimensions, this creates rotation in x1-x2 plane
                Real cos_theta = std::cos(m_rotation_angle);
                Real sin_theta = std::sin(m_rotation_angle);
                
                m_rotation_matrix[0][0] = cos_theta;
                m_rotation_matrix[0][1] = -sin_theta;
                m_rotation_matrix[1][0] = sin_theta;
                m_rotation_matrix[1][1] = cos_theta;
                
                // Apply ill-conditioning scaling if specified
                if (m_condition_number > 1.0) {
                    // Scale axes to create ill-conditioning
                    Real scale_factor = std::sqrt(m_condition_number);
                    for (int i = 0; i < numVar; ++i) {
                        Real axis_scale = (i == 0) ? scale_factor : (i == 1) ? 1.0/scale_factor : 1.0;
                        for (int j = 0; j < numVar; ++j) {
                            m_rotation_matrix[i][j] *= axis_scale;
                        }
                    }
                }
            }

            // Set shift to center of search space
            m_shift.resize(numVar, 0.5);
        }

        void MapXPartyBias::transfer(std::vector<Real>& x, const std::vector<Real>& var) {
            // Make a copy to work with
            std::vector<Real> transformed = x;

            // Step 1: Apply party-specific bias (position offset)
            applyBias(transformed);

            // Step 2: Apply rotation (creates non-separability)
            if (m_rotation_angle != 0.0) {
                applyRotation(transformed);
            }

            // Step 3: Apply ill-conditioning (different scaling per dimension)
            if (m_condition_number > 1.0) {
                applyIllConditioning(transformed);
            }

            // Swap result back to output
            swap(transformed, x);
        }

        void MapXPartyBias::applyBias(std::vector<Real>& x) {
            // Create party-specific offset (different for each party)
            // This implements the core idea: each party sees peaks at different locations
            Real party_offset = static_cast<Real>(m_party_id) * m_magnitude;

            for (size_t i = 0; i < x.size(); ++i) {
                // Dimension-specific bias: creates more complex landscape differences
                Real dim_bias = party_offset * static_cast<Real>(i + 1);
                x[i] += dim_bias;
                
                // Wrap to [0, 1] range to maintain feasibility
                wrapToUnitInterval(x[i]);
            }
        }

        void MapXPartyBias::applyRotation(std::vector<Real>& x) {
            if (m_rotation_matrix.row() == 0) return;

            int n = static_cast<int>(x.size());
            std::vector<Real> rotated(n, 0.0);

            // First, shift to center (rotate around center point)
            for (int i = 0; i < n; ++i) {
                x[i] -= m_shift[i];
            }

            // Apply rotation matrix multiplication: y = R * x
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    rotated[i] += m_rotation_matrix[i][j] * x[j];
                }
            }

            // Shift back
            for (int i = 0; i < n; ++i) {
                rotated[i] += m_shift[i];
                wrapToUnitInterval(rotated[i]);
            }

            swap(rotated, x);
        }

        void MapXPartyBias::applyIllConditioning(std::vector<Real>& x) {
            // Apply dimensional scaling to create ill-conditioning
            // This makes optimization harder by stretching the landscape
            Real scale_factor = std::sqrt(m_condition_number);
            
            for (size_t i = 0; i < x.size(); ++i) {
                // Center, scale, then re-center
                x[i] -= m_shift[i];
                
                // Different scaling per dimension (geometric progression)
                Real dim_scale = std::pow(scale_factor, static_cast<Real>(2 * static_cast<int>(i) / static_cast<int>(x.size()) - 1));
                x[i] *= dim_scale;
                
                x[i] += m_shift[i];
                wrapToUnitInterval(x[i]);
            }
        }

        void MapXPartyBias::wrapToUnitInterval(Real& val) {
            // Wrap value to [0, 1] range
            while (val < static_cast<Real>(0.0)) val += static_cast<Real>(1.0);
            while (val > static_cast<Real>(1.0)) val -= static_cast<Real>(1.0);
        }

        std::string MapXPartyBias::getName() const {
            return "MapXPartyBias";
        }

    } // namespace free_peaks
} // namespace ofec