#ifndef OFEC_FREE_PEAKS_MAP_X_PARTY_BIAS_H
#define OFEC_FREE_PEAKS_MAP_X_PARTY_BIAS_H

/*
 * MapXPartyBias - Party-Specific Variable Transformation
 *
 * In Multiparty Multimodal Optimization (MPMMO), each decision-maker
 * (party) perceives the shared fitness landscape from a different perspective.
 * This transform applies three sequential, party-indexed operations:
 *   1. Bias shift    - translates the coordinate frame per party
 *   2. Rotation      - introduces non-separability (axis coupling)
 *   3. Ill-conditioning - stretches axes to create varying curvature
 */

#include "transform_x_base.h"
#include "../../../../../../../utility/linear_algebra/matrix.h"

#include <cmath>
#include <vector>

namespace ofec {
    namespace free_peaks {

        class MapXPartyBias : public X_TransformBase {
            
            OFEC_CONCRETE_INSTANCE(MapXPartyBias)

        protected:
            // Parameters registered with the input-parameter subsystem
            int  m_party_id = 0;
            Real m_magnitude = 5.0;   // shift amplitude (domain units)
            Real m_rotation_angle = 0.0;    // plane-rotation angle (radians)
            Real m_condition_number = 1.0;  // axis-stretch ratio (>= 1)

            // Derived state
            bool              m_initialized = false;
            std::vector<Real> m_shift;           // party-specific bias vector
            Matrix            m_rotation_matrix; // dim x dim orthogonal rotation
            Matrix            m_scaling_matrix;  // dim x dim diagonal scaling

        public:
            // Called automatically by the macro-generated constructor.
            // Registers parameters so the framework can persist / restore them.
            void addInputParameters();

            // TransformBase::initialize() first (sets m_pro, m_subspace_name,
            // runs m_input_parameters.input(param)), then build derived matrices.
            void initialize(Problem* pro,
                const std::string& subspace_name,
                const ParameterMap& param) override;

            // Core interface: maps x in-place before the fitness evaluation.
            void transfer(std::vector<Real>& x,
                const std::vector<Real>& var) override;

        private:
            // Sub-steps of transfer(), applied in this order:
            void applyBias(std::vector<Real>& x) const;
            void applyRotation(std::vector<Real>& x) const;
            void applyIllConditioning(std::vector<Real>& x) const;

        };

    } // namespace free_peaks
} // namespace ofec

#endif // OFEC_FREE_PEAKS_MAP_X_PARTY_BIAS_H