#pragma once
#define _USE_MATH_DEFINES
#include "transform_x_base.h"
#include "../transform_base.h" 
#include "../../../../../../../utility/linear_algebra/matrix.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ofec {
    class Problem;
    namespace free_peaks {
        class FreePeaks;
        class Subproblem;

        /// Each party (decision-maker) sees a different fitness landscape due to:
        /// 1. Position bias: Peaks appear at slightly different locations for each party
        /// 2. Dimensional scaling: Different conditioning per party
        /// 3. Rotation: Non-separable interactions between variables
        /// 
        /// Mathematical formulation: y = T_party(x) where T_party is party-specific
        /// This creates two distinct fitness landscapes: f_1(x) = f(T_1(x)) and f_2(x) = f(T_2(x))
        class MapXPartyBias : public X_TransformBase {
            OFEC_CONCRETE_INSTANCE(MapXPartyBias)
            
        protected:
            std::vector<Real> m_bias;
            Matrix m_rotation_matrix;
            Real m_condition_number = 1.0;
            Real m_rotation_angle = 0.0;
            std::vector<Real> m_shift;
            
            
            std::list<Real> m_recommended_magnitudes = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.1};
            std::list<Real> m_recommended_angles = {0.0, M_PI/12, M_PI/8, M_PI/6, M_PI/4, M_PI/3};
            std::list<Real> m_recommended_conditions = {1.0, 10.0, 50.0, 100.0, 500.0, 1000.0};

        public:
            void addInputParameters() override;
            void initialize(Problem* prob, const std::string& subspace_name, const ParameterMap& param) override;
            void transfer(std::vector<Real>& x, const std::vector<Real>& var) override;
            std::string getName() const override { return "MapXPartyBias"; }  

        private:
            int m_party_id = 0;
            Real m_magnitude = 0.0;
            
            // Helper methods
            void applyBias(std::vector<Real>& x);
            void applyRotation(std::vector<Real>& x);
            void applyIllConditioning(std::vector<Real>& x);
            void wrapToUnitInterval(Real& val);
        };
    } // namespace free_peaks
} // namespace ofec
