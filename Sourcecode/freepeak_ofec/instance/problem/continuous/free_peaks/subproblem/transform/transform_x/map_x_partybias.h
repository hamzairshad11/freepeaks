#pragma once
#include "transform_x_base.h"
#include "../transform_base.h" 

namespace ofec {
    class Problem;
    namespace free_peaks {
        class FreePeaks;
        class Subproblem;

        class MapXPartyBias : public X_TransformBase {
        protected:
            std::vector<Real> m_bias;

        public:
            MapXPartyBias();
            virtual ~MapXPartyBias();

            void initialize(Problem* prob, const std::string& subspace_name, const ParameterMap& param);
            void transform(VariableVector<Real>& x, const std::string& subspace_name);
            std::string getName() const;  

        private:
            int m_party_id = 0;
            double m_magnitude = 0.0;
        };
    } // namespace free_peaks
} // namespace ofec