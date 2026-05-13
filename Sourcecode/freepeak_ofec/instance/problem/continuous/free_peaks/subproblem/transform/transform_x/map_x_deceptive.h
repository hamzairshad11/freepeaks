#ifndef OFEC_MAP_X_DECEPTIVE_H
#define OFEC_MAP_X_DECEPTIVE_H

#include "transform_x_base.h"

namespace ofec {
    namespace free_peaks {

        class MapXDeceptive : public X_TransformBase {
            OFEC_CONCRETE_INSTANCE(MapXDeceptive)
        protected:
            Real m_alpha = 0.05;
            bool m_initialized = false;
        public:
            void addInputParameters();
            void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) override;
            void transfer(std::vector<Real>& x, const std::vector<Real>& var) override;
        };

    }
}

#endif