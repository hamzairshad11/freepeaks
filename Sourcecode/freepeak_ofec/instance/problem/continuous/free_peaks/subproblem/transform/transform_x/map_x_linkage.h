#ifndef OFEC_MAP_X_LINKAGE_H
#define OFEC_MAP_X_LINKAGE_H

#include "transform_x_base.h"

namespace ofec {
    namespace free_peaks {

        class MapXLinkage : public X_TransformBase {
            OFEC_CONCRETE_INSTANCE(MapXLinkage)
        protected:
            Real m_beta = 0.1;
            bool m_initialized = false;
        public:
            void addInputParameters();
            void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) override;
            void transfer(std::vector<Real>& x, const std::vector<Real>& var) override;
        };

    }
}
#endif