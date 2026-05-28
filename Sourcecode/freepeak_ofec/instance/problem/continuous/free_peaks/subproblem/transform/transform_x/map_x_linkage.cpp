#include "map_x_linkage.h"
#include "../../../free_peaks.h"
#include <cmath>
#include <algorithm>

namespace ofec {
    namespace free_peaks {

        void MapXLinkage::addInputParameters() {
            m_input_parameters.add("beta", new InputParameterValueTypeReal(m_beta));
        }

        void MapXLinkage::initialize(Problem* pro,
            const std::string& subspace_name,
            const ParameterMap& param)
        {
            TransformBase::initialize(pro, subspace_name, param);
            m_initialized = true;
        }

        void MapXLinkage::transfer(std::vector<Real>& x,
            const std::vector<Real>& /*var*/)
        {
            if (!m_initialized || x.size() < 2) return;
            std::vector<Real> x_orig = x;
            const double scaled_beta = static_cast<double>(m_beta) / 40000.0;
            for (size_t i = 0; i < x.size(); ++i) {
                size_t j = (i + 1) % x.size();
                double xi = static_cast<double>(x_orig[i]);
                double xj = static_cast<double>(x_orig[j]);
                x[i] = static_cast<Real>(xi + scaled_beta * xj * xj);
                // NO clamping — residuals must remain in function domain
            }
        }

    } // namespace free_peaks
} // namespace ofec