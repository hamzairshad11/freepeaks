#include "map_x_deceptive.h"
#include "../../../free_peaks.h"
#include <cmath>
#include <algorithm>

namespace ofec {
    namespace free_peaks {

        void MapXDeceptive::addInputParameters() {
            m_input_parameters.add("alpha", new InputParameterValueTypeReal(m_alpha));
        }

        void MapXDeceptive::initialize(Problem* pro,
            const std::string& subspace_name,
            const ParameterMap& param)
        {
            TransformBase::initialize(pro, subspace_name, param);
            m_initialized = true;
        }

        void MapXDeceptive::transfer(std::vector<Real>& x,
            const std::vector<Real>& /*var*/)
        {
            if (!m_initialized || x.empty()) return;
            for (size_t i = 0; i < x.size(); ++i) {
                double xi = static_cast<double>(x[i]);
                double yi = xi + static_cast<double>(m_alpha) * std::sin(2.0 * 3.14159265358979323846 * xi);
                x[i] = static_cast<Real>(std::max(0.0, std::min(1.0, yi)));
            }
        }
    } // namespace free_peaks
    
} // namespace ofec
