#include "map_x_asymmetry.h"
#include "../../../free_peaks.h"
#include <cmath>

namespace ofec::free_peaks {

    void MapXAsymmetry::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
        MapXBase::initialize(pro, subspace_name, param);

        if (param.has("beta")) {
            m_beta = std::get<Real>(param.at("beta"));
        }
    }

    void MapXAsymmetry::transfer(std::vector<Real>& x_, const std::vector<Real>& var) {
        if (x_.size() <= 1) return;

        const size_t n = x_.size();
        for (size_t i = 0; i < n; ++i) {
            if (x_[i] > Real(0.0)) {
                // Cast to Real so the coefficient scales properly across all variables.
                Real idx_ratio = static_cast<Real>(i) / static_cast<Real>(n - 1);
                Real exponent = Real(1.0) + m_beta * idx_ratio * std::log(x_[i] + Real(1.0));
                x_[i] = std::pow(x_[i], exponent);
            }
        }
    }

} // namespace ofec::free_peaks