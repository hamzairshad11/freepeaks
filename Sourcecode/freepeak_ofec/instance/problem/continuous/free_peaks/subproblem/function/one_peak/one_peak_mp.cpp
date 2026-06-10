#include "one_peak_mp.h"
#include <cmath>

namespace ofec::free_peaks {
    void OnePeakMP::addInputParameters() {
        m_input_parameters.add("radius", new RangedReal(m_radius, 1e-6, 1e6, 32));
    }

    Real OnePeakMP::evaluate_(Real dummy, size_t) {
        const Real r = std::max<Real>(m_radius, 1e-6);
        return m_height * std::exp(-0.5 * (dummy / r) * (dummy / r));
    }
}
