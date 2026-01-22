#include "one_peak_s9.h"

namespace ofec::free_peaks {
	OnePeakS9::OnePeakS9(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		OnePeakBase(pro, subspace_name, param) {}

	Real OnePeakS9::evaluate_(Real dummy, size_t var_size) {
		return m_height - pow(m_height, -0.4) * pow(dummy, 1.4);
	}
}