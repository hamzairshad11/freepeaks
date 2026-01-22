#include "one_peak_s5.h"

namespace ofec::free_peaks {
	OnePeakS5::OnePeakS5(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		OnePeakBase(pro, subspace_name, param) {}

	Real OnePeakS5::evaluate_(Real dummy, size_t var_size) {
		return m_height - pow(dummy, 2) / m_height;
	}
}