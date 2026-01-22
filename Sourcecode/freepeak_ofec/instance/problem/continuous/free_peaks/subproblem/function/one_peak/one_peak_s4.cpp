#include "one_peak_s4.h"

namespace ofec::free_peaks {
	OnePeakS4::OnePeakS4(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		OnePeakBase(pro, subspace_name, param) {}

	Real OnePeakS4::evaluate_(Real dummy, size_t var_size) {
		return m_height / (dummy + 1);
	}
}