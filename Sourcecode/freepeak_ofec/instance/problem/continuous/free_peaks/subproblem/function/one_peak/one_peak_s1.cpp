#include "one_peak_s1.h"
#include "../../../free_peaks.h"

namespace ofec::free_peaks {
	OnePeakS1::OnePeakS1(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		OnePeakBase(pro, subspace_name, param) {}

	Real OnePeakS1::evaluate_(Real dummy, size_t var_size) {
		return m_height - dummy;
	}
}