#include "one_peak_s10.h"

namespace ofec::free_peaks {
	void OnePeakS10::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		OnePeakBase::initialize(pro, subspace_name, param);
	}

	Real OnePeakS10::evaluate_(Real dummy, size_t var_size) {
		return m_height - pow(m_height, -2) * pow(dummy, 3);
	}
}