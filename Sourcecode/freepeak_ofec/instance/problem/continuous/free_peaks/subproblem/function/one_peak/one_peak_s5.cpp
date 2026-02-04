#include "one_peak_s5.h"

namespace ofec::free_peaks {
	void OnePeakS5::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		OnePeakBase::initialize(pro, subspace_name, param);
	}

	Real OnePeakS5::evaluate_(Real dummy, size_t var_size) {
		return m_height - pow(dummy, 2) / m_height;
	}
}