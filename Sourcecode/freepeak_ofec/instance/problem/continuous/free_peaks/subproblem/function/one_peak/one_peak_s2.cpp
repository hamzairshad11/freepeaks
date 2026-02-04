#include "one_peak_s2.h"

namespace ofec::free_peaks {
	void OnePeakS2::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		OnePeakBase::initialize(pro, subspace_name, param);
	}

	Real OnePeakS2::evaluate_(Real dummy, size_t var_size) {
		return m_height * exp(-dummy);
	}
}