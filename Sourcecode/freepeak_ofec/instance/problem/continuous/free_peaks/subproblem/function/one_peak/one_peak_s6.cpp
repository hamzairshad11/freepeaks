#include "one_peak_s6.h"

namespace ofec::free_peaks {
	void OnePeakS6::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		OnePeakBase::initialize(pro, subspace_name, param);
	}

	Real OnePeakS6::evaluate_(Real dummy, size_t var_size) {
		return m_height - exp(2 * sqrt(dummy / sqrt(var_size))) + 1;
	}
}