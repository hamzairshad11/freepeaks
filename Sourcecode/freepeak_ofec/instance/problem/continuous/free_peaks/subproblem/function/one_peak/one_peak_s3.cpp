#include "one_peak_s3.h"

namespace ofec::free_peaks {
	void OnePeakS3::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		OnePeakBase::initialize(pro, subspace_name, param);
	}


	Real OnePeakS3::evaluate_(Real dummy, size_t var_size) {
		return m_height - sqrt(m_height * dummy);
	}
}