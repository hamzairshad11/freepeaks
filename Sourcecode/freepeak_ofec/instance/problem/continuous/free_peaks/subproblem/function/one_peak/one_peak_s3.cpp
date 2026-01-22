#include "one_peak_s3.h"

namespace ofec::free_peaks {
	OnePeakS3::OnePeakS3(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		OnePeakBase(pro, subspace_name, param) {}

	Real OnePeakS3::evaluate_(Real dummy, size_t var_size) {
		return m_height - sqrt(m_height * dummy);
	}
}