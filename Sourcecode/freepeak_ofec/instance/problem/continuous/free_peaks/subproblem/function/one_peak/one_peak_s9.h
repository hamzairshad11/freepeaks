#ifndef OFEC_FREE_PEAKS_ONE_PEAK_S9_H
#define OFEC_FREE_PEAKS_ONE_PEAK_S9_H

#include "one_peak_base.h"

namespace ofec::free_peaks {
	class OnePeakS9 final : public OnePeakBase {
	public:
		OnePeakS9(Problem *pro, const std::string &subspace_name, const ParameterMap &param);
		Real evaluate_(Real dummy, size_t var_size) override;
	};
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_S9_H
