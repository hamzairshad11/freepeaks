#ifndef OFEC_FREE_PEAKS_ONE_PEAK_S2_H
#define OFEC_FREE_PEAKS_ONE_PEAK_S2_H

#include "one_peak_base.h"

namespace ofec::free_peaks {
	class OnePeakS2 final : public OnePeakBase {
	public:
		OnePeakS2(Problem *pro, const std::string &subspace_name, const ParameterMap &param);
		Real evaluate_(Real dummy, size_t var_size) override;
	};
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_S2_H
