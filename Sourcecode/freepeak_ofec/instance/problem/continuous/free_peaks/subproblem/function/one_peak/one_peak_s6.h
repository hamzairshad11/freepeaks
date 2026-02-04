#ifndef OFEC_FREE_PEAKS_ONE_PEAK_S6_H
#define OFEC_FREE_PEAKS_ONE_PEAK_S6_H

#include "one_peak_base.h"

namespace ofec::free_peaks {
	class OnePeakS6 final : public OnePeakBase {
		OFEC_CONCRETE_INSTANCE(OnePeakS6)
	public:
		void addInputParameters() {}	
		virtual void initialize(Problem *pro, const std::string &subspace_name, const ParameterMap &param) override;
		Real evaluate_(Real dummy, size_t var_size) override;
	};
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_S6_H
