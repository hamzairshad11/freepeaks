#ifndef OFEC_FREE_PEAKS_ONE_PEAK_S1_H
#define OFEC_FREE_PEAKS_ONE_PEAK_S1_H

#include "one_peak_base.h"

namespace ofec::free_peaks {
	class OnePeakS1 final : public OnePeakBase {

		OFEC_CONCRETE_INSTANCE(OnePeakS1)
	public:

		void addInputParameters() {}
		virtual void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param)override;
		Real evaluate_(Real dummy, size_t var_size) override;
	};
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_S1_H