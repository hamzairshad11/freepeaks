#ifndef OFEC_FREE_PEAKS_ONE_PEAK_S7_H
#define OFEC_FREE_PEAKS_ONE_PEAK_S7_H

#include "one_peak_base.h"

namespace ofec::free_peaks {
	class OnePeakS7 final : public OnePeakBase {
	private:
		Real m_radius;

	public:
		OnePeakS7(Problem *pro, const std::string &subspace_name, const ParameterMap &param);
		Real evaluate_(Real dummy, size_t var_size) override;
		void bindData() override;
	};
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_S7_H
