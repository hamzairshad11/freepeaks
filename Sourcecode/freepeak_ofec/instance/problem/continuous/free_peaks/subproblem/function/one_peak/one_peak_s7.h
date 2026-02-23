#ifndef OFEC_FREE_PEAKS_ONE_PEAK_S7_H
#define OFEC_FREE_PEAKS_ONE_PEAK_S7_H

#include "one_peak_base.h"

namespace ofec::free_peaks {
	class OnePeakS7 final : public OnePeakBase {
		OFEC_CONCRETE_INSTANCE(OnePeakS7)
	private:
		Real m_radius;
		Real m_r_ratio;

	public:
		void addInputParameters();
		virtual void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) override;
		Real evaluate_(Real dummy, size_t var_size) override;
		virtual void bindData()	override;
	};
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_S7_H