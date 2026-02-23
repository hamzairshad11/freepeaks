#ifndef OFEC_FREE_PEAKS_ONE_PEAK_S8_H
#define OFEC_FREE_PEAKS_ONE_PEAK_S8_H

#include "one_peak_base.h"

namespace ofec::free_peaks {
	class OnePeakS8 final : public OnePeakBase {
		OFEC_CONCRETE_INSTANCE(OnePeakS8)
	private:
		Real m_radius, m_eta;
		Real m_r_ratio;

		int m_m;



	public:
		void addInputParameters();
		virtual void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) override;
		Real evaluate_(Real dummy, size_t var_size) override;
		virtual void bindData()override;
	};
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_S8_H