/********* Begin Register Information **********
{
	"name": "free_peaks_init_pop_imbalance",
	"identifier": "FreePeaksInitPopImbalance",
	"tags": [ "continuous"]
}
*********** End Register Information **********/

#ifndef OFEC_FREE_PEAKS_INIT_POP_IMBALANCE_H
#define OFEC_FREE_PEAKS_INIT_POP_IMBALANCE_H

#include "free_peaks.h"

namespace ofec {
	class FreePeaksInitPopImbalance : public FreePeaks {
		OFEC_CONCRETE_INSTANCE(FreePeaksInitPopImbalance)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void initializeVariables(VariableBase &vars, Random *rnd) const override;
		std::string m_prob_ratio_init;
		std::vector<Real> m_each_prob;
		Real m_sum_prob;
	};
}

#endif // !OFEC_FREE_PEAKS_INIT_POP_IMBALANCE_H
