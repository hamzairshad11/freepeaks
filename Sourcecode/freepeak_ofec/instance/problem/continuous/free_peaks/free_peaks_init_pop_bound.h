/********* Begin Register Information **********
{
	"name": "free_peaks_init_pop_bounded",
	"identifier": "FreePeaksInitPopBounded",
	"tags": [ "continuous"]
}
*********** End Register Information **********/

#ifndef OFEC_FREE_PEAKS_INIT_POP_BOUNDED_H
#define OFEC_FREE_PEAKS_INIT_POP_BOUNDED_H

#include "free_peaks.h"
#include "../init_pop_bounded.h"

namespace ofec {
	class FreePeaksInitPopBounded : public FreePeaks, public InitPopBounded {
		OFEC_CONCRETE_INSTANCE(FreePeaksInitPopBounded)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		size_t m_id_peak_to_init;
	};
}

#endif // !OFEC_FREE_PEAKS_INIT_POP_BOUNDED_H
