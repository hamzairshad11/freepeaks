#include "register_one_peak.h"
#include "../../../factory.h"
#include "one_peak_s1.h"
#include "one_peak_s2.h"
#include "one_peak_s3.h"
#include "one_peak_s4.h"
#include "one_peak_s5.h"
#include "one_peak_s6.h"
#include "one_peak_s7.h"
#include "one_peak_s8.h"
#include "one_peak_s9.h"
#include "one_peak_s10.h"

//#include "one_peak_map.h"

namespace ofec::free_peaks {
	void registerOnePeak() {
		REGISTER_FP(OnePeakBase, OnePeakS1, "s1");
		REGISTER_FP(OnePeakBase, OnePeakS2, "s2");
		REGISTER_FP(OnePeakBase, OnePeakS3, "s3");
		REGISTER_FP(OnePeakBase, OnePeakS4, "s4");
		REGISTER_FP(OnePeakBase, OnePeakS5, "s5");
		REGISTER_FP(OnePeakBase, OnePeakS6, "s6");
		REGISTER_FP(OnePeakBase, OnePeakS7, "s7");
		REGISTER_FP(OnePeakBase, OnePeakS8, "s8");
		REGISTER_FP(OnePeakBase, OnePeakS9, "s9");
		REGISTER_FP(OnePeakBase, OnePeakS10, "s10");
		//REGISTER_FP(OnePeakBBOB_Rug, OnePeakBBOB_Rug, "sBBOB");
	}
}