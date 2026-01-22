#include "register_function.h"
#include "../../factory.h"
#include "one_peak_function.h"

namespace ofec::free_peaks {
	void registerFunction() {
		REGISTER_FP(FunctionBase, OnePeakFunction, "one_peak");
	}
}