#ifndef OFEC_FREE_PEAKS_PARAM_FUN_H
#define OFEC_FREE_PEAKS_PARAM_FUN_H

#include <vector>
#include <string>



namespace ofec::free_peaks {

	extern void changeParameter(double alpha);
	extern void changeParameter(double alpha, const std::string&filename);
	extern void changeParameterRatio(
		std::string filename, 
		std::vector<std::pair<double, std::string>>& ratio_name);
	

}

#endif