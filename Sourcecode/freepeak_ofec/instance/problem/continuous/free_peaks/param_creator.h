#ifndef OFEC_FREE_PEAKS_PARAM_CREATOR_H
#define OFEC_FREE_PEAKS_PARAM_CREATOR_H

#include "../../../../core/parameter/parameter_map.h"

namespace ofec::free_peaks {
	struct ParamIO {
		static void writeTreeData(std::ofstream &out, const std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>> &tree_name);
		static void writeNameRatio(std::ofstream &out, const std::vector<std::pair<std::string, double>> &data);
		static void writeParameter(const ParameterMap &v, std::ofstream &out);

		static void readTreeData(std::ifstream &in, std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>> &tree_name);
		static void readNameRatio(std::ifstream &in, std::vector<std::pair<std::string, double>> &data);
		static void readParameter(std::ifstream &in, ParameterMap &v);
		static void readSizes(std::ifstream &in, size_t &num_vars, size_t &num_objs, size_t &num_cons);
	};
}

#endif // !OFEC_FREE_PEAKS_PARAM_CREATOR_H
