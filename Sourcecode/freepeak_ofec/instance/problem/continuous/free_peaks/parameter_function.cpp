#include "parameter_function.h"
#include "parameter_function.h"
#include "param_creator.h"
#include "../../../../core/global.h"

void  ofec::free_peaks::changeParameter(double alpha) {
	using namespace ofec;



	std::string file_name = "sop/3_s5_1_s2.txt";
	std::string file_path = ofec::g_working_directory + "/instance/problem/continuous/free_peaks/" + file_name;


	std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>  treeData;
	std::ifstream in(file_path);
	free_peaks::ParamIO::readTreeData(in, treeData);


	std::vector<ParameterMap> v_subpro_param;

	int number_objectives = -1, num_cons = -1;
	for (size_t id_subpro = 0; id_subpro < treeData.front().second.size(); ++id_subpro) {
		ParameterMap cur_param;
		free_peaks::ParamIO::readParameter(in, cur_param);
		v_subpro_param.push_back(cur_param);
	}



	std::vector<double> values(treeData.front().second.size(), 1);

	if (values.size() > 1) {
		for (int idx(0); idx < values.size() - 1; ++idx) {
			values[idx] = (1.0 - alpha) / (values.size() - 1);
		}
		values.back() = alpha;
	}

	{
		std::string file_name = "sop/3_s5_s2.txt";
		std::string file_path = ofec::g_working_directory + "/instance/problem/continuous/free_peaks/" + file_name;


		for (int idx(0); idx < values.size(); ++idx) {
			treeData.front().second[idx].second = values[idx];
		}

		//std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>  treeData;
		std::ofstream out(file_path);
		free_peaks::ParamIO::writeTreeData(out, treeData);






		//std::vector<ParameterMap> v_subpro_param;

		int number_objectives = -1, num_cons = -1;
		for (size_t id_subpro = 0; id_subpro < treeData.front().second.size(); ++id_subpro) {
			auto& cur_param = v_subpro_param[id_subpro];
			free_peaks::ParamIO::writeParameter(cur_param, out);

		}
	}


	/*{
		std::string file_name = "sop/3_s5_s2.txt";
		std::string file_path = ofec::g_working_directory + "/instance/problem/continuous/free_peaks/" + file_name;


		std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>  treeData2;
		std::ifstream in(file_path);
		free_peaks::ParamIO::readTreeData(in, treeData2);


		std::vector<ParameterMap> v_subpro_param2;

		int number_objectives = -1, num_cons = -1;
		for (size_t id_subpro = 0; id_subpro < treeData2.front().second.size(); ++id_subpro) {
			ParameterMap cur_param;
			free_peaks::ParamIO::readParameter(in, cur_param);
			v_subpro_param2.push_back(cur_param);
		}
	}*/
}


void changeParameter(double alpha, const std::string& filename) {
	using namespace ofec;



	std::string file_name = "sop/3_s5_1_s2.txt";
	std::string file_path = ofec::g_working_directory + "/instance/problem/continuous/free_peaks/" + file_name;


	std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>  treeData;
	std::ifstream in(file_path);
	free_peaks::ParamIO::readTreeData(in, treeData);


	std::vector<ParameterMap> v_subpro_param;

	int number_objectives = -1, num_cons = -1;
	for (size_t id_subpro = 0; id_subpro < treeData.front().second.size(); ++id_subpro) {
		ParameterMap cur_param;
		free_peaks::ParamIO::readParameter(in, cur_param);
		v_subpro_param.push_back(cur_param);
	}



	std::vector<double> values(treeData.front().second.size(), 1);

	if (values.size() > 1) {
		for (int idx(0); idx < values.size() - 1; ++idx) {
			values[idx] = (1.0 - alpha) / (values.size() - 1);
		}
		values.back() = alpha;
	}

	{
		std::string file_name = filename;
		std::string file_path = ofec::g_working_directory + "/instance/problem/continuous/free_peaks/" + file_name;


		for (int idx(0); idx < values.size(); ++idx) {
			treeData.front().second[idx].second = values[idx];
		}

		//std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>  treeData;
		std::ofstream out(file_path);
		free_peaks::ParamIO::writeTreeData(out, treeData);

		//std::vector<ParameterMap> v_subpro_param;

		int number_objectives = -1, num_cons = -1;
		for (size_t id_subpro = 0; id_subpro < treeData.front().second.size(); ++id_subpro) {
			auto& cur_param = v_subpro_param[id_subpro];
			free_peaks::ParamIO::writeParameter(cur_param, out);

		}
	}

}

void ofec::free_peaks::changeParameterRatio(std::string filename, 
	std::vector<std::pair<double, std::string>>& ratio_name)
{

	std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>  treeData(1);

	treeData.front().first = "root";
	treeData.front().second.resize(ratio_name.size());
	
	std::vector<std::string> subspace_names(ratio_name.size());
	for (int idx(0); idx < subspace_names.size(); ++idx) {
		subspace_names[idx] = "s" + std::to_string(idx);
		treeData.front().second[idx].first = subspace_names[idx];
		treeData.front().second[idx].second = ratio_name[idx].first;
	}


	std::vector<ParameterMap> v_subpro_param;

	for (size_t id_subpro = 0; id_subpro < treeData.front().second.size(); ++id_subpro) {
		ParameterMap curPara;
		
		curPara["subspace"] = subspace_names[id_subpro];
		curPara["generation"] = std::string("read_file");
		curPara["file"] = "sop/" + ratio_name[id_subpro].second;

		v_subpro_param.push_back(curPara);
	}

	


	std::string file_name = filename + ".txt";
	std::string file_path = ofec::g_working_directory + "/instance/problem/continuous/free_peaks/" + file_name;

	//std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>  treeData;
	std::ofstream out(file_path);
	free_peaks::ParamIO::writeTreeData(out, treeData);


	for (size_t id_subpro = 0; id_subpro < treeData.front().second.size(); ++id_subpro) {
		auto& cur_param = v_subpro_param[id_subpro];
		free_peaks::ParamIO::writeParameter(cur_param, out);

	}


}
