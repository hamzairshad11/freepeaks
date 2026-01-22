#include "param_creator.h"
#include "../../../../utility/functional.h"

namespace ofec::free_peaks {
	void ParamIO::writeTreeData(std::ofstream &out, const std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>> &tree_name) {
		out << "subtree_number \t" << tree_name.size() << std::endl;
		for (auto &it : tree_name) {

			out << "subtree_name \t" << it.first << std::endl;
			writeNameRatio(out, it.second);
		}
	}

	void ParamIO::writeNameRatio(std::ofstream &out, const std::vector<std::pair<std::string, double>> &data) {
		out << "num\t" << data.size() << std::endl;
		for (auto &it : data) {
			out << "[ " << it.first << " , " << it.second << " ]" << "\t";
		}
		out << std::endl;
	}

	void ParamIO::readTreeData(std::ifstream &in, std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>> &tree_name) {
		std::string tmp;
		int size(0);
		in >> tmp >> size;
		tree_name.resize(size);
		for (auto &it : tree_name) {
			in >> tmp >> it.first;
			readNameRatio(in, it.second);
		}
	}

	void ParamIO::readNameRatio(std::ifstream &in, std::vector<std::pair<std::string, double>> &data) {
		std::string tmp;
		int size(0);
		in >> tmp >> size;
		data.resize(size);
		for (auto &it : data) {
			in >> tmp >> it.first >> tmp >> it.second >> tmp;
		}
	}

	void ParamIO::writeParameter(const ParameterMap &v, std::ofstream &out) {
		std::string svalue;
		out << "_%BEGIN%" << std::endl;
		for (auto it : v) {
			if (it.first.find(char(2)) != std::string::npos) {
				continue;
			}
			// TODO
			//out << it.first << "\t" << it.second << std::endl;
			if (static_cast<ParameterType>(it.second.index()) == ParameterType::kString) {
				svalue = std::get<std::string>(it.second);
				if (svalue == "_%vd%") {
					std::vector<double> tmp(paramToVec<double>(v, it.first));
					outputFileVec(out, tmp);
				}
				else if (svalue == "_%vt%") {
					std::vector<int> tmp(paramToVec<int>(v, it.first));
					outputFileVec(out, tmp);
				}
				else if (svalue == "_%vs%") {
					std::vector<std::string> tmp(paramToVec<std::string>(v, it.first));
					outputFileVec(out, tmp);
				}
				else if (svalue == "_%vpdd%") {
					outputFileVecPair(out, paramToVecPair<double, double>(v, it.first));
				}
				else if (svalue == "_%vpii%") {
					outputFileVecPair(out, paramToVecPair<int, int>(v, it.first));
				}
				else if (svalue == "_%vpss%") {
					outputFileVecPair(out, paramToVecPair<std::string, std::string>(v, it.first));
				}
			}
		}
		out << "_%END%" << std::endl;
	}

	void ParamIO::readParameter(std::ifstream &in, ParameterMap &v) {
		std::string name;
		int ival(0);
		double dvalue(0);
		std::string svalue;
		std::string letter;
		for (auto i = 'a'; i <= 'z'; i++) {
			letter += i;
			letter += i - 32;
		}
		letter += "\\/:";
		while (in) {
			in >> name;
			if (name == "_%BEGIN%")continue;
			if (name == "_%END%")break;
			in >> svalue;
			if (svalue.compare("true") == 0) {
				v[name] = true;
			}
			else if (svalue.compare("false") == 0) {
				v[name] = false;
			}
			else if (svalue.find_first_of(letter) != std::string::npos) {
				if (svalue.size() == 1)
					v[name] = svalue[0];
				else {
					v[name] = svalue;
					if (svalue == "_%vd%") {
						std::vector<double> tmp(inputFileVec<double>(in));
						vecToParam(v, name, tmp);
					}
					else if (svalue == "_%vt%") {
						std::vector<int> tmp(inputFileVec<int>(in));
						vecToParam(v, name, tmp);
					}
					else if (svalue == "_%vs%") {
						std::vector<std::string> tmp(inputFileVec<std::string>(in));
						vecToParam(v, name, tmp);
					}
					else if (svalue == "_%vpdd%") {
						std::vector<std::pair<double, double>> tmp(inputFileVecPair<double, double>(in));
						vecPairToParam(v, name, tmp);
					}
					else if (svalue == "_%vpii%") {
						std::vector<std::pair<int, int>> tmp(inputFileVecPair<int, int>(in));
						vecPairToParam(v, name, tmp);
					}
					else if (svalue == "_%vpss%") {
						std::vector<std::pair<std::string, std::string>> tmp(inputFileVecPair<std::string, std::string>(in));
						vecPairToParam(v, name, tmp);
					}
				}
			}
			else if (svalue.find_first_of('.') != std::string::npos) {
				double val = atof(svalue.c_str());
				v[name] = val;
			}
			else {
				int val = atoi(svalue.c_str());
				v[name] = val;
			}

		}
	}

	void ParamIO::readSizes(std::ifstream &in, size_t &num_vars, size_t &num_objs, size_t &num_cons) {
		std::string name;
		std::string svalue;
		for (size_t i = 0; i < 3; ++i) {
			in >> name;
			if (name == "number_of_variables") {
				in >> svalue;
				num_vars = std::stoull(svalue);
			}
			if (name == "number_of_objectives") {
				in >> svalue;
				num_objs = std::stoull(svalue);
			}
			if (name == "number_of_constraints") {
				in >> svalue;
				num_cons = std::stoull(svalue);
			}
		}
	}
}
