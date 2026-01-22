#include "free_peaks.h"
#include "subproblem/distance/register_distance.h"
#include "subproblem/function/register_function.h"
#include "subproblem/function/one_peak/register_one_peak.h"
#include "subproblem/constraint/register_constraint.h"
#include "subproblem/transform/register_transform.h"
#include "../../../../core/global.h"
#include "param_creator.h"
#include "../../../../core/problem/solution.h"
#include "../../../../utility/nondominated_sorting/filter_sort.h"
#include "../../../../utility/functional.h"

namespace ofec {
	bool FreePeaks::ms_registered = false;

	void FreePeaks::addInputParameters() {
		std::list<std::string> group = { "read_file", "random", "default" };
		m_input_parameters.add("generation type", new EnumeratedString(m_generation_type, group, "read_file"));
		m_input_parameters.add("dataFile1", new FileName(m_file_name,
			"instance/problem/continuous/free_peaks",
			"config (*.txt)", "sop/2d_s1_to_s8.txt"));	
	}

	void FreePeaks::initialize_(Environment *env) {
		Continuous::initialize_(env);
		if (!ms_registered) {
			free_peaks::registerDistance();
			free_peaks::registerFunction();
			free_peaks::registerOnePeak();
			free_peaks::registerConstraint();
			free_peaks::registerTransform();
			ms_registered = true;
		}
		if (m_generation_type == "read_file") {
			readFromFile(g_working_directory + "/instance/problem/continuous/free_peaks/" + m_file_name);
		}
		else if (m_generation_type == "random") {

		}
		else {

		}
		for (auto &it : m_subspace_tree.name_box_subproblem) {
			if (it.second.second == nullptr)
				continue;
			it.second.second->bindData();
		}
	}

	void FreePeaks::evaluate(const VariableBase &vars, std::vector<Real> &objs, std::vector<Real> &cons) const {
		auto &x = dynamic_cast<const VariableVector<>&>(vars);
		std::vector<Real> vars_(x.begin(), x.end()); //for parallel running
		auto name = m_subspace_tree.tree->getRegionName(vars_);
		m_subspace_tree.name_box_subproblem.at(name).second->evaluate(vars_, objs, cons);
	}

	void FreePeaks::updateOptima(Environment *env) {
		auto new_optima = new Optima<>();
		if (m_number_objectives == 1) {
			for (auto &it : m_subspace_tree.name_box_subproblem) {
				if (it.second.second == nullptr)
					continue;
				auto &subpro_optima = it.second.second->optima();
				for (size_t i = 0; i < subpro_optima.numberSolutions(); ++i) {
					new_optima->appendSolution(subpro_optima.solution(i));
				}
			}
			m_objective_accuracy = 1e-3;
		}
		else if(m_number_objectives > 1) {
			//依次在每个子空间的峰之间进行采样，再进行排序
			for (auto& it : m_subspace_tree.name_box_subproblem) {
				if (it.second.second == nullptr)
					continue;
				auto& box = m_subspace_tree.tree->getBox(it.first);
				auto& optima = it.second.second->optima();//optima的个数为：num_obj*num_subspace
				//size_t number_objectives = m_param->get<int>("number of objectives");
				size_t number_objectives = m_number_objectives;
				////在一组点内逐个加入进行线采样
				std::vector<std::vector<Real>> subspace_peak_pos;
				for (size_t j = 0; j < number_objectives; ++j) {
					subspace_peak_pos.emplace_back(optima.solution(j).variable().vector());
				}
				std::vector<std::vector<Real>> sampling_sols;
				sampleAmongPoints(subspace_peak_pos, sampling_sols);
			}
			////依次在每个子空间的峰之间进行采样，再进行排序
			//for (auto& it : m_subspace_tree.name_box_subproblem) {
			//	if (it.second.second == nullptr)
			//		continue;
			//	auto& box = m_subspace_tree.tree->getBox(it.first);
			//	auto& optima = it.second.second->optima();//optima的个数为：num_obj*num_subspace
			//	size_t number_objectives = m_param->get<int>("number of objectives");
			//	////在一组点内逐个加入进行线采样
			//	std::vector<std::vector<Real>> subspace_peak_pos;
			//	for (size_t j = 0; j < number_objectives; ++j) {
			//		subspace_peak_pos.emplace_back(optima.variable(j).vect());
			//	}
			//	std::vector<std::vector<Real>> sampling_sols;
			//	sampleAmongPoints(subspace_peak_pos, sampling_sols);

			//	//对samplings评价，局部排序，加入前排
			//	std::vector<std::vector<Real>> temp_objs;
			//	for (size_t i = 0; i < sampling_sols.size(); ++i) {
			//		Solution<> sol(number_objectives, 0);
			//		sol.variable() = sampling_sols[i];
			//		evaluate_(sol);
			//		temp_objs.emplace_back(sol.objective());
			//	}

			//	//排序
			//	std::vector<std::vector<Real>*> objs;
			//	for (size_t i = 0; i < temp_objs.size(); ++i) {
			//		objs.emplace_back(&temp_objs[i]);
			//	}
			//	std::vector<int> rank;
			//	ofec::nd_sort::filterSortP<Real>(objs, rank, m_optimize_mode);
			//	for (size_t i = 0; i < rank.size(); ++i) {
			//		if (rank[i] == 0) {
			//			optima()->appendVar(sampling_sols[i]);
			//			m_optima->appendObj(temp_objs[i]);
			//		}
			//	}
			//}
			//m_optima->setVariableGiven(true);
			//m_optima->setObjectiveGiven(true);
		}
		m_optima.reset(new_optima);
	}

	void FreePeaks::readFromFile(const std::string &file_path) {
		std::ifstream in(file_path);
		if (in.fail()) {
			throw Exception("Free Peaks config file not exist.");
		}

		size_t number_variables, number_objectives, number_constraints;
		free_peaks::ParamIO::readSizes(in, number_variables, number_objectives, number_constraints);
		resizeVariable(number_variables);
		for (size_t i = 0; i < number_variables; ++i) {
			m_domain.setRange(0.0, 1.0, i);
		}
		resizeObjective(number_objectives);
		for (size_t i = 0; i < number_objectives; ++i) {
			m_optimize_mode[i] = OptimizeMode::kMaximize;
		}
		resizeConstraint(number_constraints);

		m_subspace_tree.clear();
		std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>  treeData;
		free_peaks::ParamIO::readTreeData(in, treeData);
		std::vector<std::pair<Real, Real>> ranges(m_number_variables);
		for (size_t i = 0; i < m_number_variables; ++i) {
			ranges[i] = m_domain[i].limit;
		}
		m_subspace_tree.createTree(ranges, treeData);
		ParameterMap subpro_param;
		for (size_t id_subpro = 0; id_subpro < m_subspace_tree.tree->size(); ++id_subpro) {
			free_peaks::ParamIO::readParameter(in, subpro_param);
			auto new_subpro = new free_peaks::Subproblem(subpro_param, this);
			std::string subspace_name;
			if (static_cast<ParameterType>(subpro_param.at("subspace").index()) == ParameterType::kChar) {
				subspace_name.push_back(subpro_param.get<char>("subspace"));
			}
			else if (static_cast<ParameterType>(subpro_param.at("subspace").index()) == ParameterType::kString) {
				subspace_name = subpro_param.get<std::string>("subspace");
			}
			m_subspace_tree.name_box_subproblem.at(subspace_name).second.reset(new_subpro);
		}	
	}

	void FreePeaks::updateCandidates(const SolutionBase &sol, std::vector<std::unique_ptr<SolutionBase>> &candidates) const {
		if (m_number_objectives == 1) {
			if (candidates.empty()) {
				for (size_t i = 0; i < m_subspace_tree.tree->size(); i++) {
					candidates.emplace_back(new Solution<>(dynamic_cast<const Solution<>&>(sol)));
				}
			}
			else {
				size_t id_peak = m_subspace_tree.tree->getRegionIdx(
					dynamic_cast<const Solution<>&>(sol).variable().vector());
				if (sol.objective(0) > candidates[id_peak]->objective(0)) {
					dynamic_cast<Solution<>&>(*candidates[id_peak]) = dynamic_cast<const Solution<>&>(sol);
				}
			}
		}
	}

	bool FreePeaks::isSolved(const std::vector<std::unique_ptr<SolutionBase>> &candidates) const {
		bool result = false;
		if (m_number_objectives == 1) {
			std::list<size_t> global_optima = { 0 };
			Real best_obj = optima()->solution(0).objective(0);
			for (size_t i = 1; i < optima()->numberSolutions(); ++i) {
				if (optima()->solution(i).objective(0) < best_obj) {
					global_optima = { i };
					best_obj = optima()->solution(i).objective(0);
				}
				else if (optima()->solution(i).objective(0) == best_obj) {
					global_optima.push_back(i);
				}
				
			}
			result = true;
			for (size_t i : global_optima) {
				if (candidates[i]->objectiveDistance(optima()->solution(i)) > m_objective_accuracy) {
					result = false;
					break;
				}
			}
		}
		return result;
	}

	std::vector<Real> FreePeaks::errorToPeaks(const std::vector<std::unique_ptr<SolutionBase>> &candidates) const {
		std::vector<Real> result(candidates.size());
		if (m_number_objectives == 1) {
			for (size_t i = 0; i < candidates.size(); i++) {
				result[i] = candidates[i]->objectiveDistance(optima()->solution(i));
			}
		}
		return result;
	}

	void FreePeaks::writeToFile(const std::string &file_path) const {
		std::ofstream out(file_path);
		free_peaks::ParamIO::writeTreeData(out, m_subspace_tree.treeData());
		out << "\n\n";
		for (auto &it : m_subspace_tree.name_box_subproblem) {
			if (it.second.second != nullptr) {
				free_peaks::ParamIO::writeParameter(it.second.second->parameter(), out);
				out << '\n';
			}
		}
	}
	
	free_peaks::Subproblem* FreePeaks::subproblem(const std::string &subspace_name) const {
		return m_subspace_tree.name_box_subproblem.at(subspace_name).second.get();
	}
}