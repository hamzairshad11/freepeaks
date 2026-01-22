#include "subproblem.h"
#include <fstream>
#include "../factory.h"
#include "../../../../../core/global.h"
#include "../param_creator.h"

namespace ofec::free_peaks {
	Subproblem::Subproblem(const ParameterMap &param, Problem *pro) :
		m_param(param), m_pro(pro)
	{
		auto generation = m_param.get<std::string>("generation");
		if (generation == "read_file") {
			auto file_name = m_param.get<std::string>("file") + ".txt";
			readFromFile(g_working_directory + "/instance/problem/continuous/free_peaks/subproblem/" + file_name);
		}
		else if (generation == "random") {

		}
		else {

		}
	}

	void Subproblem::setOptima() {
		Solution<> opt_sol(m_pro->numberObjectives(), m_pro->numberConstraints(), m_pro->numberVariables());
		if (m_pro->numberObjectives() == 1) {
			m_function->mapOptimalVars();
			for (auto &opt_var : m_function->optimalVars()) {
				opt_sol.variable() = opt_var;
				evaluate(opt_sol.variable().vector(), opt_sol.objective(), opt_sol.constraint());
				m_optima.appendSolution(opt_sol);
			}
		}
		else if(m_pro->numberObjectives()>1){
			m_function->mapOptimalVars();
			for (auto&opt_var : m_function->optimalVars()) {
				opt_sol.variable() = opt_var;
				evaluate(opt_sol.variable().vector(), opt_sol.objective(), opt_sol.constraint());
				m_optima.appendSolution(opt_sol);
			}
		}
	}

	void Subproblem::bindData() {
		m_distance->bindData();
		m_function->bindData();
		for (auto& it : m_constraints) {
			it->bindData();
		}
		for (auto& it : m_trans_vars) {
			it->bindData();
			it->bindSubPro(this);
		}
		for (auto& it : m_trans_objs) {
			it->bindData();
			it->bindSubPro(this);
		}
		setOptima();
	}

	void Subproblem::evaluate(const std::vector<Real> &var, std::vector<Real> &obj, std::vector<Real> &cons) const {
		std::vector<Real> var_(var);
		m_function->evaluate(var_, obj);
		for (auto &it : m_trans_objs)
			it->transfer(obj,var);
		int cons_id(0);
		if (cons.size() != m_pro->numberConstraints())
			cons.resize(m_pro->numberConstraints());
		for (auto &it : m_constraints)
			cons[cons_id++] = it->evaluate(var_);
	}

	void Subproblem::transferX(std::vector<Real>& vecX, const std::vector<Real>& var) {
		for (auto& it : m_trans_vars)
			it->transfer(vecX, var);
	}

	void Subproblem::readFromFile(const std::string &file_path) {
		std::string subspace_name;
		if (static_cast<ParameterType>(m_param.at("subspace").index()) == ParameterType::kChar)
			subspace_name.push_back(m_param.get<char>("subspace"));
		else if (static_cast<ParameterType>(m_param.at("subspace").index()) == ParameterType::kString)
			subspace_name = m_param.get<std::string>("subspace");

		std::ifstream in(file_path);

		ParameterMap param_dis;
		ParamIO::readParameter(in, param_dis);
		m_distance.reset(FactoryFP<DistanceBase>::produce(param_dis.get<std::string>("distance"),
			m_pro, subspace_name, param_dis));

		ParameterMap param_fun;
		ParamIO::readParameter(in, param_fun);
		m_function.reset(FactoryFP<FunctionBase>::produce(param_fun.get<std::string>("function"),
			m_pro, subspace_name, param_fun));

		std::string name;
		int size(0);
		in >> name >> size;
		while (size--) {
			ParameterMap param_con;
			ParamIO::readParameter(in, param_con);
			m_constraints.emplace_back(FactoryFP<ConstraintBase>::produce(param_con.get<std::string>("constraint"), 
				m_pro, subspace_name, param_con));
		}

		in >> name >> size;
		while (size--) {
			ParameterMap param_trv;
			ParamIO::readParameter(in, param_trv);
			m_trans_vars.emplace_back(FactoryFP<TransformBase>::produce(param_trv.get<std::string>("transform"), 
				m_pro, subspace_name, param_trv));
		}

		in >> name >> size;
		while (size--) {
			ParameterMap param_tro;
			ParamIO::readParameter(in, param_tro);
			m_trans_objs.emplace_back(FactoryFP<TransformBase>::produce(param_tro.get<std::string>("transform"), 
				m_pro, subspace_name, param_tro));
		}
	}

	void Subproblem::writeToFile(const std::string &file_path) const {
		std::ofstream out(file_path);
		ParamIO::writeParameter(m_distance->parameter(), out);
		out << "\n\n";
		ParamIO::writeParameter(m_function->parameter(), out);
		out << "\n\n";
		out << "num_constraints\t" << m_constraints.size() << '\n';
		for (auto &it : m_constraints) {
			out << '\n';
			ParamIO::writeParameter(it->parameter(), out);
		}
		out << "\n\n";
		out << "num_trans_vars\t" << m_trans_vars.size() << '\n';
		for (auto &it : m_trans_vars) {
			out << '\n';
			ParamIO::writeParameter(it->parameter(), out);
		}
		out << "\n\n";
		out << "num_trans_objs\t" << m_trans_objs.size() << '\n';
		for (auto &it : m_trans_objs) {
			out << '\n';
			ParamIO::writeParameter(it->parameter(), out);
		}
		out << std::endl;
	}
}
