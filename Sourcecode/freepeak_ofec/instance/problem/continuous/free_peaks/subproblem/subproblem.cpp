#include "subproblem.h"
#include <fstream>
#include "../factory.h"
#include "../../../../../core/global.h"
#include "../param_creator.h"

namespace ofec::free_peaks {

	const std::string Subproblem::ms_file_dir = "instance/problem/continuous/free_peaks/subproblem/";
	std::string Subproblem::directory() {
		return g_working_directory + ms_file_dir;
	}
	void Subproblem::addInputParameters() {
		std::list<std::string> group = { "read_file", "random", "assigned"};
		m_input_parameters.add("generation_type", new EnumeratedString(m_generation_type, group, "read_file"));
		m_input_parameters.add("dataFile1", new FileName(m_file_name,
			ms_file_dir,
			"config (*.txt)", "sop/s1.txt"));

		m_input_parameters.add("subspace", new InputParameterValueTypeString(m_subspace_name));

	}

	void Subproblem::initialize(const ParameterMap &param, Problem *pro)

	{
		m_pro = (pro);
		m_input_parameters.input(param);
		if (m_generation_type == "read_file") {
			readFromFile(g_working_directory + ms_file_dir + m_file_name);
			//setOptima();
		}
		else if (m_generation_type == "random") {
			// to do 

		}
		else if (m_generation_type == "assigned") {
			// do nothing
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

		for (auto& it : m_trans_vars) {
			it->bindData();
		}


		m_function->bindData();

		for (auto& it : m_constraints) {
			it->bindData();
		}
		for (auto& it : m_trans_objs) {
			it->bindData();
		}
		setOptima();
	}

	void Subproblem::evaluate(const std::vector<Real> &var, std::vector<Real> &obj, std::vector<Real> &cons) const {
		std::vector<Real> var_(var);
		m_function->evaluate(var_, obj);
		for (auto& it : m_trans_objs) {
			it->transfer(obj, var);
		}
		int cons_id(0);
		if (cons.size() != m_pro->numberConstraints())
			cons.resize(m_pro->numberConstraints());
		for (auto &it : m_constraints)
			cons[cons_id++] = it->evaluate(var_);
	}

	void Subproblem::transferX(std::vector<Real>& vecX, const std::vector<Real>& var) {
		for (auto& it : m_trans_vars) {
			it->transfer(vecX, var);
			//std::cout << "Subproblem::transferX" << std::endl;
			//for (int idx(0); idx < vecX.size(); ++idx)
			//	std::cout << " " << vecX[idx];
			//std::cout << std::endl;
		}

	}

	void Subproblem::readFromFile(const std::string &file_path) {
		//std::string subspace_name;
		//if (static_cast<ParameterType>(m_param.at("subspace").index()) == ParameterType::kChar)
		//	subspace_name.push_back(m_param.get<char>("subspace"));
		//else if (static_cast<ParameterType>(m_param.at("subspace").index()) == ParameterType::kString)
		//	subspace_name = m_param.get<std::string>("subspace");
		std::string name;
		std::ifstream in(file_path);

		{
			//in >> name;
			ParameterMap param_dis;
			in >> param_dis;
			//ParamIO::readParameter(in, param_dis);
			setDistance(param_dis);
		}


		{
			//in >> name;
			ParameterMap param_fun;
			in >> param_fun;
			//ParamIO::readParameter(in, param_fun);
			setFunction(param_fun);

		}


		int size(0);
		in >> name >> size;
		clearConstraints();
		while (size--) {
			ParameterMap param_con;
			in >> param_con;
			//ParamIO::readParameter(in, param_con);
			addConstraint(param_con);

		}

		in >> name >> size;
		clearVariableTransforms();
		while (size--) {
			ParameterMap param_trv;
			in >> param_trv;
			//ParamIO::readParameter(in, param_trv);
			addVariableTransform(param_trv);
		}

		in >> name >> size;
		clearObjectiveTransforms();
		while (size--) {
			ParameterMap param_tro;
			in >> param_tro;
		//	ParamIO::readParameter(in, param_tro);
			addObjectiveTransform(param_tro);
		}
	}

	void Subproblem::writeToFile() const {
		auto file_path = g_working_directory + ms_file_dir + m_file_name;
		
		std::ofstream out(file_path);
		//out << "distance\n";
		out << m_distance->archivedParameters();
		//ParamIO::writeParameter(m_distance->archivedParameters(), out);
		out << "\n\n";

		//out << "function\n";
		out << m_function->archivedParameters();
		//ParamIO::writeParameter(m_function->archivedParameters(), out);
		out << "\n\n";
		out << "num_constraints\t" << m_constraints.size() << '\n';
		for (auto &it : m_constraints) {
			out << '\n';
			out << it->archivedParameters();
			//ParamIO::writeParameter(it->archivedParameters(), out);
		}
		out << "\n\n";
		out << "num_trans_vars\t" << m_trans_vars.size() << '\n';
		for (auto &it : m_trans_vars) {
			out << '\n';
			out << it->archivedParameters();
			//ParamIO::writeParameter(it->archivedParameters(), out);
		}
		out << "\n\n";
		out << "num_trans_objs\t" << m_trans_objs.size() << '\n';
		for (auto &it : m_trans_objs) {
			out << '\n';
			out << it->archivedParameters();
			//ParamIO::writeParameter(it->archivedParameters(), out);
		}
		out << std::endl;

		out.close();
		
	}

	void Subproblem::setDistance(const ParameterMap& param_dis)
	{
		m_distance.reset(FactoryFP<DistanceBase>::produce(param_dis.get<std::string>("distance")));
		m_distance->initialize(m_pro, m_subspace_name, param_dis);
	//	m_distance->recordInputParameters();
	}

	void Subproblem::setFunction(const ParameterMap& param_fun)
	{
		m_function.reset(FactoryFP<FunctionBase>::produce(param_fun.get<std::string>("function")));
		m_function->initialize(m_pro, m_subspace_name, param_fun);
	//	m_function->recordInputParameters();
	}

	void Subproblem::addConstraint(const ParameterMap& param_con)
	{
		m_constraints.emplace_back(FactoryFP<ConstraintBase>::produce(param_con.get<std::string>("constraint")));
		m_constraints.back()->initialize(m_pro, m_subspace_name, param_con);
	//	m_constraints.back()->recordInputParameters();
	}

	void Subproblem::addVariableTransform(const ParameterMap& param_trv)
	{
		m_trans_vars.emplace_back(FactoryFP<X_TransformBase>::produce(param_trv.get<std::string>("transformX")));
		m_trans_vars.back()->initialize(m_pro, m_subspace_name, param_trv);
//		m_trans_vars.back()->recordInputParameters();
	}

	void Subproblem::addObjectiveTransform(const ParameterMap& param_tro)
	{
		m_trans_objs.emplace_back(FactoryFP<Y_TransformBase>::produce(param_tro.get<std::string>("transformY")));
		m_trans_objs.back()->initialize(m_pro, m_subspace_name, param_tro);
	//	m_trans_objs.back()->recordInputParameters();
	}

	void Subproblem::outputTotalFile() {
		writeToFile();
		m_function->writeToFile();	
	}

	void Subproblem::recordInputParameters()
	{
		Instance::recordInputParameters();
		m_function->recordInputParameters();	
		m_distance->recordInputParameters();
		for (auto& it : m_constraints) it->recordInputParameters();
		for (auto& it : m_trans_objs) it->recordInputParameters();
		for (auto& it : m_trans_vars) it->recordInputParameters();
	}

	void Subproblem::restoreInputParameters()
	{
		Instance::restoreInputParameters();
		m_function->restoreInputParameters();
		m_distance->restoreInputParameters();
		for (auto& it : m_constraints) it->restoreInputParameters();
		for (auto& it : m_trans_objs) it->restoreInputParameters();
		for (auto& it : m_trans_vars) it->restoreInputParameters();
	}

}
