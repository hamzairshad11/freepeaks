#include "one_peak_function.h"
#include "../../free_peaks.h"
#include "../../param_creator.h"
#include "../../factory.h"
#include "../../../../../../core/global.h"
#include <fstream>

namespace ofec::free_peaks {
	OnePeakFunction::OnePeakFunction(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		FunctionBase(pro, subspace_name, param), m_dis(nullptr)
	{
		auto generation = m_param.get<std::string>("generation");
		if (generation == "read_file") {
			auto file_name = m_param.get<std::string>("file") + ".txt";
			readFromFile(g_working_directory + "/instance/problem/continuous/free_peaks/subproblem/function/one_peak/" + file_name);
		}
		else if (generation == "random") {

		}
		else {

		}

	}

	void OnePeakFunction::bindData() {
		//FunctionBase::bindData();
		
		m_var_ranges.assign(CAST_CONOP(m_pro)->numberVariables(), { -100, 100 });

		m_domain_size = 0;
		for (auto &r : m_var_ranges)
			m_domain_size += pow(r.second - r.first, 2);
		m_domain_size = sqrt(m_domain_size);
		m_dis = CAST_FPs(m_pro)->subproblem(m_subspace_name)->distance();
		for (auto it = m_onepeaks.begin(); it != m_onepeaks.end(); ++it) {
			(*it)->bindData();
		}

		//Real theta = 2 * OFEC_PI / m_onepeaks.size();
		//for (auto it = m_onepeaks.begin(); it != m_onepeaks.end();++it) {
		//	(*it)->bindData();
		//	auto inx = std::distance(m_onepeaks.begin(), it);
		//	Real rotate_angle = theta * inx;
		//	auto first_center = (*it)->center();
		//	//rotate the first center to obtain the next center of objs
		//	Real v1 = (*it)->center()[0] * std::sqrt(2) * std::sin(OFEC_PI/4+rotate_angle);
		//	Real v2 = (*it)->center()[0] * std::sqrt(2) * std::cos(OFEC_PI / 4 + rotate_angle);
		//	(*it)->setCenter(0, v1);
		//	(*it)->setCenter(1, v2);
		//}
			



		{
			m_transfered_centers.resize(m_onepeaks.size());
			for (int idx(0); idx < m_onepeaks.size(); ++idx) {
				m_transfered_centers[idx] = m_onepeaks[idx]->center();
				transferX(m_transfered_centers[idx], m_transfered_centers[idx], m_onepeaks[idx]->center());
				///m_transfered_centers.push_back(transfer)
			}
		}

		m_obj_ranges.resize(m_onepeaks.size());
		for (size_t j = 0; j < m_obj_ranges.size(); ++j) {
			auto max_dis = computeMaxDis(m_onepeaks[j]->center());
			m_obj_ranges[j].first = m_onepeaks[j]->evaluate(max_dis, m_var_ranges.size());
			m_obj_ranges[j].second = m_onepeaks[j]->height();
		}
		if (m_onepeaks.size() == 1) {
			m_optimal_vars.push_back(m_onepeaks.front()->center());
			m_optimal_vars_mapped = false;
			m_optimal_vars_inFunction = m_optimal_vars;
		}
		else if (m_onepeaks.size() > 1) {
			for (auto& it:m_onepeaks) {
				m_optimal_vars.push_back(it->center());
			}
			m_optimal_vars_mapped = false;
			m_optimal_vars_inFunction = m_optimal_vars;
		}

		
	}

	void OnePeakFunction::readFromFile(const std::string &file_path) {
		std::ifstream in(file_path);
		std::string name;
		int size(0);
		in >> name >> size;
		while (size--) {
			ParameterMap param;
			ParamIO::readParameter(in, param);
			m_onepeaks.emplace_back(FactoryFP<OnePeakBase>::produce(param.get<std::string>("one_peak"),
				m_pro, m_subspace_name, param));
		}
	}

	void OnePeakFunction::evaluate_(const std::vector<double> &var, std::vector<double> &obj) {
		//
		//for (size_t j = 0; j < obj.size(); ++j) {
		//	obj[j] = m_onepeaks[j]->evaluate((*m_dis)(var, m_onepeaks[j]->center()), var.size());
		//}
		
		
		for (size_t j = 0; j < m_pro->numberObjectives(); ++j) {
			auto x_ = var;
			transferX(x_, var, m_onepeaks[j]->center());
			obj[j] = m_onepeaks[j]->evaluate((*m_dis)(x_, m_transfered_centers[j]), var.size());
		}
	}

	void OnePeakFunction::transferX(std::vector<Real>& vecX, const std::vector<Real>& var, const std::vector<ofec::Real>& loc)const {
		for (int idx(0); idx < vecX.size(); ++idx) {
			vecX[idx] -= loc[idx];
		}
		FunctionBase::transferX(vecX, var);
	}



	Real OnePeakFunction::computeMaxDis(const std::vector<Real> &center) {

		std::vector<Real> from(m_var_ranges.size()), to(m_var_ranges.size());
		for (size_t i = 0; i < from.size(); ++i) {
			
			double l = m_var_ranges[i].first, u = m_var_ranges[i].second;
			from[i] = l;
			to[i] = u;
		}

		transferX(from, from, center);
		transferX(to, to, center);
		auto center_x = center;


		std::vector<Real> fpoint(m_var_ranges.size());
		for (size_t i = 0; i < fpoint.size(); ++i) {
			double l = from[i], u = to[i];
			if (u - center[i] > center[i] - l)
				fpoint[i] = u;
			else
				fpoint[i] = l;
		}


		transferX(center_x, center_x, center);
		return (*m_dis)(fpoint, center_x);


		//std::vector<Real> fpoint(m_var_ranges.size());
		//for (size_t i = 0; i < fpoint.size(); ++i) {
		//	double l = m_var_ranges[i].first, u = m_var_ranges[i].second;
		//	if (u - center[i] > center[i] - l)
		//		fpoint[i] = u;
		//	else
		//		fpoint[i] = l;
		//}

		//transferX(fpoint, fpoint, center);
		//auto center_x = center;

		//transferX(center_x, center_x, center);
		//return (*m_dis)(fpoint, center_x);
	}
	
	Real OnePeakFunction::computeMinDis(const std::vector<Real> &center) {
		Real min_dis = m_dis->infiniy(); 
		size_t dim_min;
		bool is_lower_bnd;
		for (size_t i = 0; i < center.size(); ++i) {
			double l = m_var_ranges[i].first, u = m_var_ranges[i].second;
			if (u - center[i] > center[i] - l && min_dis > center[i] - l) {
				min_dis = center[i] - l; dim_min = i; is_lower_bnd = true;
			}
			else if (min_dis > u - center[i]) {
				min_dis = u - center[i]; dim_min = i; is_lower_bnd = false;
			}
		}
		std::vector<Real> cpoint(center);
		cpoint[dim_min] = is_lower_bnd ? m_var_ranges[dim_min].first : m_var_ranges[dim_min].second;

		transferX(cpoint, cpoint, center);

		auto center_x = center;
		transferX(center_x, center_x, center);

		return (*m_dis)(cpoint, center_x);
		return (*m_dis)(cpoint, center);
	}
}
