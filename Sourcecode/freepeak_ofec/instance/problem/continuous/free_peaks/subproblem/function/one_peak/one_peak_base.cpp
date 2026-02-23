#include "one_peak_base.h"
#include "../../../type.h"
#include "../../../free_peaks.h"
#include <fstream>

namespace ofec::free_peaks {


	const std::pair<Real, Real> OnePeakBase::ms_x_range = { -100, 100 };
	//	const std::pair<Real, Real> OnePeakBase::ms_y_range = {0, 100};
	const std::string OnePeakBase::ms_register_key = "one_peak";
	void OnePeakBase::addInputParameters() {
		m_input_parameters.add(ms_register_key, new InputParameterValueTypeString(m_register_name));

		std::list<std::string> group = { "default", "random", "assigned", "designated" };
		m_input_parameters.add("center_type", new EnumeratedString(m_center_type, group, "default"));
		// center position in the range [-100, 100]^D
		m_input_parameters.add("center_postion", new InputParameterValueTypeVectorReal(m_center));
		m_input_parameters.add("height", new RangedReal(m_height, 0, 100, 100));

		std::list<int> map_flag = { 0,1 };
		m_input_parameters.add("map_flag", new EnumeratedInt(m_map_flag, map_flag, 0));
		m_input_parameters.add("alpha", new RangedReal(m_alpha, 1e-6, 1e6, 1.0));
	}


	void OnePeakBase::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param)

	{

		m_pro = (pro);
		m_subspace_name = (subspace_name);
		m_input_parameters.input(param);
		if (m_pro->numberObjectives() == 1) {
			if (m_center_type == "default")
				m_center.assign(CAST_CONOP(m_pro)->numberVariables(), 0);
			else if (m_center_type == "random") {
				m_center.resize(CAST_CONOP(m_pro)->numberVariables());
				for (size_t i = 0; i < m_center.size(); ++i) {
					m_center[i] = m_pro->random()->uniform.nextNonStd<Real>(ms_x_range.first, ms_x_range.second);// todo
				}
				m_center_type = "assigned";
			}
			else if (m_center_type == "assigned") {
				if (m_center.size() != CAST_CONOP(m_pro)->numberVariables())
					throw Exception("OnePeakBase::initialize() center size incorrect");
			}
		}
		else if (m_pro->numberObjectives() > 1) {
			if (m_center_type == "default")
				m_center.assign(CAST_CONOP(m_pro)->numberVariables(), 0);
			else if (m_center_type == "random") {
				m_center.resize(CAST_CONOP(m_pro)->numberVariables());
				for (size_t i = 0; i < m_center.size(); ++i) {
					m_center[i] = m_pro->random()->uniform.nextNonStd<Real>(ms_x_range.first, ms_x_range.second);// todo
				}
			}
			else if (m_center_type == "designated") {
				m_center.resize(CAST_CONOP(m_pro)->numberVariables());
				for (size_t i = 0; i < m_center.size(); ++i) {
					m_center[i] = 10;// todo
				}
			}
		}
	}


}