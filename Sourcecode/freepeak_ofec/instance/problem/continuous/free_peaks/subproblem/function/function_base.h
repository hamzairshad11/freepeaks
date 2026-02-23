#ifndef OFEC_FREE_PEAKS_FUNCTION_BASE_H
#define OFEC_FREE_PEAKS_FUNCTION_BASE_H

//#include "../../../../../../utility/parameter/param_map.h"
#include "../../../../../../core/problem/encoding.h"
#include "../../../../../../core/problem/problem.h"

namespace ofec {
	namespace free_peaks {
		class FunctionBase : virtual public Instance {
			OFEC_ABSTRACT_INSTANCE(FunctionBase)

		private:
			std::string m_register_name;
		protected:

			// to do : Change members, m_subspace_name and m_dis to private, and provide accessor functions that return const pointers.

			Problem* m_pro = nullptr;
			std::string m_subspace_name;
			//const ParameterMap m_param;


			std::vector<std::pair<Real, Real>> m_var_ranges;
			Real m_domain_size;
			std::vector<std::pair<Real, Real>> m_obj_ranges;
			std::vector<VariableVector<Real>> m_optimal_vars;
			std::vector<VariableVector<Real>> m_optimal_vars_inFunction;
			bool m_optimal_vars_mapped;
			virtual void evaluate_(const std::vector<double>& var, std::vector<double>& obj) = 0;
			virtual void transferX(std::vector<Real>& vecX, const std::vector<Real>& var)const;

		public:

			void setRegisterName(const std::string& name) {
				m_register_name = name;
			}
			void addInputParameters();
			virtual void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param);
			virtual ~FunctionBase() = default;
			virtual void bindData() {};


			void evaluate(const std::vector<Real>& var, std::vector<double>& obj);
			//const ParameterMap& parameter() const { return m_param; }
			const std::vector<std::pair<Real, Real>>& varRanges() const { return m_var_ranges; }
			Real domainSize() const { return m_domain_size; }
			const std::vector<std::pair<Real, Real>>& objRanges() const { return m_obj_ranges; }
			virtual size_t numObjs() const = 0;
			void mapOptimalVars();
			const std::vector<VariableVector<Real>>& optimalVars() const { return m_optimal_vars; };
			const std::vector<VariableVector<Real>>& optimalVarsFunction() const { return m_optimal_vars_inFunction; };

			virtual void writeToFile() const {}


		};
	}
}
#endif // !OFEC_FREE_PEAKS_FUNCTION_BASE_H