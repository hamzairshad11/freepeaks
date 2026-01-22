#ifndef OFEC_FREE_PEAKS_FUNCTION_BASE_H
#define OFEC_FREE_PEAKS_FUNCTION_BASE_H

//#include "../../../../../../utility/parameter/param_map.h"
#include "../../../../../../core/problem/encoding.h"
#include "../../../../../../core/problem/problem.h"

namespace ofec {
	namespace free_peaks {
		class FunctionBase {
		protected:
			Problem *m_pro = nullptr;
			const std::string m_subspace_name;
			const ParameterMap m_param;
			std::vector<std::pair<Real, Real>> m_var_ranges;
			Real m_domain_size;
			std::vector<std::pair<Real, Real>> m_obj_ranges;
			std::vector<VariableVector<Real>> m_optimal_vars;
			std::vector<VariableVector<Real>> m_optimal_vars_inFunction;
			bool m_optimal_vars_mapped;
			virtual void evaluate_(const std::vector<double> &var, std::vector<double> &obj) = 0;
			virtual void transferX(std::vector<Real> &vecX, const std::vector<Real> &var)const;

		public:
			FunctionBase(Problem *pro, const std::string &subspace_name, const ParameterMap &param);
			virtual ~FunctionBase() = default;
			void evaluate(const std::vector<Real> &var, std::vector<double> &obj);
			virtual void bindData() = 0;
			const ParameterMap& parameter() const { return m_param; }
			const std::vector<std::pair<Real, Real>> &varRanges() const { return m_var_ranges; }
			Real domainSize() const { return m_domain_size; }
			const std::vector<std::pair<Real, Real>> &objRanges() const { return m_obj_ranges; }
			virtual size_t numObjs() const = 0;
			void mapOptimalVars();
			const std::vector<VariableVector<Real>> &optimalVars() const { return m_optimal_vars; };
			const std::vector<VariableVector<Real>> &optimalVarsFunction() const { return m_optimal_vars_inFunction; };

		};
	}
}
#endif // !OFEC_FREE_PEAKS_FUNCTION_BASE_H