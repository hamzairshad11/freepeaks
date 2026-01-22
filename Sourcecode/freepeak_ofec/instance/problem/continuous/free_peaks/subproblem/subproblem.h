#ifndef OFEC_FREE_PEAKS_SUBPROBLEM_H
#define OFEC_FREE_PEAKS_SUBPROBLEM_H

#include "function/function_base.h"
#include "constraint/constraint_base.h"
#include "distance/distance_base.h"
#include "transform/transform_base.h"
#include "../../../../../core/problem/optima.h"

namespace ofec {
	namespace free_peaks {
		class Subproblem {
		private:
			const ParameterMap m_param;
			Problem* m_pro = nullptr;
			std::unique_ptr<FunctionBase> m_function;
			std::unique_ptr<DistanceBase> m_distance;
			std::vector<std::unique_ptr<ConstraintBase>> m_constraints;
			std::vector<std::unique_ptr<TransformBase>> m_trans_vars;
			std::vector<std::unique_ptr<TransformBase>> m_trans_objs;
			Optima<> m_optima;

			void readFromFile(const std::string& file_path);
			void setOptima();

		public:
			Subproblem(const ParameterMap& param, Problem *pro);
			void evaluate(const std::vector<Real>& var, std::vector<Real>& obj, std::vector<Real>& cons) const;
			void bindData();
			const ParameterMap& parameter() const { return m_param; }
			virtual void transferX(std::vector<Real>& vecX, const std::vector<Real>& var);
			DistanceBase* distance() const { return m_distance.get(); }
			FunctionBase* function() const { return m_function.get(); }
			ConstraintBase* constraint(size_t i) const { return m_constraints[i].get(); }
			TransformBase* transVar(size_t i) const { return m_trans_vars[i].get(); }
			TransformBase* transObj(size_t i) const { return m_trans_objs[i].get(); }
			size_t numObjs() const { return m_function->numObjs(); }
			size_t numCons() const { return m_constraints.size(); }
			size_t numTransVars() const { return m_trans_vars.size(); }
			size_t numTransObjs() const { return m_trans_objs.size(); }
			const Optima<>& optima() const { return m_optima; }
			void writeToFile(const std::string& file_path) const;
		};
	}
}

#endif // !OFEC_FREE_PEAKS_SUBPROBLEM_H
