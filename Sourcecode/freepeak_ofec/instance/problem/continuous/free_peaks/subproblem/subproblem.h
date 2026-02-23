#ifndef OFEC_FREE_PEAKS_SUBPROBLEM_H
#define OFEC_FREE_PEAKS_SUBPROBLEM_H

#include "function/function_base.h"
#include "constraint/constraint_base.h"
#include "distance/distance_base.h"
#include "transform/transform_x/transform_x_base.h"
#include "transform/transform_y/transform_y_base.h"
#include "../../../../../core/problem/optima.h"
#include "../../../../../core/instance.h"


namespace ofec {
	namespace free_peaks {
		class Subproblem : virtual public Instance {

			OFEC_CONCRETE_INSTANCE(Subproblem)
		public:


			static const std::string ms_file_dir;
			static std::string directory();

		private:
			//const ParameterMap m_param;

			std::string m_generation_type = "read_file";
			std::string m_file_name;
			std::string m_subspace_name = "s1";

			Problem* m_pro = nullptr;


			std::unique_ptr<FunctionBase> m_function;
			std::unique_ptr<DistanceBase> m_distance;
			std::vector<std::unique_ptr<ConstraintBase>> m_constraints;
			std::vector<std::unique_ptr<X_TransformBase>> m_trans_vars;
			std::vector<std::unique_ptr<Y_TransformBase>> m_trans_objs;
			Optima<> m_optima;

			void readFromFile(const std::string& file_path);
			void setOptima();



		public:


			void addInputParameters();
			void initialize(const ParameterMap& param, Problem* pro);
			void evaluate(const std::vector<Real>& var, std::vector<Real>& obj, std::vector<Real>& cons) const;
			void bindData();
			//const ParameterMap& parameter() const { return m_param; }
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
			void writeToFile() const;


			// construction related functions;
			void setDistance(const ParameterMap& param_dis);
			void setDistance(DistanceBase* dis) {
				m_distance.reset(dis);
			}
			void setFunction(const ParameterMap& param_fun);
			void setFunction(FunctionBase* fun) {
				m_function.reset(fun);
			}
			void clearConstraints() {
				m_constraints.clear();
			}
			void addConstraint(const ParameterMap& param_con);
			void addConstraint(ConstraintBase* con) {
				m_constraints.emplace_back(con);
			}
			void clearVariableTransforms() {
				m_trans_vars.clear();
			}
			void addVariableTransform(const ParameterMap& param_trans);
			void addVariableTransform(X_TransformBase* tranX) {
				m_trans_vars.emplace_back(tranX);
			}
			void clearObjectiveTransforms() {
				m_trans_objs.clear();
			}
			void addObjectiveTransform(const ParameterMap& param_objs);
			void addObjectiveTransform(Y_TransformBase* tranY) {
				m_trans_objs.emplace_back(tranY);
			}

			void outputTotalFile();

			virtual void recordInputParameters();
			virtual void restoreInputParameters();
		};
	}
}

#endif // !OFEC_FREE_PEAKS_SUBPROBLEM_H