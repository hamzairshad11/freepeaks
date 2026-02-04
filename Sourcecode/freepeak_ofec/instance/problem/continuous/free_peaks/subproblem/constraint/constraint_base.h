#ifndef OFEC_FREE_PEAKS_CONSTRAINT_BASE_H
#define OFEC_FREE_PEAKS_CONSTRAINT_BASE_H

#include "../../../../../../core/parameter/parameter_map.h"
#include "../../../../../../core/instance.h"
#include "../distance/distance_base.h"
#include "../../type.h"


namespace ofec {
	namespace free_peaks {
		class ConstraintBase : virtual public Instance {
			OFEC_ABSTRACT_INSTANCE(ConstraintBase)

		private:
			std::string m_register_name = "ConstraintBase";
		protected:



			// to do : Change members, m_subspace_name and m_dis to private, and provide accessor functions that return const pointers.
			Problem* m_pro = nullptr;
			std::string m_subspace_name;
			DistanceBase* m_dis;

			
		public:


			void setRegisterName(const std::string& name) { 
				m_register_name = name; 
			}	
			void addInputParameters() {}
			void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param);
			virtual ~ConstraintBase() = default;
			virtual double evaluate(const std::vector<Real>& var) const = 0;
			virtual void bindData();
		};
	}
}

#endif // !OFEC_FREE_PEAKS_CONSTRAINT_BASE_H
