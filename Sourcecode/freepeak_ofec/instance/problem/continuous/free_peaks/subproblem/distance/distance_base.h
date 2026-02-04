#ifndef OFEC_FREE_PEAKS_DISTANCE_BASE_H
#define OFEC_FREE_PEAKS_DISTANCE_BASE_H

#include "../../../../../../core/parameter/parameter_map.h"
#include <vector>
#include "../../../../../../core/environment/environment.h"

namespace ofec {
	class Problem;
	namespace free_peaks {
		
		class  DistanceBase : virtual public Instance {
			OFEC_ABSTRACT_INSTANCE(DistanceBase)

		private:
			std::string m_register_name = "DistanceBase";

		protected:

			// to do : Change members, m_subspace_name and m_dis to private, and provide accessor functions that return const pointers.

			Problem* m_pro = nullptr;
			std::string m_subspace_name;
			static const Real ms_infinity;
		

		public:
			void setRegisterName(const std::string& name) {
				m_register_name = name;
			}
			void addInputParameters();
			virtual void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param);
			virtual ~DistanceBase() = default;
			virtual Real operator()(const std::vector<Real>& a, const std::vector<Real>& b) const = 0;
			virtual void bindData() {}
			//const ParameterMap& parameter() const { return m_param; }
			Real infiniy() const { return ms_infinity; }
		};
	}
}

#endif // !OFEC_FREE_PEAKS_DISTANCE_BASE_H