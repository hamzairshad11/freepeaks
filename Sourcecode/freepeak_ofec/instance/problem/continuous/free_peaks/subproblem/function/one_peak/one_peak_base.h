#ifndef OFEC_FREE_PEAKS_ONE_PEAK_BASE_H
#define OFEC_FREE_PEAKS_ONE_PEAK_BASE_H

#include "../../distance/distance_base.h"
//#include "../../../../../../../utility/functional.h"

namespace ofec {
	class Problem;
	namespace free_peaks {
		class OnePeakBase : virtual public Instance {
			OFEC_ABSTRACT_INSTANCE(OnePeakBase)

		private:
			std::string m_register_name = "OnePeakBase";

		public:
			static const std::pair<Real, Real> ms_x_range;
			static const std::string ms_register_key;

			//	static const std::pair<Real, Real> ms_y_range;	

		protected:

			Problem* m_pro = nullptr;
			std::string m_subspace_name;
			//const ParameterMap m_param			
			std::string m_center_type;


			double m_height;
			std::vector<Real> m_center;
			double m_alpha = 1.0;
			const double m_maxValue = 100;
			int m_map_flag = 0;


			double mapperFun(Real dummy) {
				return std::pow(dummy / (m_maxValue + 1e-6), m_alpha) * (m_maxValue + 1e-6);
			}


			virtual Real evaluate_(Real dummy, size_t var_size) = 0;

		public:

			void setRegisterName(const std::string& name) {
				m_register_name = name;
			}
			void addInputParameters();
			virtual void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param);
			virtual ~OnePeakBase() = default;

			virtual void bindData() {}

			void setCenter(size_t i, Real v) { m_center[i] = v; }
			virtual Real evaluate(Real dummy, size_t var_size) {
				if (m_map_flag == 0) {
					return evaluate_(dummy, var_size);
				}
				else {
					return evaluate_(mapperFun(dummy), var_size);
				}
			}
			double height() const { return m_height; }
			const std::vector<Real>& center() const { return m_center; }
			//static const std::pair<Real, Real>& range()const { return ms_range; }

		};
	}
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_BASE_H