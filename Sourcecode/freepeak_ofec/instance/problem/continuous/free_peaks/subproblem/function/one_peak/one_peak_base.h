#ifndef OFEC_FREE_PEAKS_ONE_PEAK_BASE_H
#define OFEC_FREE_PEAKS_ONE_PEAK_BASE_H

#include "../../distance/distance_base.h"
//#include "../../../../../../../utility/functional.h"

namespace ofec {
	class Problem;
	namespace free_peaks {
		class OnePeakBase {
		protected:
			Problem* m_pro = nullptr;
			const std::string m_subspace_name;
			const ParameterMap m_param;
			double m_height;
			std::vector<Real> m_center;
			double m_alpha = 1.0;
			const double m_maxValue = 100;
			int m_map_flag = 0;
			//const std::pair<double, double> m_range = { -100, 100 };

			double mapperFun(Real dummy) {
				return std::pow(dummy / m_maxValue, m_alpha) * m_maxValue;
			}


			virtual Real evaluate_(Real dummy, size_t var_size) = 0;

		public:
			OnePeakBase(Problem *pro, const std::string& subspace_name, const ParameterMap& param);
			virtual ~OnePeakBase() = default;

			void setCenter(size_t i, Real v) { m_center[i] = v; }

			virtual Real evaluate(Real dummy, size_t var_size) {
				if (m_map_flag == 0) {
					return evaluate_(dummy, var_size);
				}
				else {
					return evaluate_(mapperFun(dummy), var_size);
				}
			}
			virtual void bindData();
			double height() const { return m_height; }
			const std::vector<Real>& center() const { return m_center; }
		};
	}
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_BASE_H
