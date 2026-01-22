#ifndef OFEC_FREE_PEAKS_ONE_PEAK_FUNCTION_H
#define OFEC_FREE_PEAKS_ONE_PEAK_FUNCTION_H

#include "function_base.h"
#include "one_peak/one_peak_base.h"
#include <memory>

namespace ofec::free_peaks {
	class OnePeakFunction : public FunctionBase {
	protected:
		virtual void evaluate_(const std::vector<double> &var, std::vector<double> &obj) override;
		void transferX(std::vector<Real>& vecX, const std::vector<Real>& var, const std::vector<ofec::Real>& loc)const ;
	
	public:
		OnePeakFunction(Problem *pro, const std::string &subspace_name, const ParameterMap &param);
		void bindData() override;
		size_t numObjs() const override { return m_onepeaks.size(); }
		Real computeMaxDis(const std::vector<Real> &center);
		Real computeMinDis(const std::vector<Real> &center);

		const std::vector<std::unique_ptr<OnePeakBase>>& onepeaksFun() {
			return m_onepeaks;
		}

	private:
		std::vector<std::unique_ptr<OnePeakBase>> m_onepeaks;
		const DistanceBase *m_dis;
		void readFromFile(const std::string &file_path);
	
		std::vector<std::vector<ofec::Real>> m_transfered_centers;
	};
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_FUNCTION_H
