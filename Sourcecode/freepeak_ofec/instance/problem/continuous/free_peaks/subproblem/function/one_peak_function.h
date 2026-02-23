#ifndef OFEC_FREE_PEAKS_ONE_PEAK_FUNCTION_H
#define OFEC_FREE_PEAKS_ONE_PEAK_FUNCTION_H

#include "function_base.h"
#include "one_peak/one_peak_base.h"
#include <memory>

namespace ofec::free_peaks {

#define CAST_OnePeakFunction(fun) dynamic_cast<OnePeakFunction*>(fun)


	class OnePeakFunction : public FunctionBase {

		OFEC_CONCRETE_INSTANCE(OnePeakFunction)

	protected:
		virtual void evaluate_(const std::vector<double>& var, std::vector<double>& obj) override;
		void transferX(std::vector<Real>& vecX, const std::vector<Real>& var, const std::vector<ofec::Real>& loc)const;


	public:

		static const std::string ms_file_dir;

		static std::string directory();

		void addInputParameters();
		virtual void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) override;

		size_t numObjs() const override { return m_onepeaks.size(); }
		Real computeMaxDis(const std::vector<Real>& center);
		Real computeMinDis(const std::vector<Real>& center);

		const std::vector<std::unique_ptr<OnePeakBase>>& onepeaksFun() {
			return m_onepeaks;
		}
		virtual void writeToFile() const override;

		// construt functions
		void clearOnePeaks() {
			m_onepeaks.clear();
		}
		void addOnePeaks(const ParameterMap& param);
		void addOnePeaks(OnePeakBase* one_peak) {
			m_onepeaks.emplace_back(one_peak);
		}


		virtual void bindData();
		virtual void recordInputParameters()override;
		virtual void restoreInputParameters()override;
	private:


		std::vector<std::unique_ptr<OnePeakBase>> m_onepeaks;
		const DistanceBase* m_dis;
		void readFromFile(const std::string& file_path);
		std::vector<std::vector<ofec::Real>> m_transfered_centers;

		std::string m_generation_type;
		std::string m_file_name;

	};
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_FUNCTION_H