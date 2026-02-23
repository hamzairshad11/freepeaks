#ifndef OFEC_FREE_PEAKS_X_ILL_CONDITIONING_H
#define OFEC_FREE_PEAKS_X_ILL_CONDITIONING_H

#include "transform_x_base.h"
#include "../../../../../../../utility/linear_algebra/matrix.h"
namespace ofec::free_peaks {
	class MapXIllConditioning : public X_TransformBase {
		OFEC_CONCRETE_INSTANCE(MapXIllConditioning)
	private:
		Matrix mp_rotationMatrix;
		Real m_condition = 1.0;

		// 完整测试序列（覆盖所有难度）
		std::list<Real> m_recommended_values = {
			1.0,      // 基线
			10.0,     // 轻微
			50.0,     // 轻度
			100.0,    // 中等
			500.0,    // 中高
			1000.0,   // 高
			5000.0,   // 严重
			10000.0,  // 很严重
			50000.0,  // 极度
			100000.0  // 极限
		};
	public:

		void addInputParameters();
		virtual void initialize(Problem *pro, const std::string& subspace_name, const ParameterMap& param)override;
		void transfer(std::vector<Real>& x, const std::vector<Real>& var) override;
	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
