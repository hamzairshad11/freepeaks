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

		std::list<Real> m_recommended_values = {
			1.0,
			10.0,
			50.0,
			100.0,
			500.0,
			1000.0,
			5000.0,
			10000.0,
			50000.0,
			100000.0
		};
	public:

		void addInputParameters();
		virtual void initialize(Problem *pro, const std::string& subspace_name, const ParameterMap& param)override;
		void transfer(std::vector<Real>& x, const std::vector<Real>& var) override;
	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
