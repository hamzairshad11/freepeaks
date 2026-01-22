#include "map_x_ill_conditioning.h"


#include "../../free_peaks.h"


namespace ofec::free_peaks {
	MapXIllConditioning::MapXIllConditioning(Problem *pro, const std::string& subspace_name, const ParameterMap& param) :
		TransformBase(pro, subspace_name, param) {
		
		if (param.has("condition")) {
			m_condition = param.get<double>("condition");
		}


		int numVar = pro->numberVariables();
		mp_rotationMatrix.resize(numVar, numVar);
		;
		auto& uniform = pro->random()->uniform;
		auto& normal = pro->random()->normal;

		mp_rotationMatrix.randomize(&uniform);
		mp_rotationMatrix.generateRotationClassical(&normal, m_condition);
	}

	void MapXIllConditioning::transfer(std::vector<Real>& x_, const std::vector<Real>& var) {
		std::vector<double> x(x_);
		for (int i = 0; i < x.size(); i++) {
			x[i] = 0;

			for (int j = 0; j < x.size(); j++) {
				x[i] += mp_rotationMatrix[j][i] * x_[j];
			}
		}
		swap(x, x_);
	}

	void MapXIllConditioning::bindData() {

	}
}