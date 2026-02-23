#include "map_x_ill_conditioning.h"


#include "../../../free_peaks.h"


namespace ofec::free_peaks {

	void MapXIllConditioning::addInputParameters() {
		m_input_parameters.add("condition", new EnumeratedReal(m_condition, m_recommended_values, 100.0));
	}

	void MapXIllConditioning::initialize(Problem *pro, const std::string& subspace_name, const ParameterMap& param) {
		TransformBase::initialize(pro, subspace_name, param);
		int numVar = pro->numberVariables();
		mp_rotationMatrix.resize(numVar, numVar);
		
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


}