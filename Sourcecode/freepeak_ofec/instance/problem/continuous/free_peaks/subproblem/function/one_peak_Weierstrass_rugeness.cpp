#include "map_Weierstrass_rugeness.h"
#include "../../../free_peaks.h"
#include "../../../../../../../utility/functional.h"

namespace ofec::free_peaks {


	void WeierstrassRugenessMap::addInputParameters() {


// 数学要求：与 b 共同满足 ab > 1 + 3π/2 ≈ 5.71 以保证无处可导
// 但数值实现中通常固定为 0.5 , a → 0：高频项迅速衰减，函数更平滑 a → 1：高频项衰减慢，函数更崎岖
		m_input_parameters.add("alpha", new RangedReal(m_alpha, 0, 1, 0.5));
        //取值范围：b > 1，通常为奇数, b 太大时，b^k 会迅速溢出（需要控制 kmax）
		m_input_parameters.add("b", new RangedInt(m_b, 1, 20, 3));
        //平衡精度与计算成本（CEC标准）
        m_input_parameters.add("kmax", new RangedInt(m_kmax, 1, 100, 20));  
	}

	void WeierstrassRugenessMap::initialize(Problem *pro, const std::string& subspace_name, const ParameterMap& param) {
		TransformBase::initialize(pro, subspace_name, param);
		//if (param.has("alpha")) {
		//	m_alpha = param.get<double>("alpha");
		//}
	}

    void WeierstrassRugenessMap::transfer(std::vector<Real>& obj, const std::vector<Real>& var) {
        auto m_number_variables = m_pro->numberVariables();
        auto m_number_objectives = m_pro->numberObjectives();


		auto& range = m_sub_pro->function()->varRanges();   

        // 标准 Weierstrass：变量定义域为 [-0.5, 0.5]
        std::vector<Real> x = var;
        for (int i = 0; i < m_number_variables; ++i) {
            x[i] /= m_sub_pro->function()->domainSize();
            x[i] = x[i] - 0.5;  // 将 [0,1] 映射到 [-0.5, 0.5]
        }

        double value = 0.0;

        // 标准 Weierstrass 计算公式
        for (int i = 0; i < m_number_variables; ++i) {
            double tmp = 0.0;
            // CEC2005 标准使用 20 项求和
            for (int k = 0; k < m_kmax; ++k) {  // k = 0 到 kmax-1
                tmp += m_aK[k] * cos(2 * OFEC_PI * m_bK[k] * (x[i] + 0.5));
            }
            value += tmp;
        }

        // 减去偏移项：D * sum_{k=0}^{kmax-1} a^k * cos(π * b^k)
        value -= m_number_variables * m_F0;

        // 赋值给目标函数
        for (int idO(0); idO < m_number_objectives; ++idO) {
            obj[idO] = value;  // 直接赋值，不要立方运算
        }


    }

    void WeierstrassRugenessMap::bindData() {
        X_TransformBase::bindData();
        // 标准 Weierstrass 参数设置


        m_aK.resize(m_kmax);
        m_bK.resize(m_kmax);
        m_F0 = 0.0;

        // 预计算系数和偏移常数
        for (size_t k = 0; k < m_kmax; ++k) {
            m_aK[k] = pow(m_alpha, (Real)k);          // a^k
            m_bK[k] = pow(m_b, (Real)k);              // b^k
            // F0 = sum_{k=0}^{kmax-1} a^k * cos(π * b^k)
            m_F0 += m_aK[k] * cos(OFEC_PI * m_bK[k]);
        }
    }
}