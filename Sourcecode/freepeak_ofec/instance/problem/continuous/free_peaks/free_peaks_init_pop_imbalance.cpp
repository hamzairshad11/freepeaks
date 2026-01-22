#include "free_peaks_init_pop_imbalance.h"

namespace ofec {
	void FreePeaksInitPopImbalance::addInputParameters() {
		m_input_parameters.add("ratio of probabilities", new InputString(m_prob_ratio_init, "1"));
	}

	void FreePeaksInitPopImbalance::initialize_(Environment *env) {
		FreePeaks::initialize_(env);
		if (m_prob_ratio_init != "1") {
			m_each_prob.clear();
			std::stringstream ss(m_prob_ratio_init);
			std::string temp;
			while (getline(ss, temp, ',')) {
				m_each_prob.push_back(std::stod(temp));
			}
			if (m_each_prob.size() != m_subspace_tree.tree->size()) {
				throw Exception("Wrong ratio of probabilities.");
			}
			m_sum_prob = 0;
			for (auto prob : m_each_prob) {
				m_sum_prob += prob;
			}
		}
	}

	void FreePeaksInitPopImbalance::initializeVariables(VariableBase &vars, Random *rnd) const {
		if (m_prob_ratio_init == "1") {
			FreePeaks::initializeVariables(vars, rnd);
		}
		else {
			Real rand_pos = m_sum_prob * rnd->uniform.next();
			Real accum = 0;
			size_t id_peak_to_init = 0;
			for (auto prob : m_each_prob) {
				accum += prob;
				if (rand_pos <= accum) {
					break;
				}
				id_peak_to_init++;
			}
			auto &bnd = m_subspace_tree.tree->getBox(id_peak_to_init);
			auto &x = dynamic_cast<VariableVector<Real>&>(vars);
			for (int i = 0; i < m_number_variables; ++i) {
				x[i] = rnd->uniform.nextNonStd(bnd[i].first, bnd[i].second);
			}
		}
	}
}
