#include "free_peaks_init_pop_bound.h"

namespace ofec {
	void FreePeaksInitPopBounded::addInputParameters() {
		m_input_parameters.add("index of peak to initialize", new RangedSizeT(m_id_peak_to_init, 0, 10, 0));
	}

	void FreePeaksInitPopBounded::initialize_(Environment *env) {
		FreePeaks::initialize_(env);
		if (m_id_peak_to_init < m_subspace_tree.tree->size()) {
			setDomainInitPop(m_subspace_tree.tree->getBox(m_id_peak_to_init));
		}
		else {
			throw Exception("The index of peak is out of range.");
		}
	}
}
