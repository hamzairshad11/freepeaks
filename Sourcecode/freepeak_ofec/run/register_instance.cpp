#include "interface.h"

#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../instance/problem/continuous/free_peaks/free_peaks_init_pop_bound.h"
#include "../instance/problem/continuous/free_peaks/free_peaks_init_pop_imbalance.h"
#include "../instance/problem/continuous/free_peaks/free_peaks_multiparty.h"


namespace ofec {
	void registerInstance() {

		REGISTER(Problem, FreePeaks, "free_peaks", std::set<std::string>({ "continuous" }));
		REGISTER(Problem, FreePeaksInitPopBounded, "free_peaks_init_pop_bounded", std::set<std::string>({ "continuous" }));
		REGISTER(Problem, FreePeaksInitPopImbalance, "free_peaks_init_pop_imbalance", std::set<std::string>({ "continuous" }));
		REGISTER(Problem, FreePeaksMultiParty, "free_peaks_multiparty", std::set<std::string>({ "continuous" }));


		for (auto &pro : g_factory_problem.get()) {
			for (auto &alg : g_factory_algorithm.get()) {
				for (auto it1 = pro.second.second.begin();;) {
					if (!alg.second.second.count(*it1)) {
						break;
					}
					if (++it1 == pro.second.second.end()) {
						if (g_factory_environment.get().count(pro.first)) {
							for (auto it2 = g_factory_environment.get().at(pro.first).second.begin();;) {
								if (!alg.second.second.count(*it2)) {
									break;
								}
								if (++it2 == g_factory_environment.get().at(pro.first).second.end()) {
									g_algorithm_for_environment_problem[pro.first].insert(alg.first);
								}
							}
						}
						else {
							g_algorithm_for_environment_problem[pro.first].insert(alg.first);
						}
						break;
					}
				}
			}
		}
	}
}
