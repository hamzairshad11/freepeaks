#ifndef OFEC_CUSTOM_METHOD_HPP
#define OFEC_CUSTOM_METHOD_HPP

#include "../core/global.h"
#include "../instance/problem/continuous/free_peaks/free_peaks_init_pop_bound.h"
#include <fstream>

namespace ofec {
	void run() {
		ofec::g_working_directory = "D:/Document/GitLab/ofec-data";
		std::shared_ptr<Environment> env(Environment::create());
		env->recordInputParameters();
		env->initialize();

		auto free_peaks = FreePeaksInitPopBounded::create();
		free_peaks->inputParameters().at("dataFile1")->setValue("sop/2d_s9_s1_s10_s3.txt");
		free_peaks->recordInputParameters();

		env->setProblem(free_peaks);
		env->initializeProblem();

		size_t num_samples = 10000;
		auto &kdtree = free_peaks->subspaceTree().tree;
		std::unique_ptr<SolutionBase> sol(free_peaks->createSolution());
		Random rnd(0.5);
		std::vector<std::vector<Real>> data(num_samples, std::vector<Real>(kdtree->size()));
		for (size_t i = 0; i < kdtree->size(); i++) {
			free_peaks->inputParameters().at("index of peak to initialize")->setValue(std::to_string(i));
			free_peaks->recordInputParameters();
			env->initializeProblem();
			for (size_t j = 0; j < num_samples; j++) {
				sol->initialize(env.get(), &rnd);
				sol->evaluate(env.get(), false);
				data[j][i] = sol->objective(0);
			}
		}
		std::stringstream ss;
		for (size_t i = 0; i < kdtree->size() - 1; i++) {
			ss << "BoA" << std::to_string(i + 1) << ',';
		}
		ss << "BoA" << std::to_string(kdtree->size()) << '\n';
		for (auto &row : data) {
			for (size_t i = 0; i < row.size() - 1; ++i) {
				ss << row[i] << ',';
			}
			ss << row.back() << '\n';
		}
		auto data_file_name = free_peaks->fileName();
		for (auto &l : data_file_name) {
			if (l == '/' || l == '\\') {
				l = '_';
			}
		}
		size_t last_dot = data_file_name.find_last_of('.');
		data_file_name = data_file_name.substr(0, last_dot);
		std::ofstream ofs(ofec::g_working_directory + "/result/free_peaks__fit_dist__" + data_file_name + ".csv");
		ofs << ss.str();
	}
}

#endif // !OFEC_CUSTOM_METHOD_HPP