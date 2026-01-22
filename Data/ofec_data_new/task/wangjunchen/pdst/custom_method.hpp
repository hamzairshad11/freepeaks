#ifndef OFEC_STATISTICS_CUSTOM_METHOD_HPP
#define OFEC_STATISTICS_CUSTOM_METHOD_HPP

#include <run/interface.h>
#include <core/global.h>
#include "../statistics/global.h"
#include <list>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include "../record/record_pdst.h"
#include <filesystem>
#include <mutex>

namespace ofec_statistics {
	std::set<std::string> g_parameters_omit;
	std::map<std::string, std::string> g_parameter_abbreviation;
	std::mutex g_cout_mutex;

	void run();
	void runAlgorithms(std::list<std::shared_ptr<ofec::Environment>> &environments);
	void registerParameterAbbreviation();
	std::map<std::string, std::string> abbreviationToParameter(const std::string &abbreviations);
	void omitParameterInFileName();
	std::string parametersToFileName(const std::map<std::string, std::string> &parameters);

	void run() {
		ofec::registerInstance();
		registerParameterAbbreviation();
		omitParameterInFileName();
		std::string argvs;
		ofec::g_working_directory = "E:/Document/GitLab/ofec-data"; // path of ofec-data
		std::ifstream task_file(ofec::g_working_directory + "/task/wangjunchen/pdst/example.sh");
		while (std::getline(task_file, argvs)) {
			if (argvs.empty() || argvs.front() == '#') {
				continue;
			}
			auto parameters = abbreviationToParameter(argvs);
			if (!parameters.count("problem name") || !parameters.count("algorithm name")) {
				continue;
			}
			auto problem_name = parameters.at("problem name");
			auto algorithm_name = parameters.at("algorithm name");
			if (ofec::checkValidation(problem_name, algorithm_name)) {
				/* Start mark */
				std::cout << "[Task start] Parameters: " << argvs << std::endl;
				auto start = std::chrono::steady_clock::now();
				/* Set records */
				g_record_task.reset(new RecordTask);
				auto record = new RecordPDST;
				record->inputParameters().input(parameters);
				record->setFileName(parametersToFileName(parameters));
				g_record_task->add(record);
				/* Generate instances by factory */
				std::vector<ofec::Real> seeds;
				std::vector<std::shared_ptr<ofec::Environment>> environments;
				size_t num_runs = std::stoull(parameters.at("number of runs"));
				try {
					for (size_t i = 0; i < num_runs; i++) {
						ofec::Real seed = ofec::Real(i + 1) / (num_runs + 1);
						environments.emplace_back(ofec::generateEnvironmentByFactory(problem_name));
						environments.back()->setProblem(ofec::generateProblemByFactory(problem_name));
						environments.back()->setAlgorithm(ofec::generateAlgorithmByFactory(algorithm_name));
						seeds.push_back(seed);
					}
				}
				catch (const ofec::Exception &e) {
					std::cout << "Failed to generate instance by factory: " << e.what() << std::endl;
					break;
				}
				/* Initialize instances */
				try {
					for (size_t i = 0; i < num_runs; i++) {
						environments[i]->inputParameters().input(parameters);
						environments[i]->recordInputParameters();
						environments[i]->initialize(seeds[i]);
						environments[i]->problem()->inputParameters().input(parameters);
						environments[i]->problem()->recordInputParameters();
						environments[i]->initializeProblem(seeds[i]);
						environments[i]->algorithm()->inputParameters().input(parameters);
						environments[i]->algorithm()->recordInputParameters();
						environments[i]->initializeAlgorithm(seeds[i]);
					}
				}
				catch (const ofec::Exception &e) {
					std::cout << "Failed to initialize instance: " << e.what() << std::endl;
					break;
				}
				/* Run algorithms */
				std::vector<std::thread> threads;
				size_t num_threads = parameters.count("number of threads") ?
					std::stoull(parameters.at("number of threads")) :
					std::min<size_t>(num_runs, std::thread::hardware_concurrency());
				std::vector<std::list<std::shared_ptr<ofec::Environment>>> environments_each_thread(num_threads);
				size_t id_thread = 0;
				for (auto &environment : environments) {
					environments_each_thread[id_thread].push_back(environment);
					id_thread++;
					if (id_thread == num_threads) {
						id_thread = 0;
					}
				}
				for (auto &environments : environments_each_thread) {
					threads.emplace_back(runAlgorithms, std::ref(environments));
				}
				for (auto &t : threads) {
					t.join();
				}
				/* Output record */
				std::cout << "############ Output data ############" << std::endl;
				try {
					g_record_task->outputProgress();
					g_record_task->outputFinal();
				}
				catch (const ofec::Exception &e) {
					std::cout << "Failed to output record: " << e.what() << std::endl;
				}
				/* End mark */
				auto end = std::chrono::steady_clock::now();
				float diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
				std::cout << "[Task finished] Time used: " << diff / 1000 << " seconds.\n" << std::endl;
			}
		}
	}

	void registerParameterAbbreviation() {
		g_parameter_abbreviation.insert({ "AID", "algId" });
		g_parameter_abbreviation.insert({ "al", "alpha" });
		g_parameter_abbreviation.insert({ "AN", "algorithm name" });
		g_parameter_abbreviation.insert({ "be", "beta" });
		g_parameter_abbreviation.insert({ "BPF", "mean basin potential" });
		g_parameter_abbreviation.insert({ "C1", "accelerator1" });
		g_parameter_abbreviation.insert({ "C2", "accelerator2" });
		g_parameter_abbreviation.insert({ "ca", "case" });
		g_parameter_abbreviation.insert({ "CAN", "cluster add neighbors" });
		g_parameter_abbreviation.insert({ "CDBGFID", "comDBGFunID" });
		g_parameter_abbreviation.insert({ "CE", "crossover eta" });
		g_parameter_abbreviation.insert({ "CF", "changeFre" });
		g_parameter_abbreviation.insert({ "CoF", "convFactor" });
		g_parameter_abbreviation.insert({ "CPR", "changePeakRatio" });
		g_parameter_abbreviation.insert({ "CR", "crossover rate" });
		g_parameter_abbreviation.insert({ "CS", "changeSeverity" });
		g_parameter_abbreviation.insert({ "CT", "changeType" });
		g_parameter_abbreviation.insert({ "CT2", "changeType2" });
		g_parameter_abbreviation.insert({ "CTH", "convThreshold" });
		g_parameter_abbreviation.insert({ "DD", "diversity degree" });
		g_parameter_abbreviation.insert({ "DD1", "dataDirectory1" });
		g_parameter_abbreviation.insert({ "DF1", "dataFile1" });
		g_parameter_abbreviation.insert({ "DF2", "dataFile2" });
		g_parameter_abbreviation.insert({ "DF3", "dataFile3" });
		g_parameter_abbreviation.insert({ "DFR", "dataFileResult" });
		g_parameter_abbreviation.insert({ "DG", "durationGeneration" });
		g_parameter_abbreviation.insert({ "DM", "divisionMode" });
		g_parameter_abbreviation.insert({ "DR", "detectRatio" });
		g_parameter_abbreviation.insert({ "EBP", "evolve by potential" });
		g_parameter_abbreviation.insert({ "ECF", "evalCountFlag" });
		g_parameter_abbreviation.insert({ "ep", "epsilon" });
		g_parameter_abbreviation.insert({ "ER", "exlRadius" });
		g_parameter_abbreviation.insert({ "et", "eta" });
		g_parameter_abbreviation.insert({ "exitPS", "exploit pop size" });
		g_parameter_abbreviation.insert({ "exrePS", "explore pop size" });
		g_parameter_abbreviation.insert({ "F", "scaling factor" });
		g_parameter_abbreviation.insert({ "FA", "flagAsymmetry" });
		g_parameter_abbreviation.insert({ "FI", "flagIrregular" });
		g_parameter_abbreviation.insert({ "FN", "flagNoise" });
		g_parameter_abbreviation.insert({ "FNDC", "flagNumDimChange" });
		g_parameter_abbreviation.insert({ "FNPC", "flagNumPeakChange" });
		g_parameter_abbreviation.insert({ "FR", "flagRotation" });
		g_parameter_abbreviation.insert({ "FTL", "flagTimeLinkage" });
		g_parameter_abbreviation.insert({ "ga", "gamma" });
		g_parameter_abbreviation.insert({ "GM", "god mode" });
		g_parameter_abbreviation.insert({ "gOptFlag", "gOptFlag" });
		g_parameter_abbreviation.insert({ "GS", "glstrture" });
		g_parameter_abbreviation.insert({ "GT", "generation type" });
		g_parameter_abbreviation.insert({ "HCM", "heightConfigMode" });
		g_parameter_abbreviation.insert({ "HR", "hibernatingRadius" });
		g_parameter_abbreviation.insert({ "IC", "initial clustering" });
		g_parameter_abbreviation.insert({ "IT1", "interTest1" });
		g_parameter_abbreviation.insert({ "IT2", "interTest2" });
		g_parameter_abbreviation.insert({ "IT3", "interTest3" });
		g_parameter_abbreviation.insert({ "IT4", "interTest4" });
		g_parameter_abbreviation.insert({ "ITV", "interval number of solutions" });
		g_parameter_abbreviation.insert({ "JH", "jumpHeight" });
		g_parameter_abbreviation.insert({ "KF", "kFactor" });
		g_parameter_abbreviation.insert({ "LF", "lFactor" });
		g_parameter_abbreviation.insert({ "ME", "maximum evaluations" });
		g_parameter_abbreviation.insert({ "MF", "model file" });
		g_parameter_abbreviation.insert({ "MI", "maxIter" });
		g_parameter_abbreviation.insert({ "MNPS", "minNumPopSize" });
		g_parameter_abbreviation.insert({ "MP", "mutProbability" });
		g_parameter_abbreviation.insert({ "MR", "mutation rate" });
		g_parameter_abbreviation.insert({ "MRT", "maxRunTime" });
		g_parameter_abbreviation.insert({ "MSDE", "mutation strategy" });
		g_parameter_abbreviation.insert({ "MSI", "maxSucIter" });
		g_parameter_abbreviation.insert({ "MSS", "max subpop size" });
		g_parameter_abbreviation.insert({ "MuE", "mutation eta" });
		g_parameter_abbreviation.insert({ "NAS", "number of atomspaces" });
		g_parameter_abbreviation.insert({ "NB", "numBox" });
		g_parameter_abbreviation.insert({ "NC", "numChange" });
		g_parameter_abbreviation.insert({ "NCus", "numCus" });
		g_parameter_abbreviation.insert({ "ND", "number of variables" });
		g_parameter_abbreviation.insert({ "NDi", "number of divisions" });
		g_parameter_abbreviation.insert({ "NER", "number of solutions for exploration" });
		g_parameter_abbreviation.insert({ "NF", "noiseFlag" });
		g_parameter_abbreviation.insert({ "NGO", "numGOpt" });
		g_parameter_abbreviation.insert({ "NN", "number of neighbors" });
		g_parameter_abbreviation.insert({ "NNCus", "numNewCus" });
		g_parameter_abbreviation.insert({ "NO", "number of objectives" });
		g_parameter_abbreviation.insert({ "NoN", "number of neighbors" });
		g_parameter_abbreviation.insert({ "NP", "numPeak" });
		g_parameter_abbreviation.insert({ "NPa", "numPartition" });
		g_parameter_abbreviation.insert({ "NPR", "numParetoRegion" });
		g_parameter_abbreviation.insert({ "NPS", "numPS" });
		g_parameter_abbreviation.insert({ "NR", "number of runs" });
		g_parameter_abbreviation.insert({ "NRP", "number of reference points" });
		g_parameter_abbreviation.insert({ "NS", "noiseSeverity" });
		g_parameter_abbreviation.insert({ "NSP", "number of subpopulations" });
		g_parameter_abbreviation.insert({ "NSSP", "number of subspaces" });
		g_parameter_abbreviation.insert({ "NT", "number of threads" });
		g_parameter_abbreviation.insert({ "NumGPS", "number of global PS shapes" });
		g_parameter_abbreviation.insert({ "NumPri", "vector of private variables" });
		g_parameter_abbreviation.insert({ "NumPriPeak", "vector of private peaks" });
		g_parameter_abbreviation.insert({ "NumPS", "number of PS shapes" });
		g_parameter_abbreviation.insert({ "NumPub", "number of public variables" });
		g_parameter_abbreviation.insert({ "OA", "objective accuracy" });
		g_parameter_abbreviation.insert({ "OLD", "overlapDgre" });
		g_parameter_abbreviation.insert({ "ONP", "number of obj subspaces" });
		g_parameter_abbreviation.insert({ "PC", "peakCenter" });
		g_parameter_abbreviation.insert({ "PF", "predicFlag" });
		g_parameter_abbreviation.insert({ "PFtype", "global PF shape type" });
		g_parameter_abbreviation.insert({ "ph", "phi" });
		g_parameter_abbreviation.insert({ "PID", "proId" });
		g_parameter_abbreviation.insert({ "PIM", "populationInitialMethod" });
		g_parameter_abbreviation.insert({ "PkS", "peakShape" });
		g_parameter_abbreviation.insert({ "PN", "problem name" });
		g_parameter_abbreviation.insert({ "PNCM", "peakNumChangeMode" });
		g_parameter_abbreviation.insert({ "POS", "peakOffset" });
		g_parameter_abbreviation.insert({ "PPB", "peaksPerBox" });
		g_parameter_abbreviation.insert({ "PPS", "potential prediction strategy" });
		g_parameter_abbreviation.insert({ "PriDimNorm", "private dim norm" });
		g_parameter_abbreviation.insert({ "PriPeakRatio", "private peak height ratio" });
		g_parameter_abbreviation.insert({ "PriSlopeRange", "private peak slope range" });
		g_parameter_abbreviation.insert({ "PRS", "PS radiate slope" });
		g_parameter_abbreviation.insert({ "PS", "population size" });
		g_parameter_abbreviation.insert({ "PSDR", "PS decay rate" });
		g_parameter_abbreviation.insert({ "PSspan", "vector of global PS span" });
		g_parameter_abbreviation.insert({ "PStype", "global PS shape type" });
		g_parameter_abbreviation.insert({ "PT", "proTag" });
		g_parameter_abbreviation.insert({ "PubDimNorm", "public dim norm" });
		g_parameter_abbreviation.insert({ "ra", "radius" });
		g_parameter_abbreviation.insert({ "RBP", "resource4BestPop" });
		g_parameter_abbreviation.insert({ "RID", "runId" });
		g_parameter_abbreviation.insert({ "RP", "repulsion probability" });
		g_parameter_abbreviation.insert({ "RPS", "radiusPS" });
		g_parameter_abbreviation.insert({ "RR", "replaceRatio" });
		g_parameter_abbreviation.insert({ "RS", "repulsion size" });
		g_parameter_abbreviation.insert({ "SB", "sample in basin" });
		g_parameter_abbreviation.insert({ "SBS", "select by subspace" });
		g_parameter_abbreviation.insert({ "SF", "sample frequency" });
		g_parameter_abbreviation.insert({ "SI", "stepIndi" });
		g_parameter_abbreviation.insert({ "SL", "shiftLength" });
		g_parameter_abbreviation.insert({ "SPS", "subpopulation size" });
		g_parameter_abbreviation.insert({ "SVM", "solutionValidationMode" });
		g_parameter_abbreviation.insert({ "THD", "threshold number of solutions" });
		g_parameter_abbreviation.insert({ "TI", "testItems" });
		g_parameter_abbreviation.insert({ "TLS", "timelinkageSeverity" });
		g_parameter_abbreviation.insert({ "TT", "trainingTime" });
		g_parameter_abbreviation.insert({ "TW", "time window" });
		g_parameter_abbreviation.insert({ "UAM", "use acceleration mode" });
		g_parameter_abbreviation.insert({ "UEHO", "use exploratory history only" });
		g_parameter_abbreviation.insert({ "ULP", "use lstm proportion" });
		g_parameter_abbreviation.insert({ "UseLSTM", "use LSTM" });
		g_parameter_abbreviation.insert({ "USPL", "updateSchemeProbabilityLearning" });
		g_parameter_abbreviation.insert({ "UVP", "useVariablePartition" });
		g_parameter_abbreviation.insert({ "VRa", "validRadius" });
		g_parameter_abbreviation.insert({ "VRe", "variableRelation" });
		g_parameter_abbreviation.insert({ "W", "weight" });
		g_parameter_abbreviation.insert({ "XP", "xoverProbability" });
	}

	std::map<std::string, std::string> abbreviationToParameter(const std::string &abbreviations) {
		std::map<std::string, std::string> parameters;
		std::string remove = "\r\t\n\b\v";
		std::stringstream ss_argvs(abbreviations);
		std::string argv;
		while (std::getline(ss_argvs, argv, ' ')) {
			while (size_t pos = argv.find_first_of(remove)) {
				if (std::string::npos == pos) break;
				argv.erase(argv.begin() + pos);
			}
			size_t pos = argv.find('=');
			if (pos == std::string::npos) {
				continue;
			}
			std::string value = argv.substr(pos + 1, argv.size() - 1), name = argv.substr(0, pos);
			if (g_parameter_abbreviation.count(name) == 0) {
				continue;
			}
			parameters[g_parameter_abbreviation.at(name)] = value;
		}
		return parameters;
	}

	void omitParameterInFileName() {
		g_parameters_omit.insert("number of threads");
		g_parameters_omit.insert("sample frequency");
		g_parameters_omit.insert("generation type");
	}

	std::string parametersToFileName(const std::map<std::string, std::string> &parameters) {
		std::stringstream file_name;
		file_name << ofec::g_working_directory << "/result/";
		if (!std::filesystem::exists(file_name.str())) {
			try {
				std::filesystem::create_directory(file_name.str());
			}
			catch (std::filesystem::filesystem_error &e) {
				throw ofec::Exception(e.what());
			}
		}
		std::vector<std::string> abbreviations;
		for (auto &i : parameters) {
			if (g_parameters_omit.count(i.first) == 0) {
				for (auto &j : g_parameter_abbreviation) {
					if (i.first == j.second) {
						abbreviations.push_back(j.first);
						break;
					}
				}
			}
		}
		std::sort(abbreviations.begin(), abbreviations.end());
		for (const auto &abbr : abbreviations) {
			for (auto &i : parameters) {
				if (i.first == g_parameter_abbreviation.at(abbr)) {
					auto str = i.second;
					for (auto &l : str) {
						if (l == '/' || l == '\\') {
							l = '_';
						}
					}
					file_name << abbr << "(" << str << ")__";
				}
			}
		}
		return file_name.str();
	}

	void runAlgorithms(std::list<std::shared_ptr<ofec::Environment>> &environments) {
		for (auto &environment : environments) {
			//{
			//	std::unique_lock<std::mutex> lock(g_cout_mutex);
			//	std::cout << "Algorithm (environment address: " << std::right << std::setw(2) << environment.get() << ") started." << std::endl;
			//}
			try {
				environment->runAlgorithm();
			}
			catch (const ofec::Exception &e) {
				std::unique_lock<std::mutex> lock(g_cout_mutex);
				std::cout << "Algorithm (environment address: "
					<< std::right << std::setw(2) << environment.get() << ") throw except:\n  "
					<< e.what() << std::endl;
				return;
			}
			//{
			//	std::unique_lock<std::mutex> lock(g_cout_mutex);
			//	std::cout << "Algorithm (environment address: "
			//		<< std::right << std::setw(2) << environment.get() << ") terminated. "
			//		<< std::right << std::setw(7) << environment->algorithm()->evaluations() << " FEs costed." << std::endl;
			//}
		}
	}
}

#endif // !OFEC_STATISTICS_CUSTOM_METHOD_HPP
