//#include "custom_method.h"
#include "prime_method.h"
#include "../core/global.h"
#include "../instance/record/rcr_vec_real.h"
#include <fstream>
#include <list>
#include <iostream>
#include <chrono>
#include <thread>
//#include <nbnFigureOut.h>
#include "../utility/myexcept.h"
#include <iostream>
//#include "testOFECmatlab.h"
//#include "drawMeshFun.h"



#ifdef  OFEC_PYTHON
#include "../utility/python/python_caller.h"
#include "Python.h"
#endif //  OFEC_PYTHON

#include "../utility/nbn_visualization/nbn_grid_tree_division.h"
#include "../utility/nbn_visualization/nbn_visualization_data.h"

#include "../utility/nbn_visualization/tree_graph.h"
#include "../utility/nbn_visualization/tree_graph_simple.h"
#include "../utility/general_multithread/general_multithread.h"
#include "../core/problem/uncertainty/dynamic.h"
#include "../instance/algorithm/combination/eax_tsp/environment.h"
#include "../utility/local_optima_network/lon.h"
#include <filesystem>
#include "../utility/local_optima_network/lon.h"
#include "../utility/fitness_landscape_analysis/walk/random_walk.h"
#include "../utility/nbn_visualization/nbn_onemax_simple.h"
#include "../utility/function/custom_function.h"

#include "../utility/functional.h"

#include "../utility/nbn_visualization/nbn_fla.h"

#include "../instance/problem/continuous/single_objective/global/classical/weierstrass.h"


#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../instance/problem/continuous/free_peaks/param_creator.h"
#include "../instance/problem/continuous/free_peaks/parameter_function.h"

#include "../instance/algorithm/visualize/LON/sampling/lon_con_sampling.h"
#include "../instance/algorithm/visualize/LON/experiment/lon.h"
#include "../instance/problem/combination/wmodel/wmodel.h"
#include "../core/problem/continuous/continuous.h"

#include "../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include "../instance/problem/combination/travelling_salesman/tsp_generator/tsp_benchmark_generator.h"







void sampleSolsAroundTSPThreadTask(
	const ofec::TravellingSalesman::SolType& sol,
	std::vector<std::shared_ptr<ofec::SolBase>>& samples,
	int neighbor_k,
	int from, int to,
	ofec::Problem* pro,
	ofec::Random* rnd
) {

	int a(0), b(0);
	int dim = pro->numVariables();

	std::shared_ptr<ofec::SolBase> cursol;
	for (int idx(from); idx < to; ++idx) {
		cursol.reset(pro->createSolution(sol));
		auto& cursolt = dynamic_cast<ofec::TravellingSalesman::SolType&>(*cursol);
		for (int idy(0); idy < neighbor_k; ++idy) {
			a = rnd->uniform.nextNonStd<int>(0, dim);
			b = rnd->uniform.nextNonStd<int>(0, dim);
			std::swap(cursolt.variable()[a], cursolt.variable()[b]);
		}

		samples[idx] = cursol;

	}
}



void sampleSolsAroundTSPMultiThread(
	const ofec::TravellingSalesman::SolType& sol,
	std::vector<std::shared_ptr<ofec::SolBase>>& samples,
	int neighbor_k,
	int numSamples,
	ofec::Problem* pro,
	ofec::Random* rnd
) {
	samples.resize(numSamples);


	std::cout << "generate solutions by multithread" << std::endl;
	int num_task = std::thread::hardware_concurrency();
	int num_samples = numSamples;
	std::vector<std::thread> thrds;
	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);


	std::vector<std::shared_ptr<ofec::Random>> rnds(num_task);
	for (auto& it : rnds) {
		double randomSeed(rnd->uniform.nextNonStd<double>(0.01, 1.0));
		it.reset(new ofec::Random(randomSeed));
	}

	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];

		thrds.push_back(std::thread(
			sampleSolsAroundTSPThreadTask,
			std::cref(sol), std::ref(samples), neighbor_k,
			tasks[i], tasks[i + 1], pro, rnds[i].get()));
	}
	for (auto& thrd : thrds)
		thrd.join();
}







void sampleSolsRandomThreadTask(
	std::vector<std::shared_ptr<ofec::SolBase>>& samples,
	int from, int to,
	ofec::Problem* pro,
	ofec::Random* rnd
) {

	int dim = pro->numVariables();

	std::shared_ptr<ofec::SolBase> cursol;
	for (int idx(from); idx < to; ++idx) {
		cursol.reset(pro->createSolution());
		cursol->initialize(pro, rnd);
		//	auto& cursolt = dynamic_cast<ofec::TravellingSalesman::SolType&>(*cursol);
		samples[idx] = cursol;
		cursol->evaluate(pro);
	}
}



void sampleSolsRandomMultiThread(
	std::vector<std::shared_ptr<ofec::SolBase>>& samples,
	int numSamples,
	ofec::Problem* pro,
	ofec::Random* rnd
) {
	samples.resize(numSamples);


	std::cout << "generate solutions by multithread" << std::endl;
	int num_task = std::thread::hardware_concurrency();
	int num_samples = numSamples;
	std::vector<std::thread> thrds;
	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);


	std::vector<std::shared_ptr<ofec::Random>> rnds(num_task);
	for (auto& it : rnds) {
		double randomSeed(rnd->uniform.nextNonStd<double>(0.01, 1.0));
		it.reset(new ofec::Random(randomSeed));
	}

	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];

		thrds.push_back(std::thread(
			sampleSolsRandomThreadTask,
			std::ref(samples),
			tasks[i], tasks[i + 1], pro, rnds[i].get()));
	}
	for (auto& thrd : thrds)
		thrd.join();
}



struct NBNAlgInfo {
	int m_iteration = 0;
	int m_numLoop = 0;
	double m_maxChange = 0;
	double m_dev = 0;
	double m_millisecond = 0;

	static void outputHead(std::ostream& out) {
		out << "iteration\tloop\tmax\tdev\tseconds" << std::endl;
	}
	void output(std::ostream& out) {
		out << m_iteration << "\t" << m_numLoop << "\t" << m_maxChange << "\t" << m_dev << "\t" << m_millisecond << std::endl;
	}

	void updateInfo(const std::vector<double>& dis1, const std::vector<double>& dis2) {
		m_dev = 0;
		m_maxChange = 0;
		for (int idx(0); idx < dis1.size(); ++idx) {
			m_maxChange = std::max(m_maxChange, std::abs(dis1[idx] - dis2[idx]));
			m_dev += (dis1[idx] - dis2[idx]) * (dis1[idx] - dis2[idx]);
		}

		m_dev = sqrt(m_dev / dis1.size());

	}
};







void solDataToNBN(
	const std::string& outputDir,
	const std::vector<int>& solId,
	const std::vector<std::shared_ptr<ofec::SolBase>>& sols,
	const std::string& tspname,
	const std::string& filename,
	ofec::Problem* pro,
	ofec::Random* rnd
	//,
	//std::vector<NBNAlgInfo> & algInfo
) {
	using namespace ofec;

	//	algInfo.clear();
	std::cout << "running\t" << tspname << "\t" << filename << std::endl;

	std::string fileDir = outputDir + "/" + tspname;
	std::filesystem::create_directories(fileDir);
	std::string filepath = fileDir + "/" + filename;


	//ofec::ParamMap params;
	//params["problem name"] = std::string("TSP");
	//params["dataFile1"] = tspname;
	////	params["division type"] = static_cast<int>(ofec::NBN_VisualizationData::NBN_Division_Type::kEdgeDivision);
	//	//  params["continous sample size"] = int(2e3);
	////params["continous sample size"] = int(2e6 + 2e3);
	////params["init sample size"] = int(2e6);
	//params["add optimal"] = false;


	//auto param= ADD_PARAM(params);
	//auto pro = ADD_PRO(param, 0.1);
	//pro->initialize();

	//auto rnd = std::make_shared<Random>(0.5);


	ofec::NBN_DivisionData::initialize_static_info();
	ofec::NBN_DivisionData nbn_network;
	nbn_network.m_division.reset(new ofec::NBN_EdgeMultiThreadDivision);
	ofec::NBN_VisualizationData nbn_data;


	auto eval_fun =
		[](ofec::SolBase& sol, Problem* pro) {
		using namespace ofec;
		sol.evaluate(pro);
		ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
	};
	{
		auto& division(nbn_network.m_division);
		division->setMaxSampleSize(0);
		division->setMultiThread(true);
		int numThread = std::thread::hardware_concurrency();
		division->setNumberThread(numThread);
		std::cout << "number of thread" << numThread << std::endl;
		division->initialize(pro->getSharedPtr(), rnd->getSharedPtr(), eval_fun, true);
		auto& edge_division = dynamic_cast<ofec::NBN_EdgeMultiThreadDivision&>(*division);
		int numLoop = 1e3 + 50;
		//edge_division.setNumberLoop(std::ceil(double(numLoop) / double(numThread)));
		edge_division.setNumberLoop(numLoop);
		edge_division.updateSols(sols);
		//int from_size = sols.size();
		//edge_division.generateSols(sols.size() + 1e3);
		//int to_size = edge_division.size();
		//edge_division.filterSameSolutions(from_size, to_size);

		auto start = std::chrono::system_clock::now();
		edge_division.updateNetwork();
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<float> difference = end - start;
		double milliseconds = difference.count();


		std::vector<double> beforeDis2parent(sols.size(), 1);
		std::vector<int> belong;
		std::vector<double> dis2parent;
		edge_division.getResult(belong, dis2parent);

		std::ofstream algInfoOut(filepath + ".algInfo");
		NBNAlgInfo algInfo;
		NBNAlgInfo::outputHead(algInfoOut);



		algInfo.m_iteration = 0;
		algInfo.m_numLoop = numThread;
		algInfo.m_millisecond = milliseconds;
		algInfo.updateInfo(beforeDis2parent, dis2parent);
		swap(dis2parent, beforeDis2parent);
		algInfo.output(algInfoOut);
		algInfo.output(std::cout);

		/*while (true)*/ {
			edge_division.updateDivisionMultithread();
			end = std::chrono::system_clock::now();
			std::chrono::duration<float> difference = end - start;
			milliseconds = difference.count();
			edge_division.getResult(belong, dis2parent);

			algInfo.m_iteration++;
			algInfo.m_millisecond = milliseconds;
			algInfo.updateInfo(beforeDis2parent, dis2parent);
			swap(dis2parent, beforeDis2parent);
			algInfo.output(algInfoOut);
			algInfo.output(std::cout);

			//		if (algInfo.m_maxChange < 0.1 && algInfo.m_dev <= 0.002) break;
		}





		std::ofstream nbnOut(filepath + ".nbn");
		std::vector<double> fitness(sols.size());
		for (int idx(0); idx < fitness.size(); ++idx) {
			fitness[idx] = sols[idx]->fitness();
		}
		outputNBNdata(nbnOut, solId, belong, dis2parent, fitness);
		nbnOut.close();


	}
}








void generateTSPdata() {

	using namespace ofec;

	for (int idx(1); idx <= 5; ++idx) {
		int n = idx * 100;
		int rand = 25;
		{
			std::string file1 = "E" + std::to_string(n) + "-" + std::to_string(rand) + ".tsp";
			{
				PortGen tsp_gen;
				tsp_gen.init();
				std::ofstream out(file1);
				tsp_gen.generate(n, 10000 + rand, out);
				out.close();
			}
		}



		{
			std::string file1 = "C" + std::to_string(n) + "-" + std::to_string(rand) + ".tsp";
			{
				PortCGen tsp_gen;
				tsp_gen.init();
				std::ofstream out(file1);
				tsp_gen.generate(n, 10000 + rand, out);
				out.close();
			}
		}
	}

}


void testPro(const std::string& tspname) {
	using namespace ofec;
	ofec::ParamMap params;
	params["problem name"] = std::string("TSP");
	params["dataFile1"] = tspname;
	//	params["division type"] = static_cast<int>(ofec::NBN_VisualizationData::NBN_Division_Type::kEdgeDivision);
		//  params["continous sample size"] = int(2e3);
	//params["continous sample size"] = int(2e6 + 2e3);
	//params["init sample size"] = int(2e6);
	params["add optimal"] = false;

	auto param = ADD_PARAM(params);
	auto pro = ADD_PRO(param, 0.1);
	pro->initialize();
}



#include <filesystem>
#include <string>
#include <vector>
#include <iostream>

namespace ofec {

	void registerParamAbbr();
	void customizeFileName();
	void run() {
		registerParamAbbr();
		customizeFileName();
		using namespace std;



		auto outputDir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/nbn_visualization_tsp2/";
		std::string filename;
		//generateTSPdata();


		std::vector<std::string> tsp_names = { "BCL380" };
		for (int idx(5); idx >= 1; --idx) {
			int n = idx * 100;
			int rand = 25;
			{
				std::string file1 = "E" + std::to_string(n) + "-" + std::to_string(rand);
				tsp_names.push_back(file1);
			}
			{
				std::string file1 = "C" + std::to_string(n) + "-" + std::to_string(rand);
				tsp_names.push_back(file1);
			}
		}



		std::vector<std::string> tsp_names2 = { "GR666","PCB442","XQL662" };

		for (auto& it : tsp_names2) {
			tsp_names.push_back(it);
		}


		//for (auto& it : tsp_names) {
		//	testPro(it);
		//}


		int numSample = 3e6 + 100;
		std::vector<int> solId(numSample);
		for (int idx(0); idx < numSample; ++idx) {
			solId[idx] = idx;
		}
		std::vector<std::shared_ptr<ofec::SolBase>> samples;

		for (auto& tspname : tsp_names) {

			ofec::ParamMap params;
			params["problem name"] = std::string("TSP");
			params["dataFile1"] = tspname;
			params["add optimal"] = false;


			auto param = ADD_PARAM(params);
			auto pro = ADD_PRO(param, 0.1);
			pro->initialize();
			auto rnd = std::make_shared<Random>(0.5);


			sampleSolsRandomMultiThread(
				samples, numSample, pro.get(), rnd.get());

			std::shared_ptr<ofec::SolBase> best_sol(samples.front());
			for (auto& it : samples) {
				if (it->objective()[0] < best_sol->objective()[0]) {
					best_sol = it;
				}
			}

			filename = "tspname_" + tspname + "_origin_";

			solDataToNBN(outputDir, solId, samples, tspname,
				filename, pro.get(), rnd.get());


			int neighbor_k = pro->numVariables() / 2;
			while (neighbor_k >= 10) {
				sampleSolsAroundTSPMultiThread(*best_sol, samples, neighbor_k, numSample, pro.get(), rnd.get());
				filename = "tspname_" + tspname + "_neighbor_" + std::to_string(neighbor_k);
				solDataToNBN(outputDir, solId, samples, tspname,
					filename, pro.get(), rnd.get());
				neighbor_k /= 4.0;
			}


		}






		//runOneMaxTotalTaskAnalysis();

		//std::vector<std::string> filenames;
		//std::string dir = "//172.29.65.56/share/Student/2018/YiyaDiao/NBN_data/visualization_nbn_3/";

		//readDir(dir, filenames);

	//	calMultiTask();


		//calMultiTask();
		//analyseProTask();

		//filterMMOP_flacco();

		//filter_flacco_FreePeakBasinAtraction();

	}

	void registerParamAbbr() {
		g_param_abbr.insert({ "AID", "algId" });
		g_param_abbr.insert({ "al", "alpha" });
		g_param_abbr.insert({ "AN", "algorithm name" });
		g_param_abbr.insert({ "be", "beta" });
		g_param_abbr.insert({ "BPF", "mean basin potential" });
		g_param_abbr.insert({ "C1", "accelerator1" });
		g_param_abbr.insert({ "C2", "accelerator2" });
		g_param_abbr.insert({ "ca", "case" });
		g_param_abbr.insert({ "CAN", "cluster add neighbors" });
		g_param_abbr.insert({ "CDBGFID", "comDBGFunID" });
		g_param_abbr.insert({ "CE", "crossover eta" });
		g_param_abbr.insert({ "CF", "changeFre" });
		g_param_abbr.insert({ "CoF", "convFactor" });
		g_param_abbr.insert({ "CPR", "changePeakRatio" });
		g_param_abbr.insert({ "CR", "crossover rate" });
		g_param_abbr.insert({ "CS", "changeSeverity" });
		g_param_abbr.insert({ "CT", "changeType" });
		g_param_abbr.insert({ "CT2", "changeType2" });
		g_param_abbr.insert({ "CTH", "convThreshold" });
		g_param_abbr.insert({ "DD", "diversity degree" });
		g_param_abbr.insert({ "DD1", "dataDirectory1" });
		g_param_abbr.insert({ "DF1", "dataFile1" });
		g_param_abbr.insert({ "DF2", "dataFile2" });
		g_param_abbr.insert({ "DF3", "dataFile3" });
		g_param_abbr.insert({ "DFR", "dataFileResult" });
		g_param_abbr.insert({ "DG", "durationGeneration" });
		g_param_abbr.insert({ "DM", "divisionMode" });
		g_param_abbr.insert({ "DR", "detectRatio" });
		g_param_abbr.insert({ "EBP", "evolve by potential" });
		g_param_abbr.insert({ "ECF", "evalCountFlag" });
		g_param_abbr.insert({ "ep", "epsilon" });
		g_param_abbr.insert({ "ER", "exlRadius" });
		g_param_abbr.insert({ "et", "eta" });
		g_param_abbr.insert({ "exitPS", "exploit pop size" });
		g_param_abbr.insert({ "exrePS", "explore pop size" });
		g_param_abbr.insert({ "F", "scaling factor" });
		g_param_abbr.insert({ "FA", "flagAsymmetry" });
		g_param_abbr.insert({ "FI", "flagIrregular" });
		g_param_abbr.insert({ "FN", "flagNoise" });
		g_param_abbr.insert({ "FNDC", "flagNumDimChange" });
		g_param_abbr.insert({ "FNPC", "flagNumPeakChange" });
		g_param_abbr.insert({ "FR", "flagRotation" });
		g_param_abbr.insert({ "FTL", "flagTimeLinkage" });
		g_param_abbr.insert({ "ga", "gamma" });
		g_param_abbr.insert({ "gOptFlag", "gOptFlag" });
		g_param_abbr.insert({ "GS", "glstrture" });
		g_param_abbr.insert({ "GT", "generation type" });
		g_param_abbr.insert({ "HCM", "heightConfigMode" });
		g_param_abbr.insert({ "HR", "hibernatingRadius" });
		g_param_abbr.insert({ "IC", "initial clustering" });
		g_param_abbr.insert({ "IT1", "interTest1" });
		g_param_abbr.insert({ "IT2", "interTest2" });
		g_param_abbr.insert({ "IT3", "interTest3" });
		g_param_abbr.insert({ "IT4", "interTest4" });
		g_param_abbr.insert({ "ITV", "interval number of solutions" });
		g_param_abbr.insert({ "JH", "jumpHeight" });
		g_param_abbr.insert({ "KF", "kFactor" });
		g_param_abbr.insert({ "LF", "lFactor" });
		g_param_abbr.insert({ "ME", "maximum evaluations" });
		g_param_abbr.insert({ "MI", "maxIter" });
		g_param_abbr.insert({ "MNPS", "minNumPopSize" });
		g_param_abbr.insert({ "MP", "mutProbability" });
		g_param_abbr.insert({ "MR", "mutation rate" });
		g_param_abbr.insert({ "MRT", "maxRunTime" });
		g_param_abbr.insert({ "MSDE", "mutation strategy" });
		g_param_abbr.insert({ "MSI", "maxSucIter" });
		g_param_abbr.insert({ "MSS", "max subpop size" });
		g_param_abbr.insert({ "MuE", "mutation eta" });
		g_param_abbr.insert({ "NAS", "number of atomspaces" });
		g_param_abbr.insert({ "NB", "numBox" });
		g_param_abbr.insert({ "NC", "numChange" });
		g_param_abbr.insert({ "NCus", "numCus" });
		g_param_abbr.insert({ "ND", "number of variables" });
		g_param_abbr.insert({ "NER", "number of solutions for exploration" });
		g_param_abbr.insert({ "NF", "noiseFlag" });
		g_param_abbr.insert({ "NGO", "numGOpt" });
		g_param_abbr.insert({ "NNCus", "numNewCus" });
		g_param_abbr.insert({ "NO", "number of objectives" });
		g_param_abbr.insert({ "NoN", "number of neighbors" });
		g_param_abbr.insert({ "NP", "numPeak" });
		g_param_abbr.insert({ "NPa", "numPartition" });
		g_param_abbr.insert({ "NPR", "numParetoRegion" });
		g_param_abbr.insert({ "NPS", "numPS" });
		g_param_abbr.insert({ "NR", "number of runs" });
		g_param_abbr.insert({ "NS", "noiseSeverity" });
		g_param_abbr.insert({ "NSP", "number of subpopulations" });
		g_param_abbr.insert({ "NSSP", "number of subspaces" });
		g_param_abbr.insert({ "NT", "numTask" });
		g_param_abbr.insert({ "NumGPS", "number of global PS shapes" });
		g_param_abbr.insert({ "NumPri", "vector of private variables" });
		g_param_abbr.insert({ "NumPriPeak", "vector of private peaks" });
		g_param_abbr.insert({ "NumPS", "number of PS shapes" });
		g_param_abbr.insert({ "NumPub", "number of public variables" });
		g_param_abbr.insert({ "OA", "objective accuracy" });
		g_param_abbr.insert({ "OLD", "overlapDgre" });
		g_param_abbr.insert({ "ONP", "number of obj subspaces" });
		g_param_abbr.insert({ "PC", "peakCenter" });
		g_param_abbr.insert({ "PF", "predicFlag" });
		g_param_abbr.insert({ "PFtype", "global PF shape type" });
		g_param_abbr.insert({ "ph", "phi" });
		g_param_abbr.insert({ "PID", "proId" });
		g_param_abbr.insert({ "PIM", "populationInitialMethod" });
		g_param_abbr.insert({ "PkS", "peakShape" });
		g_param_abbr.insert({ "PN", "problem name" });
		g_param_abbr.insert({ "PNCM", "peakNumChangeMode" });
		g_param_abbr.insert({ "POS", "peakOffset" });
		g_param_abbr.insert({ "PPB", "peaksPerBox" });
		g_param_abbr.insert({ "PPS", "potential prediction strategy" });
		g_param_abbr.insert({ "PriDimNorm", "private dim norm" });
		g_param_abbr.insert({ "PriPeakRatio", "private peak height ratio" });
		g_param_abbr.insert({ "PriSlopeRange", "private peak slope range" });
		g_param_abbr.insert({ "PRS", "PS radiate slope" });
		g_param_abbr.insert({ "PS", "population size" });
		g_param_abbr.insert({ "PSDR", "PS decay rate" });
		g_param_abbr.insert({ "PSspan", "vector of global PS span" });
		g_param_abbr.insert({ "PStype", "global PS shape type" });
		g_param_abbr.insert({ "PT", "proTag" });
		g_param_abbr.insert({ "PubDimNorm", "public dim norm" });
		g_param_abbr.insert({ "ra", "radius" });
		g_param_abbr.insert({ "RBP", "resource4BestPop" });
		g_param_abbr.insert({ "RID", "runId" });
		g_param_abbr.insert({ "RP", "reject probability" });
		g_param_abbr.insert({ "RPS", "radiusPS" });
		g_param_abbr.insert({ "RR", "replaceRatio" });
		g_param_abbr.insert({ "RS", "reject scale" });
		g_param_abbr.insert({ "SB", "sample in basin" });
		g_param_abbr.insert({ "SBS", "select by subspace" });
		g_param_abbr.insert({ "SF", "sample frequency" });
		g_param_abbr.insert({ "SI", "stepIndi" });
		g_param_abbr.insert({ "SL", "shiftLength" });
		g_param_abbr.insert({ "SPS", "subpopulation size" });
		g_param_abbr.insert({ "SVM", "solutionValidationMode" });
		g_param_abbr.insert({ "THD", "threshold number of solutions" });
		g_param_abbr.insert({ "TI", "testItems" });
		g_param_abbr.insert({ "TLS", "timelinkageSeverity" });
		g_param_abbr.insert({ "TT", "trainingTime" });
		g_param_abbr.insert({ "TW", "time window" });
		g_param_abbr.insert({ "UAM", "use acceleration mode" });
		g_param_abbr.insert({ "UEHO", "use exploratory history only" });
		g_param_abbr.insert({ "ULP", "use lstm proportion" });
		g_param_abbr.insert({ "UseLSTM", "use LSTM" });
		g_param_abbr.insert({ "USPL", "updateSchemeProbabilityLearning" });
		g_param_abbr.insert({ "UVP", "useVariablePartition" });
		g_param_abbr.insert({ "VRa", "validRadius" });
		g_param_abbr.insert({ "VRe", "variableRelation" });
		g_param_abbr.insert({ "W", "weight" });
		g_param_abbr.insert({ "XP", "xoverProbability" });
	}

	void customizeFileName() {
		g_param_omit.insert("algId");
		g_param_omit.insert("dataDirectory1");
		g_param_omit.insert("evalCountFlag");
		g_param_omit.insert("flagNoise");
		g_param_omit.insert("flagNumPeakChange");
		g_param_omit.insert("flagTimeLinkage");
		g_param_omit.insert("gOptFlag");
		g_param_omit.insert("hibernatingRadius");
		g_param_omit.insert("minNumPopSize");
		g_param_omit.insert("numGOpt");
		g_param_omit.insert("number of runs");
		g_param_omit.insert("numTask");
		g_param_omit.insert("peakNumChangeMode");
		g_param_omit.insert("proId");
		g_param_omit.insert("sample frequency");
		g_param_omit.insert("solutionValidationMode");
		g_param_omit.insert("dataFile1");
		g_param_omit.insert("generation type");
	}
}
