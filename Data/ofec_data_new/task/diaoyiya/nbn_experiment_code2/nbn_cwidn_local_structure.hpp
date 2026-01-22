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
#include <filesystem>
#include "../utility/nbn_visualization/nbn_onemax_simple.h"
#include "../utility/function/custom_function.h"
#include "../utility/functional.h"
#include "../instance/problem/realworld/csiwdn/csiwdn.h"
#include <shared_mutex>

#include <filesystem>
#include <string>
#include <vector>
#include <iostream>

#include "../utility/nbn_visualization/nbn_grid_tree_division_allSize.h"






class ThreadInfoNBN :public utility::GeneralMultiThreadInfo {
public:
	//std::ofstream& m_errOut;
	//std::mutex& m_errOut_mtx;
	std::string m_filename;
	std::string m_readfilename;
	std::string m_taskname;

	int m_dim = 0;
	std::string m_proname;
	int m_curRunId = 0;
	std::string m_save_dir;
	bool m_flagOpt = false;



	// for wmodel

	int m_instance_id = 0;
	int m_function_id = 0;
	int m_sample_size = 0;

	int m_curSampleSize = 0;



	double m_dummy_select_rate = 0.0;
	int m_epistasis_block_size = 0;
	int m_neutrality_mu = 0;
	int m_ruggedness_gamma = 0;



	std::string getConTaskName() {
		return m_proname + "_Dim_" + std::to_string(m_dim);
	}


	std::string getConTaskSampleSizeName() {
		return m_proname + "_Dim_" + std::to_string(m_dim) + "_Sample_" + std::to_string(m_curSampleSize);

	}

	std::string getBBOBfilename() {
		return std::string("BBOB_IOH") + "_Dim_" + std::to_string(m_dim)
			+ "_SampleSize_" + std::to_string(m_sample_size)
			+ "_FunctionId_" + std::to_string(m_function_id)
			+ "_InstanceId_" + std::to_string(m_instance_id)
			;
	}
};
class ThreadGlobalNBN :public utility::GeneralMultiThreadInfo {
public:


	std::string m_read_dir;
	std::string m_save_dir;
	std::string m_nbn_subfix = "_nbn.txt";
	int m_totalRun = 0;
	//	std::string m_save_path;
	int m_maxSample = 2e6;


	std::ofstream m_out;
	std::mutex m_out_mtx;

};




struct NBNdataVecInt {
	std::vector<std::shared_ptr<ofec::SolBase>> m_sols;

	std::vector<int> m_solIds;
	std::vector<int> m_belongs;
	std::vector<double> m_dis2parent;
	std::vector<double> m_fitness;

	int m_bestSolId = -1;


	std::vector<std::vector<int>> m_sol_vecInt;

	void initialize(int size) {
		m_sols.resize(size);
		m_solIds.resize(size);
		m_belongs.resize(size);
		m_dis2parent.resize(size);
		m_fitness.resize(size);


		m_sol_vecInt.resize(size);

		m_bestSolId = -1;

		for (int idx(0); idx < m_solIds.size(); ++idx) {
			m_solIds[idx] = idx;
			m_belongs[idx] = idx;
			m_dis2parent[idx] = std::numeric_limits<double>::max();
			m_fitness[idx] = -std::numeric_limits<double>::max();
		}
	}

	void assignSolMem(ofec::Problem* pro) {
		for (auto& it : m_sols) {
			it.reset(pro->createSolution());
		}
	}



	bool compareTwoSol(int a, int b) {

		if (m_sols[a]->fitness() == m_sols[b]->fitness()) {
			{
				bool flag = false;
				for (int idx(0); idx < m_sol_vecInt[a].size(); ++idx) {
					if (m_sol_vecInt[a][idx] != m_sol_vecInt[b][idx]) {
						return m_sol_vecInt[a][idx] < m_sol_vecInt[b][idx];
					}
				}
				return m_solIds[a] < m_solIds[b];
			}
		}
		else return m_sols[a]->fitness() < m_sols[b]->fitness();
	}


	void outputNBN(std::ostream& out) {
		for (int idx(0); idx < m_solIds.size(); ++idx) {
			out << m_solIds[idx] << "\t" << m_belongs[idx] << "\t" << m_dis2parent[idx] << "\t" << m_fitness[idx] << std::endl;
		}
	}

	void outputSol_CWIDN(std::ostream& out,
		ofec::Problem* pro,
		const ofec::SolBase& sol) {
		using namespace ofec;

		auto& var = dynamic_cast<const Solution<VarCSIWDN> &>(sol).variable();
		out << sol.fitness() << "\t";

		out << std::endl;
		for (size_t z = 0; z < CAST_CSIWDN(pro)->numberSource(); z++) {
			out << var.index(z) << "\t" << var.multiplier(z).size() << "\t";
			for (auto& j : var.multiplier(z)) {
				out << j << "\t";
				//j = rnd->uniform.nextNonStd<float>(m_min_multiplier, m_max_multiplier);
			}
			out << std::endl;
		}
	}


	void outputSols_CSWIDN(std::ofstream& out, ofec::Problem* pro) {
		for (int idx(0); idx < m_solIds.size(); ++idx) {
			out << m_solIds[idx] << "\t";
			outputSol_CWIDN(out, pro, *m_sols[idx]);
		}
	}

	void updateFitness() {
		for (int idx(0); idx < m_solIds.size(); ++idx) {
			m_fitness[idx] = m_sols[idx]->fitness();
		}
	}


};


// global structure

void updateRandSol(
	ofec::Problem* pro, ofec::Random* rnd,
	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)>& eval_fun,
	int subsols,
	int from, int to,
	std::vector<std::shared_ptr<ofec::SolBase>>& sols) {
	int id(0);
	for (int idx(from); idx < to; ++idx) {
		id = idx * subsols;
		sols[id].reset(pro->createSolution());
		sols[id]->initialize(pro, rnd);
		//sols[id]->setFitness(rnd->uniform.next());

		eval_fun(*sols[id], pro);

	}
}



void updateRandSolSubRegion(
	ofec::Problem* pro, ofec::Random* rnd,
	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)>& eval_fun,
	int subsols,
	int from, int to,
	std::vector<std::shared_ptr<ofec::SolBase>>& sols) {
	using namespace ofec;
	int id(0);

	for (int idx(from); idx < to; ++idx) {

		id = std::floor(idx / subsols) * subsols;
		if (id != idx) {
			sols[idx].reset(pro->createSolution());
			auto& centerSol = dynamic_cast<ofec::CSIWDN::solutionType&>(*sols[id]);
			ofec::CSIWDN::solutionType& cursol = dynamic_cast<ofec::CSIWDN::solutionType&>(*sols[idx]);
			cursol.variable() = centerSol.variable();

			CAST_CSIWDN(pro)->initSolutionMultiplier(cursol, rnd);

			//sols[idx]->setFitness(rnd->uniform.next());
			//sols[id]->initialize(pro, rnd);
			eval_fun(*sols[idx], pro);
		}
		//for(int idx())
	}
}


void generateDefaultPro(std::shared_ptr<ofec::Problem>& pro) {

	using namespace ofec;
	ParamMap params;
	params["problem name"] = std::string("CSIWDN");
	params["dataFile1"] = std::string("Net2");
	params["dataFile2"] = std::string("case1");
	params["dataFile3"] = std::string("case11");
	params["use LSTM"] = false;
	auto param = std::shared_ptr<const ParamMap>(new ParamMap(params));
	pro = ofec::ADD_PRO(param, 0.1);
	pro->initialize();
	CAST_CSIWDN(pro.get())->setPhase(CAST_CSIWDN(pro.get())->totalPhase());
}



// generate multiSols
void generateSols(ofec::Problem* pro, ofec::Random* rnd,
	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)>& eval_fun,
	NBNdataVecInt& data,
	int numSols = 2e4,
	int subSols = 100
) {
	using namespace ofec;

	data.initialize(numSols * subSols);
	auto& sols(data.m_sols);

	// generate sols
	//std::vector<std::shared_ptr<ofec::SolBase>> sols;
	//sols.resize(numSols * subSols);

	//	CAST_CSIWDN(pro)->setInitType(ofec::CSIWDN::InitType::kRandom);

	int num_task = std::thread::hardware_concurrency();
	int num_samples = numSols;
	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);

	std::vector<ofec::Random> rnds;
	std::vector<std::shared_ptr<ofec::Problem>> pros(num_task);
	for (int idx(0); idx < num_task; ++idx) {
		rnds.emplace_back(rnd->uniform.next());
		auto& cur_pro = pros[idx];
		generateDefaultPro(cur_pro);
	}
	std::pair<int, int> from_to;
	{
		std::vector<std::thread> thrds;

		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			thrds.push_back(std::thread(
				updateRandSol, pros[i].get(), &rnds[i], std::ref(eval_fun), subSols,
				from_to.first, from_to.second, std::ref(sols)));
		}
		for (auto& thrd : thrds)
			thrd.join();
		thrds.clear();
	}


	num_samples = sols.size();
	UTILITY::assignThreads(num_samples, num_task, tasks);

	{
		std::vector<std::thread> thrds;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			thrds.push_back(std::thread(
				updateRandSolSubRegion, pros[i].get(), &rnds[i], std::ref(eval_fun), subSols,
				from_to.first, from_to.second, std::ref(sols)));
		}
		for (auto& thrd : thrds)
			thrd.join();
		thrds.clear();
	}
	data.updateFitness();

}


void evaluateSolTask(
	ofec::Problem* pro,
	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)>& eval_fun,
	int from, int to,
	std::vector<std::shared_ptr<ofec::SolBase>>& sols) {
	for (int idx(from); idx < to; ++idx) {
		eval_fun(*sols[idx], pro);
	}
}



void evaluateSolMultiTask(
	ofec::Problem* pro,
	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)>& eval_fun,
	std::vector<std::shared_ptr<ofec::SolBase>>& sols) {
	std::vector<int> tasks;
	auto num_samples = sols.size();
	int num_task = std::thread::hardware_concurrency();
	UTILITY::assignThreads(num_samples, num_task, tasks);

	{

		std::vector<std::shared_ptr<ofec::Problem>> pros(num_task);
		for (int idx(0); idx < num_task; ++idx) {
			auto& cur_pro = pros[idx];
			generateDefaultPro(cur_pro);
		}
		std::pair<int, int> from_to;
		std::vector<std::thread> thrds;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			thrds.push_back(std::thread(
				evaluateSolTask, pros[i].get(), std::ref(eval_fun),
				from_to.first, from_to.second, std::ref(sols)));
		}
		for (auto& thrd : thrds)
			thrd.join();
		thrds.clear();
	}
}






void calNBN_iterater(ofec::Problem* pro,
	const std::vector<int>& subArea, int& bestId,
	NBNdataVecInt& nbn_data) {

	if (subArea.empty()) {
		bestId = -1;
		return;
	}
	auto& sols = nbn_data.m_sols;
	auto& belong = nbn_data.m_belongs;
	auto& dis2parent = nbn_data.m_dis2parent;

	std::vector<int> sortIdx(subArea);
	std::sort(sortIdx.begin(), sortIdx.end(), [&](int a, int b) {
		return nbn_data.compareTwoSol(a, b);
	});
	bestId = sortIdx.back();
	for (int idx(0); idx < sortIdx.size(); ++idx) {
		int ida = sortIdx[idx];
		for (int idy(idx + 1); idy < sortIdx.size(); ++idy) {

			int idb = sortIdx[idy];
			/*if (sols[ida]->fitness() < sols[idb]->fitness())*/ {
				double dis = sols[ida]->variableDistance(*sols[idb], pro);
				if (dis < dis2parent[ida]) {
					belong[ida] = idb;
					dis2parent[ida] = dis;
				}
			}
		}
	}

}

void calNBN_iterater_ThreadTask(
	int from, int to,
	ofec::Problem* pro,
	const std::vector<std::vector<int>>& subAreas, std::vector<int>& bestId,
	NBNdataVecInt& nbn_data
) {
	for (int idSub(from); idSub < to; ++idSub) {
		calNBN_iterater(pro, subAreas[idSub], bestId[idSub], nbn_data);
	}
}




void getNodeState(int m, std::vector<std::vector<bool>>& selected, std::vector<bool>& cur, int from) {
	if (from >= cur.size())return;

	cur[from] = true;
	--m;
	if (m == 0) {
		selected.push_back(cur);
	}
	else getNodeState(m, selected, cur, from + 1);

	++m;
	cur[from] = false;
	getNodeState(m, selected, cur, from + 1);
}


int getSolDivisionId(int stateId, int numDivision, int numNode,
	const std::vector<bool>& state, const std::vector<int>& solVecInt) {
	int divisionId = stateId * numDivision;
	int curId = 0;
	//auto& cursol = dynamic_cast<const ofec::CSIWDN::solutionType&>(sol);

	for (int idx(0); idx < state.size(); ++idx) {
		if (state[idx]) {
			//curId = curId * numNode + cursol.variable().index(idx) - 1;
			curId = curId * numNode + solVecInt[idx];
		}
	}
	return divisionId + curId;
}

void calNBN_division(ofec::Problem* pro,
	const std::vector<int>& subArea, int& bestId,
	NBNdataVecInt& nbn_data
) {
	using namespace ofec;

	auto& sols = nbn_data.m_sols;
	auto& belong = nbn_data.m_belongs;
	auto& dis2parent = nbn_data.m_dis2parent;
	auto& solVecInts = nbn_data.m_sol_vecInt;


	int selectedSource = 0;
	int numDivision = 1;
	while (selectedSource < CAST_CSIWDN(pro)->numSource() && numDivision * CAST_CSIWDN(pro)->numberNode() < subArea.size()) {
		selectedSource++;
		numDivision *= CAST_CSIWDN(pro)->numberNode();
	}

	std::vector<bool> curState(CAST_CSIWDN(pro)->numSource(), false);
	std::vector<std::vector<bool>> totalState;
	std::vector<std::vector<int>> divisions;
	std::vector<int> bestSolIds = subArea;
	std::set<int> unique_ids;
	//	std::vector<bool> visited(subArea.size(), false);
	for (int numSelect(selectedSource); numSelect > 0; --numSelect) {


		std::cout << "calNBN division numSelect " << numSelect << std::endl;
		totalState.clear();
		getNodeState(numSelect, totalState, curState, 0);
		int numDivision = pow(CAST_CSIWDN(pro)->numberNode(), numSelect);
		int totalDivision = numDivision * totalState.size();
		divisions.clear();
		divisions.resize(totalDivision);
		for (int idSol(0); idSol < bestSolIds.size(); ++idSol) {
			int solId = bestSolIds[idSol];
			for (int stateId(0); stateId < totalState.size(); ++stateId) {
				int divisionId = getSolDivisionId(stateId, numDivision, CAST_CSIWDN(pro)->numberNode(),
					totalState[stateId], solVecInts[solId]);
				divisions[divisionId].push_back(solId);
			}
		}

		bestSolIds.resize(divisions.size());


		std::vector<std::thread> thrds;
		int num_task = std::thread::hardware_concurrency();
		int num_samples = divisions.size();
		std::vector<int> tasks;
		UTILITY::assignThreads(num_samples, num_task, tasks);
		std::pair<int, int> from_to;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			thrds.push_back(std::thread(
				calNBN_iterater_ThreadTask, from_to.first, from_to.second,
				pro, std::cref(divisions), std::ref(bestSolIds),
				std::ref(nbn_data)));
		}
		for (auto& thrd : thrds)
			thrd.join();

		unique_ids.clear();
		for (auto& it : bestSolIds) {
			if (it != -1)
				unique_ids.insert(it);
		}

		bestSolIds.clear();
		for (auto& it : unique_ids) {
			bestSolIds.push_back(it);
		}
	}

	std::cout << "calNBN_iterater final " << std::endl;
	calNBN_iterater(pro, bestSolIds, bestId, nbn_data);

}





void calNBN(ofec::Problem* pro, ofec::Random* rnd,
	NBNdataVecInt& data,
	int numSols = 2e4, int subSols = 100) {



	std::vector<std::vector<int>> subRegion(numSols);
	std::vector<int> bestSolIds(numSols);
	for (int idx(0); idx < numSols; ++idx) {
		int fromId = idx * subSols;
		for (int idy(0); idy < subSols; ++idy) {
			subRegion[idx].push_back(idx * subSols + idy);
		}
	}

	std::vector<std::thread> thrds;
	int num_task = std::thread::hardware_concurrency();
	int num_samples = subRegion.size();
	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);



	std::cout << "calNBN_iterater_ThreadTask" << std::endl;
	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];
		thrds.push_back(std::thread(
			calNBN_iterater_ThreadTask, from_to.first, from_to.second,
			pro, std::cref(subRegion), std::ref(bestSolIds),
			std::ref(data)));
	}
	for (auto& thrd : thrds)
		thrd.join();

	int bestId(0);


	std::cout << "calNBN_division" << std::endl;
	calNBN_division(pro, bestSolIds, bestId, data);

	data.m_bestSolId = bestId;

}





void generateCWIDNglobalStructure() {


	using namespace ofec;
	int numSols = 2e4, subSols = 1000;
	{
		std::string outputDir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data3/";

		std::filesystem::create_directory(outputDir);
		std::ofstream out(outputDir + "test.txt");
		out << "yes" << std::endl;
		std::cout << "output test file" << std::endl;
		out.close();
	}



	NBNdataVecInt data;


	//ParamMap params;
	//params["problem name"] = std::string("CSIWDN");
	//params["dataFile1"] = std::string("Net2");
	//params["dataFile2"] = std::string("case1");
	//params["dataFile3"] = std::string("case11");
	//params["use LSTM"] = false;
	//auto param = std::shared_ptr<const ParamMap>(new ParamMap(params));
	//auto pro = ofec::ADD_PRO(param, 0.1);
	//pro->initialize();
	//CAST_CSIWDN(pro.get())->setPhase(CAST_CSIWDN(pro.get())->totalPhase());

	std::shared_ptr<ofec::Problem> pro;
	generateDefaultPro(pro);

	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)> eval_fun =
		[](ofec::SolBase& sol, Problem* pro) {
		using namespace ofec;
		sol.evaluate(pro);
		ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
	};

	//Random rnd(0.5);

	std::shared_ptr<Random> rnd(new Random(0.5));
	std::cout << "generate solutions" << std::endl;
	generateSols(pro.get(), rnd.get(), eval_fun, data,
		numSols, subSols);

	int numSource(CAST_CSIWDN(pro.get())->numberSource());

	for (int idx(0); idx < data.m_sols.size(); ++idx) {
		ofec::CSIWDN::solutionType& cursol =
			dynamic_cast<ofec::CSIWDN::solutionType&>(*data.m_sols[idx]);
		data.m_sol_vecInt[idx].resize(numSource);
		for (int ids(0); ids < numSource; ++ids) {
			data.m_sol_vecInt[idx][ids] = cursol.variable().index(ids) - 1;
		}

	}



	std::cout << "calNBN" << std::endl;
	calNBN(pro.get(), rnd.get(), data,
		numSols, subSols);


	auto outputDir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data2/";

	//std::string dir = outputDir + std::string("/result/nbn_data/");
	//dir += "CWIDN_data_" + UTILITY::getCurrentSystemTime() + "/";
	//std::filesystem::create_directory(dir);
	std::cout << "output Data" << std::endl;

	{
		std::ofstream out(outputDir + std::string("CWIDN_bestsol") + ".txt");
		CAST_CSIWDN(pro.get())->outputSol(out, *data.m_sols[data.m_bestSolId]);
		out.close();
	}



	{
		std::ofstream out(outputDir + std::string("CWIDN_nbn") + ".txt");
		data.outputNBN(out);
		out.close();
	}

	{
		std::ofstream out(outputDir + std::string("CWIDN_nbn_sols") + ".txt");
		data.outputSols_CSWIDN(out, pro.get());
		out.close();
	}


	{
		std::vector<TreeGraphSimple::NBNdata> nbn_data;
		transferNBNData(nbn_data, data.m_solIds, data.m_belongs, data.m_dis2parent, data.m_fitness);
		TreeGraphSimple nbn_graph;
		nbn_graph.setNBNdata(nbn_data);
		nbn_graph.calNetwork();
		std::ofstream out(outputDir + std::string("CWIDN_nbn_networks") + ".txt");
		nbn_graph.outputNBNnetwork(out);
		out.close();
	}
}



void readData(ofec::Problem* pro,
	NBNdataVecInt& nbn_data,
	int numSols = 2e4, int subSols = 1000) {
	nbn_data.initialize(numSols * subSols);


}

void readDataCalNBN() {
	using namespace ofec;
	int numSols = 2e4, subSols = 1000;
	//	numSols = 2e4;
//	subSols = 1;

	std::string readfilename = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data1/CWIDN_nbn_sols_2023-10-27-07-01-50.txt";

	std::string outputDir = "";
	{
		//	std::string outputDir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data2/";

			//std::filesystem::create_directory(outputDir);
		std::ofstream out(outputDir + "test.txt");
		out << "yes" << std::endl;
		std::cout << "output test file" << std::endl;
		out.close();
	}



	NBNdataVecInt data;


	//ParamMap params;
	//params["problem name"] = std::string("CSIWDN");
	//params["dataFile1"] = std::string("Net2");
	//params["dataFile2"] = std::string("case1");
	//params["dataFile3"] = std::string("case11");
	//params["use LSTM"] = false;
	//auto param = std::shared_ptr<const ParamMap>(new ParamMap(params));
	//auto pro = ofec::ADD_PRO(param, 0.1);
	//pro->initialize();
	//CAST_CSIWDN(pro.get())->setPhase(CAST_CSIWDN(pro.get())->totalPhase());

	std::shared_ptr<ofec::Problem> pro;
	generateDefaultPro(pro);

	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)> eval_fun =
		[](ofec::SolBase& sol, Problem* pro) {
		using namespace ofec;
		sol.evaluate(pro);
		ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
	};

	//Random rnd(0.5);

	std::shared_ptr<Random> rnd(new Random(0.5));
	//std::cout << "generate solutions" << std::endl;
	//generateSols(pro.get(), rnd.get(), eval_fun, data,
	//	numSols, subSols);
		// 打开文件

	FILE* fp = fopen(readfilename.c_str(), "r");
	if (fp == NULL)
	{
		printf("open file failed!\n");
		return;
	}
	std::shared_ptr<ofec::SolBase> cursol;
	cursol.reset(CAST_CSIWDN(pro.get())->createSolution());
	data.initialize(numSols * subSols);
	auto& sols(data.m_sols);
	sols.clear();
	int solId = -1;

	ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;


	int numNode = CAST_CSIWDN(pro.get())->numberNode();
	int numSource = CAST_CSIWDN(pro.get())->numSource();
	while (CAST_CSIWDN(pro.get())->inputSol(fp, *cursol)) {
		sols.push_back(cursol);

		ofec::CSIWDN::solutionType& cursoltype =
			dynamic_cast<ofec::CSIWDN::solutionType&>(*cursol);
		//cursol->evaluate(pro.get());

		//double err = abs(cursol->fitness() - cursol ->objective()[0]*pos);
		//if (err > 1e-4) {
		//	std::cout << "error2 solution id" << solId << std::endl;
		//}


		++solId;

		cursol.reset(CAST_CSIWDN(pro.get())->createSolution());

		bool flag = true;
		//	ofec::CSIWDN::solutionType& cursoltype =
		//		dynamic_cast<ofec::CSIWDN::solutionType&>(*cursol);
		for (int ids(0); ids < numSource; ++ids) {
			if (cursoltype.variable().index(ids) < 1 || cursoltype.variable().index(ids) > numNode) {
				std::cout << "error3333 solution id" << solId << std::endl;
				std::cout << "ids " << ids << " cursoltype variable idx" << cursoltype.variable().index(ids) << std::endl;

				flag = false;
				//	break;
			}


		}

		if (!flag) {

			std::cout << "outputSol" << std::endl;
			CAST_CSIWDN(pro.get())->outputSol(std::cout, *cursol);
			return;
		}

		if (sols.size() >= numSols * subSols) break;

	}
	if (sols.size() != numSols * subSols) {
		std::cout << "error: input solutions size is wrong\t" << std::endl;
	}
	else std::cout << "correct solutions size" << std::endl;
	{
		int numNode = CAST_CSIWDN(pro.get())->numberNode();
		int numSource = CAST_CSIWDN(pro.get())->numSource();
		for (int idx(0); idx < sols.size(); ++idx) {
			ofec::CSIWDN::solutionType& cursoltype =
				dynamic_cast<ofec::CSIWDN::solutionType&>(*data.m_sols[idx]);
			for (int ids(0); ids < numSource; ++ids) {
				if (cursoltype.variable().index(ids) < 1 || cursoltype.variable().index(ids) > numNode) {
					std::cout << "error1 solution id" << idx << std::endl;
					break;
				}
			}

		}


		int testSol = 1e2;
		ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
		while (testSol--) {
			int randId = rnd->uniform.nextNonStd<int>(0, sols.size());
			sols[randId]->evaluate(pro.get());

			double err = abs(sols[randId]->fitness() - sols[randId]->objective()[0]);

			if (err > 1e-4) {
				std::cout << "error2 solution id" << randId << std::endl;
			}
		}
	}


	//int numSource(CAST_CSIWDN(pro.get())->numberSource());

	std::cout << "transfer network" << std::endl;
	for (int idx(0); idx < data.m_sols.size(); ++idx) {
		ofec::CSIWDN::solutionType& cursol =
			dynamic_cast<ofec::CSIWDN::solutionType&>(*data.m_sols[idx]);
		data.m_sol_vecInt[idx].resize(numSource);
		for (int ids(0); ids < numSource; ++ids) {
			data.m_sol_vecInt[idx][ids] = cursol.variable().index(ids) - 1;
		}

	}



	std::cout << "calNBN" << std::endl;
	calNBN(pro.get(), rnd.get(), data,
		numSols, subSols);




	//std::string dir = outputDir + std::string("/result/nbn_data/");
	//dir += "CWIDN_data_" + UTILITY::getCurrentSystemTime() + "/";
	//std::filesystem::create_directory(dir);
	std::cout << "output Data" << std::endl;

	{
		std::ofstream out(outputDir + std::string("CWIDN_bestsol") + ".txt");
		CAST_CSIWDN(pro.get())->outputSol(out, *data.m_sols[data.m_bestSolId]);
		out.close();
	}



	{
		std::ofstream out(outputDir + std::string("CWIDN_nbn") + ".txt");
		data.outputNBN(out);
		out.close();
	}

	{
		std::ofstream out(outputDir + std::string("CWIDN_nbn_sols") + ".txt");
		data.outputSols_CSWIDN(out, pro.get());
		out.close();
	}


	{
		std::vector<TreeGraphSimple::NBNdata> nbn_data;
		transferNBNData(nbn_data, data.m_solIds, data.m_belongs, data.m_dis2parent, data.m_fitness);
		TreeGraphSimple nbn_graph;
		nbn_graph.setNBNdata(nbn_data);
		nbn_graph.calNetwork();
		std::ofstream out(outputDir + std::string("CWIDN_nbn_networks") + ".txt");
		nbn_graph.outputNBNnetwork(out);
		out.close();
	}



}

void testSolData(ofec::Problem* pro, ofec::Random* rnd,
	NBNdataVecInt& data,
	int numSols = 2e4, int subSols = 1000) {
	using namespace ofec;

	auto& sols = data.m_sols;
	if (sols.size() != numSols * subSols) {
		std::cout << "error: input solutions size is wrong\t" << std::endl;
	}
	else std::cout << "correct solutions size" << std::endl;
	{
		int numNode = CAST_CSIWDN(pro)->numberNode();
		int numSource = CAST_CSIWDN(pro)->numSource();
		for (int idx(0); idx < sols.size(); ++idx) {
			ofec::CSIWDN::solutionType& cursoltype =
				dynamic_cast<ofec::CSIWDN::solutionType&>(*data.m_sols[idx]);
			for (int ids(0); ids < numSource; ++ids) {
				if (cursoltype.variable().index(ids) < 1 || cursoltype.variable().index(ids) > numNode) {
					std::cout << "error1 solution id" << idx << std::endl;
					break;
				}
			}

		}


		int testSol = 1e2;
		ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
		while (testSol--) {
			int randId = rnd->uniform.nextNonStd<int>(0, sols.size());
			sols[randId]->evaluate(pro);

			double err = abs(sols[randId]->fitness() - pos * sols[randId]->objective()[0]);

			if (err > 1e-4) {
				std::cout << "error2 solution id" << randId << std::endl;
			}
		}
	}

}

void runNBN(
	const std::string& outputDir,
	ofec::Problem* pro, ofec::Random* rnd,
	NBNdataVecInt& data,
	int numSols = 2e4, int subSols = 1000
) {
	using namespace ofec;


	int numSource = CAST_CSIWDN(pro)->numSource();

	std::cout << "transfer network" << std::endl;
	for (int idx(0); idx < data.m_sols.size(); ++idx) {
		ofec::CSIWDN::solutionType& cursol =
			dynamic_cast<ofec::CSIWDN::solutionType&>(*data.m_sols[idx]);
		data.m_sol_vecInt[idx].resize(numSource);
		for (int ids(0); ids < numSource; ++ids) {
			data.m_sol_vecInt[idx][ids] = cursol.variable().index(ids) - 1;
		}

	}



	std::cout << "calNBN" << std::endl;
	calNBN(pro, rnd, data,
		numSols, subSols);




	//std::string dir = outputDir + std::string("/result/nbn_data/");
	//dir += "CWIDN_data_" + UTILITY::getCurrentSystemTime() + "/";
	//std::filesystem::create_directory(dir);
	std::cout << "output Data" << std::endl;

	{
		std::ofstream out(outputDir + std::string("CWIDN_bestsol") + ".txt");
		CAST_CSIWDN(pro)->outputSol(out, *data.m_sols[data.m_bestSolId]);
		out.close();
	}



	{
		std::ofstream out(outputDir + std::string("CWIDN_nbn") + ".txt");
		data.outputNBN(out);
		out.close();
	}




	{
		std::vector<TreeGraphSimple::NBNdata> nbn_data;
		transferNBNData(nbn_data, data.m_solIds, data.m_belongs, data.m_dis2parent, data.m_fitness);
		TreeGraphSimple nbn_graph;
		nbn_graph.setNBNdata(nbn_data);
		nbn_graph.calNetwork();
		std::ofstream out(outputDir + std::string("CWIDN_nbn_networks") + ".txt");
		nbn_graph.outputNBNnetwork(out);
		out.close();
	}

	{
		std::ofstream out(outputDir + std::string("CWIDN_nbn_sols") + ".txt");
		data.outputSols_CSWIDN(out, pro);
		out.close();
	}



}

std::vector<std::shared_ptr<ofec::SolBase>> genSols;
std::vector<std::shared_ptr<ofec::SolBase>> inputSols;
int solId = -1;




void readSolBuffer(
	const std::string& filepath,
	ofec::Problem* pro,
	std::vector<std::shared_ptr<ofec::SolBase>>& sols) {

	std::stringstream in;
	ofec::readFileBuffer(filepath, in);



	std::ofstream out("test.txt");
	out << in.str();
	out.close();

	int stop = -1;
	using namespace ofec;
	std::shared_ptr<ofec::SolBase> curSol(pro->createSolution());
	ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	sols.clear();
	//std::vector<std::shared_ptr<ofec::SolBase>> sols;
	int solId(0);
	//while (in) 
	while (!in.eof())
	{

		in >> solId;
		if (!in) break;
		curSol.reset(pro->createSolution());
		::solId++;
		CAST_CSIWDN(pro)->inputSol(in, *curSol);

		curSol->objective()[0] = curSol->fitness() * pos;
		sols.push_back(curSol);
	}

}


void testInputSol(ofec::Problem* pro, std::ifstream& in,
	std::vector<std::shared_ptr<ofec::SolBase>>& sols) {

	using namespace ofec;
	std::shared_ptr<ofec::SolBase> curSol(pro->createSolution());
	ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	sols.clear();
	//std::vector<std::shared_ptr<ofec::SolBase>> sols;

	int solId(0);




	while (in) {

		in >> solId;
		if (!in) break;
		curSol.reset(pro->createSolution());
		::solId++;
		CAST_CSIWDN(pro)->inputSol(in, *curSol);
		curSol->objective()[0] = curSol->fitness() * pos;

		//curSol->evaluate(pro);
	//	curSolX = dynamic_cast<ofec::CSIWDN::solutionType&>(*curSol);
	//	saveSolX= 
		//if (abs(pos * curSol->objective()[0] - curSol->fitness())>0.001) {
		//	std::cout << "error info\t" << std::endl;
		//}
		sols.push_back(curSol);
		//sols.push_back(curSol);
	}


}

void testOutputSol(ofec::Problem* pro, ofec::Random* rnd, std::ofstream& out,
	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)>& eval_fun,
	std::vector<std::shared_ptr<ofec::SolBase>>& sols
) {
	using namespace ofec;
	std::shared_ptr<ofec::SolBase> cursol;
	cursol.reset(pro->createSolution());
	cursol->initialize(pro, rnd);
	eval_fun(*cursol, pro);
	out << "1" << "\t";
	CAST_CSIWDN(pro)->outputSol(out, *cursol);

	sols.push_back(cursol);
}






void testFile() {

	auto outputDir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data1/";
	outputDir = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data_paper_com/CWIDN_nbn_data/data1/CWIDN_nbn_sols_2023-10-27-07-01-50.txt";
	//	const std::string filepath= outputDir:
	//	outputDir = "cwidn_test.txt";
	//	std::ifstream in(outputDir);

	std::shared_ptr<ofec::Problem> pro;
	generateDefaultPro(pro);
	std::shared_ptr<ofec::Random> rnd(new ofec::Random(0.5));

	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)> eval_fun =
		[](ofec::SolBase& sol, ofec::Problem* pro) {
		using namespace ofec;
		sol.evaluate(pro);
		ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
	};
	{
		readSolBuffer(outputDir, pro.get(), inputSols);
	}
	std::cout << "inputSolSize\t" << inputSols.size() << std::endl;
	//while (true) {
	//	{
	//		std::ofstream out(outputDir);
	//		int numSol(10);
	//		while (numSol--) {
	//			testOutputSol(pro.get(), rnd.get(), out, eval_fun, genSols);
	//		}
	//		out.close();
	//	}
	//	{
	//		std::ifstream in(outputDir);
	//		inputSols.clear();
	//		while (in) {
	//			testInputSol(pro.get(), in, inputSols);
	//		}
	//	}
	//}
}




void readSols() {

}

bool judgeSolTest(ofec::SolBase& sol, ofec::Problem* pro) {
	ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	sol.evaluate(pro);
	double err = abs(sol.fitness() - pos * sol.objective()[0]);

	if (err > 1e-4) {
		std::cout << "error2 solution id" << std::endl;
		return false;

	}
	else return true;
}


bool judgeSolCorrect(ofec::SolBase& sol, ofec::Problem* pro) {
	using namespace ofec;
	int numNode = CAST_CSIWDN(pro)->numberNode();
	int numSource = CAST_CSIWDN(pro)->numSource();
	{
		ofec::CSIWDN::solutionType& cursoltype =
			dynamic_cast<ofec::CSIWDN::solutionType&>(sol);
		for (int ids(0); ids < numSource; ++ids) {
			if (cursoltype.variable().index(ids) < 1 || cursoltype.variable().index(ids) > numNode) {


				return false;
				break;
			}
		}

	}
	return true;
}

void readCWIDNsolThreadTask(
	ofec::Problem* pro,
	const std::vector<std::unique_ptr<char[]>>& totalData,
	int from, int to,
	bool firstFlag, NBNdataVecInt& nbn_data
) {

	using namespace ofec;
	std::string inputString(totalData[from].get());

	std::string tmpStr;
	std::stringstream tmpIn;
	size_t fromIdx = 0;
	size_t beginId(0), endId(-1);

	int totalLine = 4;
	int lastLine = 4;
	int solId = 0;

	int numSource = CAST_CSIWDN(pro)->numSource();

	if (!firstFlag) {

		while (solId < numSource) {
			fromIdx = endId;

			beginId = fromIdx;
			endId = inputString.find_first_of('\n', beginId + 1);
			tmpStr = inputString.substr(beginId, endId);
			tmpIn << tmpStr;
			//if(tmpIn.)
				//tmpStr.clear();
				//fromIdx = endId + 1;
			tmpIn >> solId;
		}

		tmpIn.clear();

		//fromIdx = beginId;
	}
	else {
		beginId = fromIdx;
		endId = inputString.find_first_of('\n', beginId + 1);
		fromIdx = endId;
		tmpStr = inputString.substr(beginId, endId);
		tmpIn << tmpStr;
		--lastLine;
	}

	std::vector<std::string> tmpInfos;
	while ((from < to) || (from == to && totalLine != lastLine)) {
		{
			beginId = fromIdx;
			endId = inputString.find_first_of('\n', beginId + 1);


			if (endId < inputString.size()) {
				fromIdx = endId;
				tmpStr = inputString.substr(beginId, endId);
				tmpIn << tmpStr;
				tmpInfos.push_back(tmpStr);
				--lastLine;
			}
			else {
				tmpStr = inputString.substr(beginId, inputString.size());
				tmpIn << tmpStr;
				fromIdx = inputString.size();


			}

			if (lastLine == 0) {
				tmpIn >> solId;
				CAST_CSIWDN(pro)->inputSol(tmpIn, *nbn_data.m_sols[solId]);
				if (!judgeSolTest(*nbn_data.m_sols[solId], pro)) {
					std::cout << "err cur sol id " << solId << std::endl;
				}
				lastLine = totalLine;
			}

			if (fromIdx == inputString.size()) {
				fromIdx = 0;
				++from;
				if (from < to) {
					inputString = totalData[from].get();
				}
				else if (from == to && totalLine != lastLine) {
					inputString = totalData[from].get();
				}
			}

		}
	}

	//std::string 


}


size_t findSolIdx(int from, ofec::Problem* pro,
	const std::vector<std::unique_ptr<char[]>>& totalData) {
	using namespace ofec;
	using namespace std;
	std::string inputString(totalData[from].get());
	size_t fromIdx = 0;
	size_t beginId(0), endId(0);
	int numSource = CAST_CSIWDN(pro)->numberNode();
	int totalLine = 4;
	int lastLine = 4;
	int solId = 0;

	std::string tmpStr;
	while (solId <= numSource) {
		fromIdx = endId;
		beginId = fromIdx;
		beginId = inputString.find_first_of(10, fromIdx);
		endId = inputString.find_first_of(10, beginId + 1);
		//	std::cout << "fromId\t" << fromIdx << "\tbeginId\t" << beginId << "\t" << endId << std::endl;
		tmpStr = inputString.substr(beginId, endId - beginId);
		std::stringstream tmpIn;
		//	std::cout << tmpStr << std::endl;
		tmpIn << tmpStr << std::endl;
		//if(tmpIn.)
			//tmpStr.clear();
			//fromIdx = endId + 1;
		tmpIn >> solId;
		tmpIn.clear();
		//	std::cout << "cur Sol Id\t" << solId << std::endl;
	}
	//std::cout << "fromId\t" << fromIdx << "\tbeginId\t" << beginId << "\t" << endId << std::endl;
	//std::cout << "substring idx" << from << std::endl;
	//std::cout << tmpStr << std::endl;
	return beginId;

}


void splitString(
	ofec::Problem* pro,
	const std::vector<std::unique_ptr<char[]>>& totalData,
	std::vector<std::stringstream>& info) {

	info.resize(totalData.size());
	size_t beginId = 0;
	std::string tmp;
	int idx(0);
	for (; idx + 1 < totalData.size(); ++idx) {
		tmp = totalData[idx].get();
		info[idx] << tmp.substr(beginId, tmp.size());
		beginId = findSolIdx(idx + 1, pro, totalData);
		if (beginId) {
			tmp = totalData[idx + 1].get();
			info[idx] << tmp.substr(0, beginId);
		}
	}
	//	++idx;
	tmp = totalData[idx].get();
	info[idx] << tmp.substr(beginId, tmp.size());
}



void readDataThreadTask(int from, int to,
	ofec::Problem* pro,
	std::vector<std::stringstream>& splitData,
	NBNdataVecInt& nbn_data) {
	using namespace ofec;
	bool firstFlag2 = true;
	int dataId = -1;
	int solIter = 0;
	for (int idx(from); idx < to; ++idx) {
		auto& tmpIn = splitData[idx];
		dataId++;
		solIter = -1;
		while (!tmpIn.eof()) {
			++solIter;
			tmpIn >> solId;
			if (tmpIn.eof())break;
			if (firstFlag2) {
				firstFlag2 = false;
			}
			else if (solId < 0) {
				break;
			}
			auto& cursol = nbn_data.m_sols[solId];
			CAST_CSIWDN(pro)->inputSol(tmpIn, *cursol);
		}
	}
}



void readSolData(ofec::Problem* pro,
	const std::string& readfilename, NBNdataVecInt& nbn_data) {
	using namespace std;
	using namespace ofec;

	ifstream bigFile(readfilename);
	constexpr size_t bufferSize = 2 * 1024 * 1024;
	unique_ptr<char[]> buffer(new char[bufferSize]);
	std::vector<std::stringstream> splitData;
	{
		std::vector<unique_ptr<char[]>> totalData;
		std::stringstream tmpIn;
		std::string tmp;
		int line = 0;
		while (bigFile)
		{
			bigFile.read(buffer.get(), bufferSize);
			//	std::cout << "read file line" << line++ << std::endl;
			totalData.emplace_back(buffer.release());
			buffer.reset(new char[bufferSize]);

		}
		std::cout << "total buffer " << totalData.size() << std::endl;
		splitString(pro, totalData, splitData);
		//	std::cout << totalData.back() << std::endl;
	}
	{
		int num_task = std::thread::hardware_concurrency();
		int num_samples = splitData.size();
		std::vector<std::thread> thrds;
		std::vector<int> tasks;
		UTILITY::assignThreads(num_samples, num_task, tasks);
		std::pair<int, int> from_to;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			thrds.push_back(std::thread(readDataThreadTask,
				tasks[i], tasks[i + 1], pro, std::ref(splitData), std::ref(nbn_data)));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}

	return;

	//bool firstFlag2 = true;
	//int dataId = -1;
	//int solIter = 0;
	//for (auto& tmpIn : splitData) {
	//	dataId++;
	//	solIter = -1;
	//	while (!tmpIn.eof()) {
	//		++solIter;
	//		tmpIn >> solId;
	//		//if (!) break;
	//		if (tmpIn.eof())break;
	//		if (firstFlag2) {
	//			firstFlag2 = false;
	//		}
	//		else if (solId < 0) {
	//			break;
	//		}
	//		auto& cursol = nbn_data.m_sols[solId];
	//		CAST_CSIWDN(pro)->inputSol(tmpIn, *cursol);
	//	}
	//}
}



void testReadFile() {
	using namespace std;

	using namespace ofec;

	int numSols = 2e4;
	int subSols = 1000;
	std::string readfilename = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data1/CWIDN_nbn_sols_2023-10-27-07-01-50.txt";
	//readfilename = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data1/CWIDN_nbn_sols_2023-10-26-19-38-21.txt";
	ifstream bigFile(readfilename);
	constexpr size_t bufferSize = 2 * 1024 * 1024;
	unique_ptr<char[]> buffer(new char[bufferSize]);
	std::shared_ptr<ofec::Problem> pro;
	std::shared_ptr<ofec::Random> rnd(new Random(0.5));
	generateDefaultPro(pro);

	std::vector<std::stringstream> splitData;
	{
		std::vector<unique_ptr<char[]>> totalData;
		std::stringstream tmpIn;
		std::string tmp;
		int line = 0;
		while (bigFile)
		{
			bigFile.read(buffer.get(), bufferSize);
			//	std::cout << "read file line" << line++ << std::endl;
			totalData.emplace_back(buffer.release());
			buffer.reset(new char[bufferSize]);
			// process data in buffer
		}

		std::cout << "total buffer " << totalData.size() << std::endl;
		splitString(pro.get(), totalData, splitData);

		//	std::cout << totalData.back() << std::endl;
	}




	NBNdataVecInt nbn_data;



	nbn_data.initialize(numSols * subSols);
	nbn_data.assignSolMem(pro.get());

	int from = 0;
	bool firstFlag = true;
	int to = 3;

	int solId(0);
	//for (int idx(0); idx < numSols * subSols;++idx) {
	//	tmpIn >> solId;
	//	auto& cursol = nbn_data.m_sols[idx];
	//	CAST_CSIWDN(pro.get())->inputSol(tmpIn, *cursol);
	////	if (!judgeSolTest(*nbn_data.m_sols[idx], pro.get())) {
	//////		std::cout << "err cur sol id " << solId << std::endl;
	////	}
	//	//int stop = -1;
	//}

	bool firstFlag2 = true;
	int dataId = -1;
	int solIter = 0;
	for (auto& tmpIn : splitData) {
		dataId++;
		solIter = -1;
		while (!tmpIn.eof()) {
			++solIter;
			tmpIn >> solId;
			//if (!) break;
			if (tmpIn.eof())break;
			if (firstFlag2) {
				firstFlag2 = false;
			}
			else if (solId < 0) {
				break;
			}
			auto& cursol = nbn_data.m_sols[solId];
			CAST_CSIWDN(pro.get())->inputSol(tmpIn, *cursol);
			//if (!judgeSolCorrect(*cursol, pro.get())) {
			//	std::cout << "dataIder\t" << dataId << "\t" << "solIter\t" << solIter << std::endl;
			//	std::cout << "solId\t" << solId << std::endl;
			//	CAST_CSIWDN(pro.get())->outputSol(std::cout, *cursol);
			//	std::cout << endl;
			//}

		//	solId = -1;
		}
	}

	testSolData(pro.get(), rnd.get(), nbn_data, numSols, subSols);
	runNBN("", pro.get(), rnd.get(), nbn_data,
		numSols, subSols);
}


void outputSolFit(const std::string& outputFile,
	ofec::Problem* pro,
	const std::vector<std::shared_ptr<ofec::SolBase>>& sols,
	int fitness_idx) {
	using namespace ofec;


	std::cout << "output idx " << fitness_idx << std::endl;
	{
		{
			std::ofstream out(outputFile + std::string("CWIDN_sol") + std::to_string(fitness_idx) + ".txt");
			CAST_CSIWDN(pro)->outputSol(out, *sols.front());
			out.close();
		}


		std::vector<double> fitness(sols.size());
		for (int idx(0); idx < fitness.size(); ++idx) {
			fitness[idx] = sols[idx]->fitness();
		}

		{
			std::ofstream out(outputFile + std::string("CWIDN_fitness") + std::to_string(fitness_idx) + ".txt");
			for (int idx(0); idx < fitness.size(); ++idx) {
				out << idx << "\t" << fitness[idx] << std::endl;
			}
			out.close();
		}

	}
}


void readDataRunNBN(const std::string readfilename,
	const std::string& outputFile,
	int numSols = 2e4, int subSols = 1000) {
	using namespace ofec;


	//	std::string readfilename = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data1/CWIDN_nbn_sols_2023-10-27-07-01-50.txt";
		//readfilename = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data1/CWIDN_nbn_sols_2023-10-26-19-38-21.txt";
	std::shared_ptr<ofec::Problem> pro;
	std::shared_ptr<ofec::Random> rnd(new Random(0.5));
	generateDefaultPro(pro);

	NBNdataVecInt nbn_data;
	nbn_data.initialize(numSols * subSols);
	nbn_data.assignSolMem(pro.get());
	readSolData(pro.get(), readfilename, nbn_data);

	testSolData(pro.get(), rnd.get(), nbn_data, numSols, subSols);



	// generate global structure
	{
		//runNBN(outputFile, pro.get(), rnd.get(), nbn_data,
		//	numSols, subSols);
	}

	// generte local structure

	std::shared_ptr<ofec::SolBase> bestSol = nbn_data.m_sols.front();
	for (auto& it : nbn_data.m_sols) {
		if (bestSol->fitness() < it->fitness()) {
			bestSol = it;
		}
	}
	{
		std::ofstream out(outputFile + std::string("CWIDN_bestsol") + ".txt");
		CAST_CSIWDN(pro.get())->outputSol(out, *bestSol);
		out.close();
	}

	int sampleSize = 1e6;
	// init continous solutions

	//NBNdataVecInt nbn_data2;
	//nbn_data2.initialize(sampleSize);

	std::vector<std::shared_ptr<ofec::SolBase>>& samples(nbn_data.m_sols);

	for (auto& it : samples) {
		it.reset(pro->createSolution(*bestSol));
		CAST_CSIWDN(pro.get())->initSolutionMultiplier(*it, rnd.get());
	}

	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)> eval_fun =
		[](ofec::SolBase& sol, ofec::Problem* pro) {
		using namespace ofec;
		sol.evaluate(pro);
		ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
	};



	evaluateSolMultiTask(pro.get(), eval_fun, samples);
	nbn_data.updateFitness();
	{
		std::ofstream out(outputFile + "CWIDN_continous_solutions.txt");
		nbn_data.outputSols_CSWIDN(out, pro.get());
		out.close();
	}
	int fitness_idx = 0;
	outputSolFit(outputFile, pro.get(), samples, fitness_idx);

	while (true) {
		++fitness_idx;
		int rndSource = rnd->uniform.nextNonStd<int>(1, CAST_CSIWDN(pro.get())->numSource() + 1);

		CAST_CSIWDN(pro.get())->initSolutionPosition(*samples.front(), rnd.get(), rndSource);

		for (int idx(1); idx < samples.size(); ++idx) {
			CAST_CSIWDN(pro.get())->setSolutionPosition(*samples[idx], *samples.front());
		}


		evaluateSolMultiTask(pro.get(), eval_fun, samples);
		nbn_data.updateFitness();

		outputSolFit(outputFile, pro.get(), samples, fitness_idx);


	}

}



namespace ofec {

	void registerParamAbbr() {};

	void customizeFileName() {};

	void run() {
		registerParamAbbr();
		customizeFileName();
		using namespace std;
		//	testReadFile();
			//readDataCalNBN();
			//generateCWIDNglobalStructure();


		std::string readfilename = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data1/CWIDN_nbn_sols_2023-10-27-07-01-50.txt";
		//readfilename = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data1/CWIDN_nbn_sols_2023-10-26-19-38-21.txt";
		std::string outputfilepath = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/data_continous_data/";
		//readfilename = "/home/lab408/share/2018/diaoyiya/paper_com_experiment_data/cwidn/read_data/CWIDN_nbn_sols_2023-10-27-07-01-50.txt";
		//outputfilepath = "/home/lab408/share/2018/diaoyiya/paper_com_experiment_data/cwidn/cwidn_continous/";

		std::filesystem::create_directory(outputfilepath);

		readDataRunNBN(readfilename, outputfilepath);


	}


}
