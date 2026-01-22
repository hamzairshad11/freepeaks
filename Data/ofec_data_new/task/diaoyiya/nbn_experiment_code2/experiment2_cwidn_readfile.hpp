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





void test() {


	using namespace std;
	using namespace ofec;


	ParamMap params;
	params["problem name"] = std::string("CSIWDN");
	params["dataFile1"] = std::string("Net2");
	params["dataFile2"] = std::string("case1");
	params["dataFile3"] = std::string("case11");
	params["use LSTM"] = false;





	auto eval_fun =
		[](ofec::SolBase& sol, Problem* pro) {
		using namespace ofec;
		sol.evaluate(pro);
		ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
	};





	//SECTION("CASE EPANET") {
	//	int id_param = ADD_PARAM(params);
	//	Problem *pro = ADD_PRO(id_param, 0.1);
	//	pro->initialize();
	//}

	//SECTION("CASE LSTM") {
	//	int id_param = ADD_PARAM(params);
	//	Problem *pro = ADD_PRO(id_param, 0.1);
	//	pro->initialize();
	//	GET_CSIWDN(pro)->setUseLSTM(true);

	//	Solution<VarCSIWDN> sol(1, 0, GET_CSIWDN(pro)->numberSource());

	//	Random &rnd = ADD_RND(0.5);
	//	sol.initialize(pro, rnd);
	//	sol.evaluate(pro, -1);

	//	GET_CSIWDN(pro)->setAlgorithmStart();
	//	sol.initialize(pro, rnd);
	//	for (size_t i = 0; i < 10000; ++i)
	//		sol.evaluate(pro, -1);

	//	DEL_PRO(pro);
	//}


	auto id_param = std::shared_ptr<const ParamMap>(new ParamMap(params));
	auto pro = ofec::ADD_PRO(id_param, 0.1);
	pro->initialize();
	//Random rnd(0.5);
	CAST_CSIWDN(pro.get())->setPhase(CAST_CSIWDN(pro.get())->totalPhase());
	std::shared_ptr<Random> rnd(new Random(0.5));
	//	pro->initialize();


	size_t num_sols = 100;

	std::vector<std::shared_ptr<ofec::SolBase>> sols(num_sols);
	for (size_t i = 0; i < num_sols; ++i) {
		sols[i].reset(pro->createSolution());
		sols[i]->initialize(pro.get(), rnd.get());
	}

	std::vector<Real> objs(num_sols);
	std::vector<std::vector<bool>> epa_domin(num_sols, std::vector<bool>(num_sols));
	for (size_t i = 0; i < num_sols; ++i)
		sols[i]->evaluate(pro.get());



	//std::vector<std::vector<bool>> epa_domin(num_sols, std::vector<bool>(num_sols));
	//for (size_t i = 0; i < num_sols; ++i)
	//	sols[i].evaluate(pro, -1);
	//for (size_t i = 0; i < num_sols; i++) {
	//	for (size_t j = i + 1; j < num_sols; j++) {
	//		epa_domin[i][j] = sols[i].dominate(sols[j], pro);
	//	}
	//}
	for (size_t i = 0; i < num_sols; ++i)
		objs[i] = sols[i]->objective(0);


}


#include <filesystem>
#include <string>
#include <vector>
#include <iostream>

#include "../utility/nbn_visualization/nbn_grid_tree_division_allSize.h"


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

		////auto param = ADD_PARAM(pro->getParam());
		//cur_pro = ofec::ADD_PRO(pro->getParam(), rnd->uniform.next());
		//cur_pro->initialize();
		////Random rnd(0.5);
		//CAST_CSIWDN(cur_pro.get())->setPhase(CAST_CSIWDN(cur_pro.get())->totalPhase());
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

		//for (endId = beginId+1; endId < inputString.size(); ++endId) {
		//	if (inputString[endId] == '\n') {
		//		break;
		//	}
		//}


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



void readDataThreadTask(int from, int to, ofec::Problem* pro,
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
			//if (!) break;
			if (tmpIn.eof())break;
			if (firstFlag2) {
				firstFlag2 = false;
			}
			else if (solId < 0) {
				break;
			}
			auto& cursol = nbn_data.m_sols[solId];
			CAST_CSIWDN(pro)->inputSol(tmpIn, *cursol);
			//if (!judgeSolCorrect(*cursol, pro.get())) {
			//	std::cout << "dataIder\t" << dataId << "\t" << "solIter\t" << solIter << std::endl;
			//	std::cout << "solId\t" << solId << std::endl;
			//	CAST_CSIWDN(pro.get())->outputSol(std::cout, *cursol);
			//	std::cout << endl;
			//}

		//	solId = -1;
		}
	}
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

			//	tmp = buffer.get();
			//	tmpIn << tmp;

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





namespace ofec {

	void registerParamAbbr();
	void customizeFileName();

	void run() {
		registerParamAbbr();
		customizeFileName();
		using namespace std;
		testReadFile();
		//readDataCalNBN();
		//generateCWIDNglobalStructure();



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
