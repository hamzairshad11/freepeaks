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
//#include "../utility/local_optima_network/lon.h"
#include <filesystem>
//#include "../utility/local_optima_network/lon.h"
//#include "../utility/fitness_landscape_analysis/walk/random_walk.h"
#include "../utility/nbn_visualization/nbn_onemax_simple.h"
#include "../utility/function/custom_function.h"

#include "../utility/functional.h"

#include "../utility/nbn_visualization/nbn_fla.h"

//#include "../instance/problem/continuous/single_objective/global/classical/weierstrass.h"


#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../instance/problem/continuous/free_peaks/param_creator.h"
#include "../instance/problem/continuous/free_peaks/parameter_function.h"

//#include "../instance/algorithm/visualize/LON/sampling/lon_con_sampling.h"
//#include "../instance/algorithm/visualize/LON/experiment/lon.h"
#include "../core/problem/continuous/continuous.h"

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







void calNBN_Analays(ofec::ParamMap& params,
	const std::string& filepath,
	std::ofstream& out, std::mutex& mtx
) {
	using namespace ofec;
	{


		auto id_param = std::shared_ptr<const ParamMap>(new ParamMap(params));
		auto pro = ofec::ADD_PRO(id_param, 0.1);
		pro->initialize();
		//Random rnd(0.5);

		std::shared_ptr<Random> rnd(new Random(0.5));
		int sampleSize = params.get<int>("sample size", 50);

		int numSample_int = sampleSize;
		{
			int dim = pro->numVariables();
			std::vector<int> m_dim_div(dim);
			int dim_div = std::max(2.0, exp(log(double(numSample_int)) / double(dim)));
			numSample_int = 1;
			for (auto& it : m_dim_div) {
				it = dim_div;
				numSample_int *= it;
			}
		}

		int rank = sampleSize / numSample_int;
		numSample_int *= rank;
		//if (numSample_int < sampleSize) {
		//	numSample_int = numSample_int * 2;
		//}

		if (sampleSize <= 10000) numSample_int = sampleSize;

		//	numSample_int = sampleSize;

		std::string proname = params.get<std::string>("problem name");
		if (proname == "free_peaks") {

			//par["taskName"] = curcurInfo.m_taskname;
			proname = params.get<std::string>("taskName");

		}

		//	int sampleSize = numDiv * numDiv;


		auto range = CAST_CONOP(pro.get())->boundary();


		{


			ofec::NBN_DivisionData nbn_network;
			nbn_network.m_division.reset(new ofec::NBN_GridTreeDivision);
			//ofec::NBN_VisualizationData nbn_data;

			auto& grid_division = dynamic_cast<ofec::NBN_GridTreeDivision&>(*nbn_network.m_division);


			auto eval_fun =
				[](ofec::SolBase& sol, Problem* pro) {
				using namespace ofec;
				sol.evaluate(pro);
				ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
				sol.setFitness(pos * sol.objective(0));
			};

			{
				auto& division(nbn_network.m_division);
				division->setMaxSampleSize(numSample_int);
				division->setMultiThread(true);
				division->setEvaluateMultiTask(true);

				int numThread = std::thread::hardware_concurrency();
				division->setNumberThread(numThread);
				//std::cout << "number of thread" << numThread << std::endl;
				auto start = std::chrono::system_clock::now();
				division->initialize(pro, rnd, eval_fun, false);
				auto end = std::chrono::system_clock::now();



				std::chrono::duration<float> difference = end - start;
				double milliseconds = difference.count();



				//std::cout << "read bufer time (s) \t" << milliseconds << std::endl;

				std::vector<std::shared_ptr<ofec::SolBase>> sols;
				std::vector<double> fitness;
				std::vector<int> belong;
				std::vector<double> dis2parent;
				std::vector<bool> flagOpt;
				std::vector<ofec::TreeGraphSimple::SolInfo> info;
				std::vector<int> solIds(sols.size());

				division->getSharedNBN(sols, fitness, belong, dis2parent, flagOpt);
				int originSize = sols.size();
				//	curInfo.m_curSampleSize = sols.size();
				//	globalInfo.outputAlgInfo(filepath, milliseconds);
				solIds.resize(sols.size());
				for (int idx(0); idx < solIds.size(); ++idx) {
					solIds[idx] = idx + 1;
				}

				int dim = params.get<int>("number of variables");
				ofec::TreeGraphSimple::SolInfo::transfer(info, sols, flagOpt);

				std::string curfilepath = filepath + "Proname_" + proname + "_dim_" + std::to_string(dim) + "_numSamples_" + std::to_string(sols.size());








				{
					std::ofstream nbnOut(curfilepath + "_nbn.txt");
					ofec::outputNBNdata(nbnOut, solIds, belong, dis2parent, fitness);
					nbnOut.close();




					//std::ofstream infoOut(curfilepath + "_nbnSolInfo.txt");
					//ofec::TreeGraphSimple::SolInfo::outputSolInfo(infoOut, info);
					//infoOut.close();


				}



				sols.clear();

				nbn_network.m_division = nullptr;


			}
		}


	}
}








void setParByThreadInfo(ofec::ParamMap& par, const ThreadInfoNBN& curcurInfo) {
	par["problem name"] = curcurInfo.m_proname;

	//params["problem name"] = std::string("Classic_Weierstrass");
	par["number of variables"] = curcurInfo.m_dim;;
	par["generation type"] = std::string("read_file");
	//3_s5_1_s2.txt;
	par["dataFile1"] = curcurInfo.m_filename;
	par["taskName"] = curcurInfo.m_taskname;
	par["sample size"] = curcurInfo.m_sample_size;
}


enum class TaskType { NBN_analyses = 0, LonContinous_analyses };

TaskType m_task_type;


void calTask(std::unique_ptr<utility::GeneralMultiThreadInfo>& curInfo,
	std::unique_ptr<utility::GeneralMultiThreadInfo>& globalInfo) {
	auto& curGlobalInfo = dynamic_cast<ThreadGlobalNBN&>(*globalInfo);
	auto& curcurInfo = dynamic_cast<ThreadInfoNBN&>(*curInfo);

	ofec::ParamMap par;
	setParByThreadInfo(par, curcurInfo);


	calNBN_Analays(par,
		curGlobalInfo.m_save_dir,
		curGlobalInfo.m_out, curGlobalInfo.m_out_mtx);




}



void getTask(const std::string& filepath, std::vector<std::pair<std::string, int>>& tasks) {
	std::ifstream in(filepath);
	std::pair<std::string, int>cur;
	while (in >> cur.first >> cur.second) {
		tasks.push_back(cur);
	}

	in.close();
}


std::string makeFixedLength(int i, int length)
{
	std::ostringstream ostr;

	if (i < 0)
		ostr << '-';

	ostr << std::setfill('0') << std::setw(length) << (i < 0 ? -i : i);

	return ostr.str();
}






void setThreadTask(utility::GeneralMultiThread& curThread, int taskId) {

	if (taskId == 0) {
		{		int sampleSize = 2e6;
		ThreadInfoNBN curInfo;
		curInfo.m_proname = "MMOP_CEC2015_F07";
		curInfo.m_dim = 16;
		curInfo.m_sample_size = sampleSize;
		curInfo.m_error_name = curInfo.m_proname + "_" + std::to_string(curInfo.m_dim);
		curThread.ms_info_buffer.emplace_back(new ThreadInfoNBN(std::move(curInfo)));
		}

	}

}



void calMultiTask() {
	//	ThreadGlobalNBN totalInfo;
	utility::GeneralMultiThread curThread;
	curThread.ms_cout_flag = true;
	curThread.ms_global_info.reset(new ThreadGlobalNBN);
	ThreadGlobalNBN& totalInfo = dynamic_cast<ThreadGlobalNBN&>(*curThread.ms_global_info);


	m_task_type = TaskType::LonContinous_analyses;
	totalInfo.m_save_dir = "//172.24.24.151/e/DiaoYiya/experiment_data/mmop_multiOpts_nbn_samllsize/";
	//totalInfo.m_save_dir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/mmop_multiOpts_nbn_samllsize/";
	totalInfo.m_save_dir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/mmop_multiOpts_lon_left/";

	totalInfo.m_save_dir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/visualization_lon/";
	totalInfo.m_save_dir = "//172.24.24.151/e/DiaoYiya/experiment_data/test_visualization_lon/";
	//totalInfo.m_save_dir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/visualization_nbn_0821/";

	//totalInfo.m_save_dir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/visualization_nbn_0821/";
	bool m_is_multiTask = true;
	curThread.setNumberThread(80);


	std::filesystem::create_directories(totalInfo.m_save_dir);
	curThread.ms_fun = std::bind(calTask, std::placeholders::_1, std::placeholders::_2);


	setThreadTask(curThread, 6);



	if (!m_is_multiTask) {
		for (auto& it : curThread.ms_info_buffer) {
			calTask(it, curThread.ms_global_info);
		}
	}
	else {
		std::ofstream out(totalInfo.m_save_dir + "error.txt");
		curThread.run(out);
		out.close();

	}


}


int maxSol = 1e6;


std::shared_mutex m_info_mtx;
std::vector<int> solId;
std::shared_mutex m_maxSolId_mtx;
int maxSolId = -1;
std::vector<std::shared_ptr<ofec::SolBase>> sols;
std::vector<double> fitness;



std::string filepath;






struct NbnUpdateInfo {
	std::vector<int> belongs;
	std::vector<double> dis2parent;


	bool update(int curId, int belongId, double dis, ofec::Random* rand) {
		if (fitness[belongId] > fitness[curId]) {
			if (dis < dis2parent[curId]) {
				dis2parent[curId] = dis;
				belongs[curId] = belongId;
				return true;
			}
			else if (dis == dis2parent[curId] && rand->uniform.next() < 0.5) {
				dis2parent[curId] = dis;
				belongs[curId] = belongId;
				return true;
			}
		}

		return false;
	}


	void resize(int size) {
		belongs.resize(size);
		dis2parent.resize(size);


	}

	void init(int size) {
		belongs.resize(size);
		for (int idx(0); idx < belongs.size(); ++idx) {
			belongs[idx] = idx;
		}

		dis2parent.resize(size);
		std::fill(dis2parent.begin(), dis2parent.end(), std::numeric_limits<double>::max());
	}
};



NbnUpdateInfo globalInfo;
std::vector<bool> updateFlag;


NbnUpdateInfo outputInfo;
int outputId = -1;
bool m_isOutputUpdate = false;
std::shared_mutex outputMtx;
int leftUpdate = 0;


void updateNeighbor(int curId,
	ofec::Problem* pro,
	ofec::Random* rand,
	NbnUpdateInfo& curInfo,
	NbnUpdateInfo& curOutputInfo) {
	std::shared_ptr<ofec::SolBase> curSol(pro->createSolution());
	curSol->initialize(pro, rand);
	curSol->evaluate(pro);
	ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	curSol->setFitness(pos * curSol->objective(0));



	{
		std::unique_lock<std::shared_mutex> lock(m_maxSolId_mtx);
		sols[curId] = curSol;
		fitness[curId] = curSol->fitness();

	}

	bool readyFlag = false;
	while (!readyFlag) {
		readyFlag = true;
		{
			std::shared_lock<std::shared_mutex> lock(m_maxSolId_mtx);
			for (int idx(maxSolId + 1); idx <= curId; ++idx) {
				if (sols[idx] == nullptr) {
					readyFlag = false;
					break;
				}
			}
		}

		if (readyFlag)break;
		std::this_thread::sleep_for(std::chrono::seconds(1));
	}



	{
		std::unique_lock<std::shared_mutex> lock(m_maxSolId_mtx);
		maxSolId = std::max(maxSolId, curId);

	}

	{
		std::shared_lock<std::shared_mutex> lock(m_info_mtx);
		curInfo = globalInfo;
	}





	std::vector<int> updateIds;
	updateIds.push_back(curId);


	for (int idx(0); idx < curId; ++idx) {
		double dis = pro->variableDistance(*sols[idx], *curSol);
		if (curInfo.update(idx, curId, dis, rand)) {
			updateIds.push_back(idx);
		}
		else if (curInfo.update(curId, idx, dis, rand)) {
			//		updateIds.push_back(idx);
		}
	}

	{
		std::unique_lock<std::shared_mutex> lock(m_info_mtx);

		for (auto id : updateIds) {
			double dis = curInfo.dis2parent[id];
			int belongId = curInfo.belongs[id];
			globalInfo.update(id, belongId, dis, rand);
		}
		updateFlag[curId] = true;

	}



	bool isOuputUpdate = false;

	{
		std::shared_lock<std::shared_mutex> lock(outputMtx);
		if (curId <= outputId) {
			isOuputUpdate = true;
			curOutputInfo = outputInfo;
		}
	}


	std::vector<int> updateOuputIds;
	updateOuputIds.push_back(curId);
	for (int idx(0); idx < curId; ++idx) {
		double dis = pro->variableDistance(*sols[idx], *curSol);
		if (isOuputUpdate) {
			if (curOutputInfo.update(idx, curId, dis, rand)) {
				updateOuputIds.push_back(idx);
			}
			else if (curOutputInfo.update(curId, idx, dis, rand)) {

			}
		}
	}


	bool isOutput = false;

	int curOutputId = 0;
	if (isOuputUpdate) {
		std::unique_lock<std::shared_mutex> lock(outputMtx);


		for (auto id : updateOuputIds) {
			double dis = curOutputInfo.dis2parent[id];
			int belongId = curOutputInfo.belongs[id];
			outputInfo.update(id, belongId, dis, rand);
		}


		if (--leftUpdate == 0) {
			isOutput = true;
			std::swap(outputInfo, curOutputInfo);


			curOutputId = outputId;
		}
	}

	if (isOutput) {
		std::vector<int> curSolId = solId;
		std::vector<double> curFitness;
		{
			std::shared_lock<std::shared_mutex> lock(m_maxSolId_mtx);
			curFitness = fitness;
		}
		std::cout << "finish outputId\t" << curOutputId << std::endl;
		curSolId.resize(curOutputId + 1);
		curFitness.resize(curOutputId + 1);
		curOutputInfo.resize(curOutputId + 1);
		std::ofstream nbnOut(filepath + "_numSamle_" + std::to_string(curOutputId + 1) + "_nbn.txt");
		ofec::outputNBNdata(nbnOut, curSolId, curOutputInfo.belongs, curOutputInfo.dis2parent, curFitness);
		nbnOut.close();


		{
			std::unique_lock<std::shared_mutex> lock(outputMtx);
			m_isOutputUpdate = false;
		}
	}


}


#include<chrono>
int solIdList = -1;
std::mutex curSolIdMtx;
std::chrono::steady_clock::time_point from;
std::shared_mutex timeMtx;

int maxSample = 3e6;



std::mutex changeOutputUpdateFlagMtx;;
bool changeOutputUpdateFlag = false;



void testAlg() {


	double seed = 0.5;

	using namespace std;
	using namespace ofec;
	ParamMap params;
	params["problem name"] = std::string("CSIWDN");
	params["dataFile1"] = std::string("Net2");
	params["dataFile2"] = std::string("case1");
	params["dataFile3"] = std::string("case11");
	params["use LSTM"] = false;
	auto id_param = std::shared_ptr<const ParamMap>(new ParamMap(params));
	auto pro_ptr = ofec::ADD_PRO(id_param, 0.1);
	auto pro = pro_ptr.get();

	pro->initialize();
	std::shared_ptr<Random> rnd_ptr(new Random(seed));
	//Random rnd(0.5);
	auto rnd = rnd_ptr.get();
	CAST_CSIWDN(pro)->setPhase(CAST_CSIWDN(pro)->totalPhase());



	auto tt = std::chrono::steady_clock::now();




	std::shared_ptr<ofec::SolBase> curSol(pro->createSolution());
	curSol->initialize(pro, rnd);
	curSol->evaluate(pro);
	ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	curSol->setFitness(pos * curSol->objective(0));


	auto to = std::chrono::steady_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(to - tt);

	cout << "花费了 "
		<< double(duration.count()) //* std::chrono::microseconds::period::num 
		<< " micorsecodes" << endl;


	//cout<<"evaluated times\t"<<
}


void runTask(double seed) {
	using namespace std;
	using namespace ofec;
	NbnUpdateInfo curInfo;
	NbnUpdateInfo curOutputInfo;
	ParamMap params;
	params["problem name"] = std::string("CSIWDN");
	params["dataFile1"] = std::string("Net2");
	params["dataFile2"] = std::string("case1");
	params["dataFile3"] = std::string("case11");
	params["use LSTM"] = false;
	auto id_param = std::shared_ptr<const ParamMap>(new ParamMap(params));
	auto pro = ofec::ADD_PRO(id_param, 0.1);
	pro->initialize();
	std::shared_ptr<Random> rnd(new Random(seed));
	//Random rnd(0.5);
	CAST_CSIWDN(pro.get())->setPhase(CAST_CSIWDN(pro.get())->totalPhase());



	while (true) {

		int curSolId = -1;

		{
			std::unique_lock<std::mutex> lock(curSolIdMtx);
			curSolId = ++solIdList;
		}

		if (curSolId >= maxSample) break;



		bool isOutpuReady = false;
		{
			std::shared_lock<std::shared_mutex> lock(timeMtx);
			auto now = std::chrono::steady_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::minutes> (now - from);
			if (duration.count() > 30) {
				isOutpuReady = true;
			}
		}
		bool isOutputFlag = false;
		if (isOutpuReady) {
			std::shared_lock<std::shared_mutex> lock(outputMtx);
			if (!m_isOutputUpdate) {
				std::unique_lock<std::mutex> lock(changeOutputUpdateFlagMtx);
				if (!changeOutputUpdateFlag) {
					isOutputFlag = true;
					changeOutputUpdateFlag = true;
				}
			}

		}


		if (isOutputFlag) {
			std::unique_lock<std::shared_mutex> lock(outputMtx);
			m_isOutputUpdate = true;
			leftUpdate = 0;

			{
				std::shared_lock<std::shared_mutex> lock(m_info_mtx);
				outputInfo = globalInfo;

				for (int idx(outputId + 1); idx <= curSolId; ++idx) {
					if (!updateFlag[idx])++leftUpdate;
				}
			}


			outputId = curSolId;


			std::cout << "begin outputId\t" << outputId << std::endl;
		}


		if (isOutputFlag) {
			std::unique_lock<std::mutex> lock(changeOutputUpdateFlagMtx);
			changeOutputUpdateFlag = false;
		}


		updateNeighbor(curSolId,
			pro.get(),
			rnd.get(),
			curInfo,
			curOutputInfo);



		if (isOutputFlag) {
			{
				std::unique_lock<std::shared_mutex> lock(timeMtx);
				from = std::chrono::steady_clock::now();
			}
		}

	}


}





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


struct NBNdata {
	std::vector<ofec::SolBase> m_sols;
	std::vector<int> m_solIds;
	std::vector<int> m_belongs;
	std::vector<double> m_dis2parent;
	std::vector<double> m_fitness;

	void initialize(int size) {
		m_sols.resize(size);
		m_solIds.resize(size);
		m_belongs.resize(size);
		m_dis2parent.resize(size);
		m_fitness.resize(size);

		for (int idx(0); idx < m_solIds.size(); ++idx) {
			m_solIds[idx] = idx;
			m_belongs[idx] = idx;
			m_dis2parent[idx] = std::numeric_limits<double>::max();
			m_fitness[idx] = -std::numeric_limits<double>::max();
		}
	}


	void outputNBN(std::ostream& out) {
		for (int idx(0); idx < m_solIds.size(); ++idx) {
			out << m_solIds[idx] << "\t" << m_belongs[idx] << "\t" << m_dis2parent[idx] << "\t" < m_fitness[idx] << std::endl;
		}
	}

	void outputSol_CWIDN(std::ostream& out,
		ofec::Problem* pro,
		const ofec::SolBase& sol) {
		
	}
	

	void outputSols_CSWIDN(std::ofstream& out) {
		VarCSIWDN& var = dynamic_cast<Solution<VarCSIWDN> &>(s).variable();
		var.flagLocation() = false;

		for (size_t z = 0; z < m_num_source; z++) {
			if (m_init_type == InitType::kCluster) {
				if (z == 0) {
					int size = m_clusters[m_pop_identifier].size();
					int random = rnd->uniform.nextNonStd<int>(0, size);
					var.index(z) = m_clusters[m_pop_identifier][random] + 1;
				}
				else
					var.index(z) = rnd->uniform.nextNonStd<int>(1, m_num_node + 1);
			}
			else if (m_init_type == InitType::kRandom) {
				var.index(z) = rnd->uniform.nextNonStd<int>(1, m_num_node + 1);
			}
			else if (m_init_type == InitType::kDistance) {
				std::vector<Real> roulette(m_num_node + 1, 0);
				for (size_t i = 0; i < roulette.size() - 1; ++i) {
					roulette[i + 1] = roulette[i] + m_distance_node[i];
				}
				Real random_value = rnd->uniform.nextNonStd<Real>(roulette[0], roulette[roulette.size() - 1]);
				for (size_t i = 0; i < roulette.size() - 1; ++i) {
					if (random_value >= roulette[i] && random_value < roulette[i + 1]) {
						var.index(z) = i + 1;
						break;
					}
				}
			}
			else if (m_init_type == CSIWDN::InitType::kBeVisited) {
				if (z <= endSourceIdx()) {
					int identifier = popIdentifier();
					int size = clusters()[identifier].size();
					std::vector<Real> pro(size, 0.0);
					for (size_t i = 0; i < size; ++i) {
						pro[i] = getProBeVisited()[z][clusters()[identifier][i]];
					}
					var.index(z) = (rnd->uniform.spinWheel(pro.begin(), pro.end()) - pro.begin()) + 1;
				}
				else {
					var.index(z) = rnd->uniform.nextNonStd<int>(1, numberNode() + 1);
				}
			}
			else throw MyExcept("No this initialization type");

			var.source(z) = 1.0;
			var.startTime(z) = const_cast<VarCSIWDN&>(dynamic_cast<Optima<VarCSIWDN>&>(*m_optima).variable(0)).startTime(z);
			var.duration(z) = const_cast<VarCSIWDN&>(dynamic_cast<Optima<VarCSIWDN>&>(*m_optima).variable(0)).duration(z);
			size_t size = var.duration(z) / m_pattern_step;
			if (var.duration(z) % m_pattern_step != 0)
				++size;
			var.multiplier(z).resize(size);
			for (auto& j : var.multiplier(z)) {
				j = rnd->uniform.nextNonStd<float>(m_min_multiplier, m_max_multiplier);
			}

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
		eval_fun(*sols[id], pro);

		//for(int idx())
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
			//sols[id]->initialize(pro, rnd);
			eval_fun(*sols[idx], pro);
		}
		//for(int idx())
	}
}



void calNBN_iterater(ofec::Problem* pro, 
	const std::vector<int>& subArea, int& bestId,
	std::vector<std::shared_ptr<ofec::SolBase>>& sols,
	std::vector<int>& belong,
	std::vector<int>& dis2parent) {
	std::vector<int> sortIdx(subArea);
	std::sort(sortIdx.begin(), sortIdx.end(), [&](int a, int b) {
		return sols[a]->fitness() < sols[b]->fitness();
	});
	bestId = sortIdx.back();
	for (int idx(0); idx < sortIdx.size(); ++idx) {
		int ida = sortIdx[idx];
		for (int idy(idx + 1); idy < sortIdx.size(); ++idy) {
			int idb = sortIdx[idy];
			double dis = sols[ida]->variableDistance(*sols[idb], pro);
			if (dis < dis2parent[ida]) {
				belong[ida] = idb;
				dis2parent[ida] = dis;
			}
		}
	}

}

void calNBN_iterater_ThreadTask(
	int from, int to,
	ofec::Problem* pro,
	const std::vector<std::vector<int>>& subAreas, std::vector<int>& bestId,
	std::vector<std::shared_ptr<ofec::SolBase>>& sols,
	std::vector<int>& belong,
	std::vector<int>& dis2parent
) {
	for (int idSub(from); idSub < to; ++idSub) {
		calNBN_iterater(pro, subAreas[idSub], bestId[idSub], sols, belong, dis2parent);
	}
}




void getNodeState(int m, std::vector<std::vector<bool>>& selected, std::vector<bool>&cur, int from) {
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
	const std::vector<bool> & state, const ofec::SolBase& sol) {
	int divisionId = stateId * numDivision;
	int curId = 0;
	auto& cursol = dynamic_cast<const ofec::CSIWDN::solutionType&>(sol);
	
	for (int idx(0); idx < state.size(); ++idx) {
		if (state[idx]) {
			curId = curId * numNode + cursol.variable().index(idx);
		}
	}
	return divisionId + curId;
}

void calNBN_division(ofec::Problem* pro,
	const std::vector<int>& subArea, int& bestId,
	std::vector<std::shared_ptr<ofec::SolBase>>& sols,
	std::vector<int>& belong,
	std::vector<int>& dis2parent) {

	
	int selectedSource = 0;
	int numDivision = 1;
	while (selectedSource < CAST_CSIWDN(pro)->numSource() &&  numDivision * CAST_CSIWDN(pro)->numberNode()<subArea.size()) {
		selectedSource++;
		numDivision *= CAST_CSIWDN(pro)->numberNode();
	}

	std::vector<bool> curState(CAST_CSIWDN(pro)->numSource(),false);
	std::vector<std::vector<bool>> totalState;
	std::vector<std::vector<int>> divisions;
	std::vector<int> bestSolIds = subArea;
	std::vector<bool> visited(subArea.size(), false);
	for (int numSelect(selectedNode); numSelect >0; --numSelect) {
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
					totalState[stateId], sols[solId]);
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
				std::ref(sols), std::ref(belong), std::ref(dis2parent)));
		}
		for (auto& thrd : thrds)
			thrd.join();
		
	}
	

}


// generate multiSols
void generateSols(ofec::Problem* pro, ofec::Random* rnd,
	std::function<void(ofec::SolBase& sol, ofec::Problem* pro)>& fun,
	NBNdata& data,
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
	std::vector<std::thread> thrds;
	int num_task = std::thread::hardware_concurrency();
	int num_samples = numSols;
	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);

	std::vector<std::shared_ptr<ofec::Random>> rnds(num_task);
	std::vector<std::shared_ptr<ofec::Problem>> pros(num_task);
	for (int idx(0); idx < rnds.size(); ++idx) {
		rnds[idx] = std::make_shared<ofec::Random>(ofec::Random(rnd->uniform.next()));

		auto& cur_pro = pros[idx];
		auto param = ADD_PARAM(pro->parameters());
		cur_pro = ofec::ADD_PRO(param, rnd->uniform.next());
		cur_pro->initialize();
		//Random rnd(0.5);
		CAST_CSIWDN(cur_pro.get())->setPhase(CAST_CSIWDN(cur_pro.get())->totalPhase());
	}


	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];
		thrds.push_back(std::thread(
			updateRandSol, pros[i].get(), rnds[i].get(), std::ref(fun), subSols,
			from_to.first, from_to.second, std::ref(sols)));
	}
	for (auto& thrd : thrds)
		thrd.join();



	num_samples = sols.size();
	UTILITY::assignThreads(num_samples, num_task, tasks);


	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];
		thrds.push_back(std::thread(
			updateRandSolSubRegion, pros[i].get(), rnds[i].get(), std::ref(fun), subSols,
			from_to.first, from_to.second, std::ref(sols)));
	}
	for (auto& thrd : thrds)
		thrd.join();

	data.updateFitness();

}



void generateSols(ofec::Problem* pro, ofec::Random* rnd,
	NBNdata& data,
	int numSols = 2e4,int subSols = 100) {



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

	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];
		thrds.push_back(std::thread(
			calNBN_iterater_ThreadTask, from_to.first, from_to.second,
			pro, std::cref(subRegion), std::ref(bestSolIds),
			std::ref(data.m_sols), std::ref(data.m_belongs), std::ref(data.m_dis2parent)));
	}
	for (auto& thrd : thrds)
		thrd.join();


	int bestId(0);

	calNBN_division(pro, subRegion, bestId, data.m_sols,
		data.m_belongs, data.m_dis2parent);
}


void testAlg22() {
	int id, belong;
	double x, y, z;



	std::vector<int> belongs;
	std::vector<int> ids;
	std::vector<std::vector<int>> sons;

	std::string filepath = "//172.24.24.151/f/nbn/GTOP-Messenger_7000000.txt";

	std::stringstream buf;
	UTILITY::readFileBuffer(filepath, buf);

	//std::ifstream in(filepath);


	std::string str = buf.str();

	int stop22 = -1;

	while (buf >> id >> belong >> x >> y >> z) {
		ids.push_back(id);
		belongs.push_back(belong);
	}


	std::cout << "data size\t" << ids.size() << std::endl;


	std::vector<int> header;

	sons.resize(belongs.size());
	for (int idx(0); idx < sons.size(); ++idx) {
		if (belongs[idx] == idx) {
			header.push_back(idx);
		}
		else {
			sons[belongs[idx]].push_back(idx);
		}
	}

	std::vector<bool> visited(belongs.size(), false);

	std::queue<int> queV;
	for (auto& it : header) {
		queV.push(it);
	}

	while (!queV.empty()) {

		int value = queV.front();
		queV.pop();
		visited[value] = true;
		for (auto& it : sons[value]) {
			queV.push(it);
		}
	}

	for (auto& it : visited) {
		if (!it) {
			std::cout << "error" << std::endl;
			break;
		}
	}

	int stop = -1;

}

void readFile() {
	FILE* fp;
	char buf[1024];
	std::string filepath = "//172.24.24.151/f/nbn/GTOP-Messenger_7000000.txt";
	// 打开文件
	fp = fopen(filepath.c_str(), "r");
	if (fp == NULL)
	{
		printf("open file failed!\n");
		return;
	}
	int id, belong;
	double x, y, z;
	std::vector<int> belongs;
	std::vector<int> ids;

	// 一行一行读取文件内容
	while (fgets(buf, 1024, fp) != NULL)
	{

		std::stringstream ss;
		ss << buf;

		ss >> id >> belong >> x >> y >> z;
		{
			ids.push_back(id);
			belongs.push_back(belong);
		}
		//printf("%s", buf);
	}

	std::cout << "data size\t" << ids.size() << std::endl;


	std::vector<int> header;
	std::vector<std::vector<int>> sons;
	sons.resize(belongs.size());
	for (int idx(0); idx < sons.size(); ++idx) {
		if (belongs[idx] == idx) {
			header.push_back(idx);
		}
		else {
			sons[belongs[idx]].push_back(idx);
		}
	}

	std::vector<bool> visited(belongs.size(), false);

	std::queue<int> queV;
	for (auto& it : header) {
		queV.push(it);
	}

	while (!queV.empty()) {

		int value = queV.front();
		queV.pop();
		visited[value] = true;
		for (auto& it : sons[value]) {
			queV.push(it);
		}
	}

	for (auto& it : visited) {
		if (!it) {
			std::cout << "error" << std::endl;
			break;
		}
	}

	int stop = -1;
}

namespace ofec {

	void registerParamAbbr();
	void customizeFileName();

	void run() {
		registerParamAbbr();
		customizeFileName();
		using namespace std;


		//readFile();
		//testAlg22();



	//	//testAlg();
	//	//int stop = -1;


	////	filepath = "//172.24.24.151/e/DiaoYiya/experiment_data/cswidn_nbn_data/cwidn_nbn";

	//	filepath = "/home/lab408/share//Student/2018/YiyaDiao/NBN_data/cswidn_nbn_data4/cwidn_nbn2";
	//	//	filepath = "E:/DiaoYiya/experiment_data/cswidn_nbn_data2/cwidn_nbn";
	//	std::filesystem::create_directories(filepath);

	//	globalInfo.init(maxSample);
	//	solId.resize(maxSample);
	//	for (int idx(0); idx < solId.size(); ++idx) {
	//		solId[idx] = idx;
	//	}
	//	sols.resize(maxSample);
	//	std::fill(sols.begin(), sols.end(), nullptr);
	//	fitness.resize(maxSample, 0);
	//	updateFlag.resize(maxSample, false);
	//	std::fill(updateFlag.begin(), updateFlag.end(), false);

	//	std::vector<std::thread> thrds;
	//	int numThread = std::thread::hardware_concurrency();

	//	//numThread = 5;
	//	for (int thrdId(0); thrdId < numThread; ++thrdId) {
	//		thrds.emplace_back(runTask, double(thrdId + 1) / double(numThread + 1));
	//	}


	//	for (auto& thrd : thrds) {
	//		thrd.join();
	//	}


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
