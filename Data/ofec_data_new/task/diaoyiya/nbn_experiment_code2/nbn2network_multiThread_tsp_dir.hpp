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





class ThreadInfoNBN :public utility::GeneralMultiThreadInfo {
public:
	//std::ofstream& m_errOut;
	//std::mutex& m_errOut_mtx;
	std::pair<std::string, std::string> m_filename;
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

std::string subfix = ".nbn";

void readDir(const std::string& dir, std::vector<std::pair<std::string, std::string>>& names) {


	std::pair<std::string, std::string> dir_name;

	const std::string filepath = dir;
	const std::filesystem::path sandbox{ filepath };
	names.clear();
	// directory_iterator can be iterated using a range-for loop
	for (auto const& dir_entry : std::filesystem::directory_iterator{ sandbox })
	{
		auto curpath = std::filesystem::path(dir_entry.path());
		//std::cout << curpath.filename().string() << std::endl;
		auto curname = curpath.filename().string();

		const std::filesystem::path sandbox2{ curpath };

		if (std::filesystem::is_directory(curpath)) {
			for (auto const& dir_entry : std::filesystem::directory_iterator{ sandbox2 })
			{
				auto curpath = std::filesystem::path(dir_entry.path());
				auto curname2 = curname + '/' + curpath.filename().string();
				if (curname2.find(subfix) != std::string::npos) {
					int idx = curname2.find(subfix);
					auto subname = curname2.substr(0, idx);
					std::cout << "subname\t" << subname << std::endl;
					dir_name.first = subname;

					curname2 = curpath.filename().string();
					idx = curname2.find(subfix);
					subname = curname2.substr(0, idx);

					dir_name.second = subname;
					names.push_back(dir_name);
				}
			}
		}
	}
}





void calTask(std::unique_ptr<utility::GeneralMultiThreadInfo>& curInfo,
	std::unique_ptr<utility::GeneralMultiThreadInfo>& globalInfo) {
	auto& curGlobalInfo = dynamic_cast<ThreadGlobalNBN&>(*globalInfo);
	auto& curcurInfo = dynamic_cast<ThreadInfoNBN&>(*curInfo);


	{

		std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
		inputNBNdata(curGlobalInfo.m_read_dir + curcurInfo.m_filename.first + subfix,
			nbn_data);
		//ofec::filterNBNDataByDensity(outputDir, filename, nbn_data);
		ofec::TreeGraphSimple nbn_graph;
		//	nbn_graph.m_task_name = curcurInfo.getBBOBfilename();

		nbn_graph.setNBNdata(nbn_data);
		nbn_graph.calNetwork();
		std::ofstream out(curGlobalInfo.m_save_dir + curcurInfo.m_filename.second + "_network.txt");
		nbn_graph.outputNBNnetwork(out);
		out.close();


		ouputNBNdata(curGlobalInfo.m_save_dir + curcurInfo.m_filename.second + "_nbn.txt", nbn_graph.get_nbn_data());
	}

}



void setThreadTask(const std::string& dir, utility::GeneralMultiThread& curThread) {
	std::vector<std::pair<std::string, std::string>> names;
	readDir(dir, names);

	ThreadInfoNBN curInfo;

	for (const auto& ssname : names) {
		const auto& filename = ssname;

		curInfo.m_filename = ssname;
		curInfo.m_taskname = filename.second;
		curInfo.m_error_name = "nbn_" + filename.second;
		curThread.ms_info_buffer.emplace_back(new ThreadInfoNBN(std::move(curInfo)));

	}
}




void calMultiTask() {
	//	ThreadGlobalNBN totalInfo;
	utility::GeneralMultiThread curThread;
	curThread.ms_cout_flag = true;
	curThread.ms_global_info.reset(new ThreadGlobalNBN);
	ThreadGlobalNBN& totalInfo = dynamic_cast<ThreadGlobalNBN&>(*curThread.ms_global_info);




	std::string dir = "//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/nbn_visualization_onemax2/";

	dir = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data_paper_com/oneMax/onemaxNeighborData/";
	//"/mnt/Data/"
	dir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data_paper_com/oneMax/onemaxNeighborData/";

	dir = "//172.24.207.203/share/2018/diaoyiya/paper_com_experiment_data/tsp/tsp2/BCL380/";

	dir = "//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/tsp_data/test2/";
	totalInfo.m_read_dir = dir;
	totalInfo.m_save_dir = "//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/nbn_visualization_onemax_network2/";
	totalInfo.m_save_dir = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data_paper_com/oneMax/onemaxNeighborFilterNetwork2/";
	totalInfo.m_save_dir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data_paper_com/oneMax/onemaxNeighborFilterNetwork2/";

	totalInfo.m_save_dir = "//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/tsp_data/tsp_figure_data2/";

	bool m_is_multiTask = true;
	curThread.setNumberThread(std::thread::hardware_concurrency());


	std::filesystem::create_directories(totalInfo.m_save_dir);
	curThread.ms_fun = std::bind(calTask, std::placeholders::_1, std::placeholders::_2);


	setThreadTask(dir, curThread);


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

		//std::vector<std::string> filenames;
		//std::string dir = "//172.29.65.56/share/Student/2018/YiyaDiao/NBN_data/visualization_nbn_3/";

		//readDir(dir, filenames);

		calMultiTask();


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
