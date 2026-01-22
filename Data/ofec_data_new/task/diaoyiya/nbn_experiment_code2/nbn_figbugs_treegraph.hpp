
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

#include <filesystem>
#include <string>
#include <vector>
#include <iostream>


void runTest(int runId) {

	std::cout << runId << std::endl;

	while (true) {
		int a;
		a++;
	}
}



void test() {

	int id, belong;
	double x, y, z;



	std::vector<int> belongs;
	std::vector<int> ids;
	std::vector<std::vector<int>> sons;

	std::string filepath = "F:/nbn/GTOP-Cassini2_7000000.txt";

	std::stringstream buf;
	UTILITY::readFileBuffer(filepath, buf);

	//std::ifstream in(filepath);


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


	filepath = "F:/nbn/GTOP-Messenger_7000000.txt";
	// 打开文件
	fp = fopen(filepath.c_str(), "r");
	if (fp == NULL)
	{
		printf("open file failed!\n");
		return;
	}
	int id, belong;
	double x, y, z;
	double dis(0);
	std::vector<int> ids;
	std::vector<int> belongs;
	std::vector<double> fitness;
	std::vector<double> dis2parent;


	// 一行一行读取文件内容
	while (fgets(buf, 1024, fp) != NULL)
	{

		std::stringstream ss;
		ss << buf;

		ss >> id >> belong >> x >> y >> z >> dis;
		{
			ids.push_back(id);
			belongs.push_back(belong);
			fitness.push_back(z);
			dis2parent.push_back(dis);
		}
		//printf("%s", buf);
	}

	std::cout << "data size\t" << ids.size() << std::endl;


	std::vector<int> header;
	std::vector<std::vector<int>> sons;
	sons.resize(belongs.size());
	int valueErrorIdx = 0;


	int exeed_data(0);


	for (int idx(0); idx < sons.size(); ++idx) {
		if (belongs[idx] == sons.size()) {
			belongs[idx] = idx;
			++exeed_data;
		}
	}





	std::cout << "exceed_data\t" << exeed_data << std::endl;


	for (int idx(0); idx < sons.size(); ++idx) {
		if (belongs[idx] == idx) {
			header.push_back(idx);
		}
		else {
			sons[belongs[idx]].push_back(idx);


			if (fitness[belongs[idx]] < fitness[idx]) {
				++valueErrorIdx;
			}
		}
	}
	if (valueErrorIdx) {
		std::cout << "error value numIdx\t" << valueErrorIdx << std::endl;
	}
	std::vector<int> visited(belongs.size(), 0);
	std::queue<int> queV;
	for (auto& it : header) {
		queV.push(it);
	}

	while (!queV.empty()) {

		int value = queV.front();
		queV.pop();
		++visited[value];
		for (auto& it : sons[value]) {
			queV.push(it);
		}
	}


	int totalNotVisited(0);

	for (auto& it : visited) {
		if (it != 1) {
			std::cout << "error" << std::endl;
			break;
		}
	}



	for (auto& it : visited) {
		if (!it) {
			++totalNotVisited;
			//break;
		}
	}

	std::cout << "not visited number\t" << totalNotVisited << std::endl;

	int stop = -1;


	std::vector<int> subset(exeed_data, 0);

	for (int idx(0); idx < header.size(); ++idx) {
		int head = header[idx];
		queV.push(head);

		while (!queV.empty()) {

			int value = queV.front();
			queV.pop();
			++subset[idx];

			//	visited[value] = true;
			for (auto& it : sons[value]) {
				queV.push(it);
			}


		}
	}


	std::cout << "subset size\t" << std::endl;
	for (auto& it : header) {
		std::cout << it << std::endl;
	}
}


void readNBNdata(
	const std::string& filepath, 
	std::vector<ofec::TreeGraphSimple::NBNdata> & nbn_data) {
	FILE* fp;
	char buf[1024];
	// 打开文件
	fp = fopen(filepath.c_str(), "r");
	if (fp == NULL)
	{
		printf("open file failed!\n");
		return;
	}
	double x(0), y(0), z(0);
	ofec::TreeGraphSimple::NBNdata curdata;
	nbn_data.clear();
	// 一行一行读取文件内容
	while (fgets(buf, 1024, fp) != NULL)
	{

		std::stringstream ss;
		ss << buf;

		ss >> curdata.m_id >> curdata.m_belong>> x >> y >> curdata.m_fitness >> curdata.m_dis2parent;
		curdata.m_originId = curdata.m_id;
		{
			nbn_data.push_back(curdata);
		}
		//printf("%s", buf);
	}

}


void filterSize(std::vector<ofec::TreeGraphSimple::NBNdata>& data) {

	bool flagAdd(false);
	for (auto& it : data) {
		if (it.m_belong == data.size()) {
			flagAdd = true;
			break;
		}
	}
	double maxFit(-std::numeric_limits<double>::max());
	double minFit(std::numeric_limits<double>::min());
	for (auto& it : data) {
		maxFit = std::max(it.m_fitness, maxFit);
		minFit = std::min(it.m_fitness, minFit);
	}
	maxFit = minFit + (maxFit - minFit) * 1.01;
	
	if (flagAdd) {
		ofec::TreeGraphSimple::NBNdata curdata;
		curdata.m_id = curdata.m_originId = data.size();
		curdata.m_belong = curdata.m_id;
		curdata.m_dis2parent = std::numeric_limits<double>::max();
		curdata.m_fitness = maxFit;
		data.push_back(curdata);
	}
}

namespace ofec {

	void registerParamAbbr();
	void customizeFileName();
	void run() {
		registerParamAbbr();
		customizeFileName();
		using namespace std;


		std::string filepath = "//172.24.24.151/f/nbn/GTOP-Messenger_7000000.txt";
		auto dir = "//172.24.207.203/share/student/2017/wangjunche/nbn/";
		dir = "E:/nbn_data/check_nbn_data/";
		std::string  filename = "GTOP-Rosetta_7000000";
		filepath = dir + filename + ".txt";
		std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
		readNBNdata(filepath,nbn_data);
		filterSize(nbn_data);
		TreeGraphSimple nbn_graph;
		nbn_graph.setNBNdata(nbn_data);
		nbn_graph.calNetwork();


		std::string savepath = dir + filename + "_network.txt";
		std::ofstream out(savepath);
		nbn_graph.outputNBNnetwork(out);
		out.close();
		//std::vector<std::thread> thrds;
		//int num_task = 160;
		////	num_task = 1;
		//for (size_t i = 0; i < num_task; ++i) {
		//	thrds.push_back(std::thread(
		//		&runTest,  i));
		//}
		//for (auto& thrd : thrds)
		//	thrd.join();


	//	calMultiTask();
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
