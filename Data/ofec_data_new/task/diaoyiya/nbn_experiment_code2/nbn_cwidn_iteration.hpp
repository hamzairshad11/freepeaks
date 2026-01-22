

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

	cout << "»¨·ÑÁË "
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



