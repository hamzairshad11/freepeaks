

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

