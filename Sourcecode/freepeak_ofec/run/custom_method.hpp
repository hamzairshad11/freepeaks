#ifndef OFEC_CUSTOM_METHOD_HPP
#define OFEC_CUSTOM_METHOD_HPP

#include "../core/global.h"
#include "interface.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>

#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../instance/problem/continuous/free_peaks/subproblem/subproblem.h"
#include "../instance/problem/continuous/free_peaks/subproblem/function/one_peak_function.h"
#include "../instance/problem/continuous/free_peaks/factory.h"
#include "../instance/problem/continuous/free_peaks/free_peaks_multiparty.h"
#include "mpmcoea_solver.hpp"

namespace fs = std::filesystem;

// Maps suite_id → "P01"-"P12" label (P01-P08 = report, P09-P12 = supplementary)
inline std::string suitePLabel(int suite_id) {
    static const std::unordered_map<int, int> m{
        {1,1},{2,2},{4,3},{5,4},{9,5},{10,6},{11,7},{12,8},
        {3,9},{6,10},{7,11},{8,12}
    };
    auto it = m.find(suite_id);
    int n = it != m.end() ? it->second : suite_id;
    return "P" + std::string(n < 10 ? "0" : "") + std::to_string(n);
}


/*void printBenchmarkStructure(ofec::FreePeaks* fp) {
    std::cout << "\n=== Benchmark Structure ===\n";
    auto tree = fp->subspaceTree().treeData();
    for (const auto& [parent, children] : tree) {
        std::cout << parent << ":\n";
        for (const auto& [child, ratio] : children) {
            std::cout << "  ├─ " << child << " (weight: " << ratio << ")\n";
            auto subpro = fp->subproblem(child);
            if (subpro) {
                auto func = subpro->function();
                if (auto opf = dynamic_cast<OnePeakFunction*>(func)) {
                    std::cout << "    └─ Peak count: " << opf->numObjs()
                        << ", center: [" << opf->optimalVars()[0][0] << ", " << opf->optimalVars()[0][1] << "]\n";
                }
            }
        }
    }
}*/


void testFreePeak() {

    using namespace ofec;

    ParameterMap params;
    params["generation"] = std::string("read_file");
    params["dataFile1"] = std::string("sop/2d_2_s1_1_1.txt");


    std::string freepeakName = "free_peaks";
    std::shared_ptr<Random> rnd(new Random(0.5));
    std::shared_ptr<Environment> env(generateEnvironmentByFactory(freepeakName));
    env->recordInputParameters();
    env->initialize();

    env->setProblem(generateProblemByFactory(freepeakName));
    env->problem()->inputParameters().input(params);
    env->problem()->recordInputParameters();
    env->initializeProblem(0.5);


    std::shared_ptr<ofec::SolutionBase> cursol(env->problem()->createSolution());;
    cursol->initialize(env.get(), rnd.get());
    cursol->evaluate(env.get(), false);



}

void outputProMesh(const std::string& filepath, ofec::Environment* env, int div = 400) {

    using namespace ofec;

    auto con_pro = CAST_CONOP(env->problem());
    std::shared_ptr<Random> rnd(new Random(0.5));
    std::vector<Solution<>> sols(div * div, Solution<>(con_pro->numberObjectives(), con_pro->numberConstraints(), con_pro->numberVariables()));
    //std::vector<std::shared_ptr<SolutionBase>> solbases;

    auto eval_fun =
        [](ofec::SolutionBase& sol, ofec::Environment* env) {
        using namespace ofec;
        sol.evaluate(env, false);
        //pro->evaluate(sol.variableBase(), sol.objective(), sol.constraint());
        ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
        sol.setFitness(pos * sol.objective(0));

        };

    for (int idx(0); idx < div; ++idx)
        for (int jdx(0); jdx < div; ++jdx) {
            sols[idx * div + jdx].variable().vector()[0] = con_pro->range(0).first + (con_pro->range(0).second - con_pro->range(0).first) * double(idx) / double(div - 1);
            sols[idx * div + jdx].variable().vector()[1] = con_pro->range(1).first + (con_pro->range(1).second - con_pro->range(1).first) * double(jdx) / double(div - 1);
            eval_fun(sols[idx * div + jdx], env);
            //fitness.push_back(it.fitness());
        }


    //auto range = CAST_CONOP(env->problem())->range();

    //std::vector<std::vector<double>*> objectives(sols.size());
    //for (int idx(0); idx < objectives.size(); ++idx) {
    //    objectives[idx] = &sols[idx].objective();
    //}
    //std::vector<int> rank(sols.size());
    ////	nd_sort::fastSort(objectives, rank, env->problem()->optimizeMode());
    //nd_sort::filterSort(objectives, rank, env->problem()->optimizeMode());

    std::stringstream out;
    out << std::fixed << std::setprecision(6);
    for (int idx(0); idx < div; ++idx)
        for (int jdx(0); jdx < div; ++jdx) {
            int curId = idx * div + jdx;
            out << sols[curId].variable().vector()[0] << " " << sols[curId].variable().vector()[1] << " " << sols[curId].fitness() << "\n";
            //sols[idx * div + jdx].variable().vector()[0] = con_pro->range(0).first + (con_pro->range(0).second - con_pro->range(0).first) * double(idx) / double(div - 1);
            //sols[idx * div + jdx].variable().vector()[1] = con_pro->range(1).first + (con_pro->range(1).second - con_pro->range(1).first) * double(jdx) / double(div - 1);
            //eval_fun(sols[idx * div + jdx], env.get());
            //fitness.push_back(it.fitness());
        }

    std::ofstream outFile(filepath + ".txt");
    outFile << out.rdbuf();
    outFile.close();
    //out.close();

    std::cout << "generated\t" << filepath << std::endl;
}


void testEqual(ofec::Environment* env1, ofec::Environment* env2, ofec::Random* rnd, int numTest = 1e3) {

    using namespace ofec;
    auto domain = CAST_CONOP(env1->problem())->domain();

    for (int idx(0); idx < domain.size(); ++idx)
        std::cout << domain[idx].limit.first << "\t" << domain[idx].limit.second << std::endl;

    std::cout << "testEqual start" << std::endl;

    std::cout << "\n\n\n";
    for (int i = 0; i < numTest; ++i) {
        std::shared_ptr<ofec::SolutionBase> sol1(env1->problem()->createSolution());
        sol1->initialize(env1, rnd);
        sol1->evaluate(env1, false);

        std::shared_ptr<ofec::SolutionBase> sol2(env2->problem()->createSolution());
        sol2->initialize(env2, rnd);
        dynamic_cast<ofec::Solution<>&>(*sol2).variable() = dynamic_cast<ofec::Solution<>&>(*sol1).variable();


        sol2->evaluate(env2, false);

        // Compare solutions
        if (sol1->objective() != sol2->objective()) {
            auto& solvec1 = dynamic_cast<ofec::Solution<>&>(*sol1).variable().vector();
            auto& solvec2 = dynamic_cast<ofec::Solution<>&>(*sol2).variable().vector();

            std::cout << "different var" << std::endl;
            for (int idx(0); idx < solvec1.size(); ++idx) {
                std::cout << solvec1[idx] << "\t";
            }
            std::cout << std::endl;


            for (int idx(0); idx < solvec2.size(); ++idx) {
                std::cout << solvec2[idx] << "\t";
            }
            std::cout << std::endl;

            //<< "\t" << sol2->variableBase().toString() << std::endl;

            std::cout << "different obj \t" << sol1->objective()[0] << "\t" << sol2->objective()[0] << std::endl;
            std::cout << "Solutions are not equal!" << std::endl;


        }
    }

    std::cout << "testEqual end" << std::endl;
}



static bool compareTreeData(
    const std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>& a,
    const std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>& b,
    double eps = 1e-9)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].first != b[i].first) return false;
        const auto& va = a[i].second;
        const auto& vb = b[i].second;
        if (va.size() != vb.size()) return false;
        for (size_t j = 0; j < va.size(); ++j) {
            if (va[j].first != vb[j].first) return false;
            if (std::fabs(va[j].second - vb[j].second) > eps) return false;
        }
    }
    return true;
}



void outputData(
    const std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>& a)
{
    for (const auto& [parent, children] : a) {
        std::cout << parent << ":\n";
        for (const auto& [child, value] : children) {
            std::cout << "  " << child << " -> " << value << "\n";
        }
        std::cout << "\n";
    }
}




void gridSampling(ofec::Environment* env, const std::string& filename,
    int grid_size = 1e2) {

    using namespace ofec;


    auto domain = CAST_CONOP(env->problem())->domain();
    int dim = domain.size();

    if (dim != 2) {
        std::cerr << "Error: This function only supports 2D problems!" << std::endl;
        std::cerr << "Current dimension: " << dim << std::endl;
        return;
    }

    std::cout << "\n=== Problem Domain Information ===" << std::endl;
    for (int idx = 0; idx < dim; ++idx) {
        std::cout << "Variable " << idx << ": ["
            << domain[idx].limit.first << ", "
            << domain[idx].limit.second << "]" << std::endl;
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    outfile << std::fixed << std::setprecision(10);

    outfile << "x1,x2,f(x1,x2)" << std::endl;

    double x1_min = domain[0].limit.first;
    double x1_max = domain[0].limit.second;
    double x2_min = domain[1].limit.first;
    double x2_max = domain[1].limit.second;

    double dx1 = (x1_max - x1_min) / (grid_size - 1);
    double dx2 = (x2_max - x2_min) / (grid_size - 1);

    std::cout << "\n=== Grid Sampling Parameters ===" << std::endl;
    std::cout << "Grid size: " << grid_size << "x" << grid_size
        << " = " << grid_size * grid_size << " points" << std::endl;
    std::cout << "x1 step: " << dx1 << std::endl;
    std::cout << "x2 step: " << dx2 << std::endl;

    std::shared_ptr<ofec::SolutionBase> sol(env->problem()->createSolution());
    auto& solvec = dynamic_cast<ofec::Solution<>&>(*sol).variable().vector();
    int total_points = grid_size * grid_size;
    int current_point = 0;

    std::cout << "\n=== Start Grid Sampling ===" << std::endl;

    for (int i = 0; i < grid_size; ++i) {
        double x1 = x1_min + i * dx1;

        for (int j = 0; j < grid_size; ++j) {
            double x2 = x2_min + j * dx2;

            //std::vector<double> vars = { x1, x2 };

            solvec = { x1,x2 };
            sol->evaluate(env, false);

            outfile << x1 << "," << x2 << "," << sol->objective().front() << std::endl;

            current_point++;
            if (current_point % (total_points / 10) == 0) {
                std::cout << "Progress: "
                    << (current_point * 100 / total_points)
                    << "% completed" << std::endl;
            }
        }
    }

    outfile.close();

    std::cout << "\n=== Grid Sampling Completed ===" << std::endl;
    std::cout << "Data saved to: " << filename << std::endl;
    std::cout << "Total points: " << total_points << std::endl;
}


bool equalFreePeaks(ofec::Problem* pro1, ofec::Problem* pro2) {
    using namespace ofec;
    using namespace free_peaks;
    auto lhs = CAST_FPs(pro2);
    auto lkh = CAST_FPs(pro1);


    // sizes
    if (lhs->numberVariables() != lkh->numberVariables()) return false;
    if (lhs->numberObjectives() != lkh->numberObjectives()) return false;
    if (lhs->numberConstraints() != lkh->numberConstraints()) return false;

    // generation type / file name
    if (lhs->generationType() != lkh->generationType()) return false;
    if (lhs->fileName() != lkh->fileName()) return false;

    // KD-tree structure (names + ratios)
    auto ta = lhs->subspaceTree().treeData();
    auto tb = lkh->subspaceTree().treeData();


    std::cout << "output TreeData\t" << std::endl;
    outputData(ta);
    outputData(tb);

    if (!compareTreeData(ta, tb)) return false;

    // Compare subproblem archived parameters for each named box
    const auto& mapA = lhs->subspaceTree().name_box_subproblem;
    const auto& mapB = lkh->subspaceTree().name_box_subproblem;

    if (mapA.size() != mapB.size()) return false;

    for (const auto& p : mapA) {
        const std::string& name = p.first;
        auto itB = mapB.find(name);
        if (itB == mapB.end()) return false;

        const auto& subA_ptr = p.second.second;
        const auto& subB_ptr = itB->second.second;

        if ((subA_ptr == nullptr) != (subB_ptr == nullptr)) return false;

        if (subA_ptr && subB_ptr) {
            // compare archived parameters (ParameterMap has operator==)
            const ofec::ParameterMap& pmA = subA_ptr->archivedParameters();
            const ofec::ParameterMap& pmB = subB_ptr->archivedParameters();
            if (!(pmA == pmB)) return false;
        }
    }

    return true;
}



void generateHandMadeFreePeak(const std::string& funname) {

    using namespace ofec;
    using namespace free_peaks;

    std::string dirname = "multiparty_multimodal/";
    int numDim = 2, numObj = 1, numCon = 0;

    std::shared_ptr<ofec::Random> rnd(new Random(0.5));


    FreePeaks::registerFP();
    std::string freepeakName = "free_peaks";
    std::shared_ptr<Environment> env(generateEnvironmentByFactory(freepeakName));
    env->recordInputParameters();
    env->initialize();
    env->setProblem(generateProblemByFactory(freepeakName));

    auto freepeak = CAST_FPs(env->problem());
    //std::shared_ptr<FreePeaks> freepeak (FreePeaks::create());

    ParameterMap freepeak_param;
    freepeak_param["generation_type"] = std::string("assigned");
    // must set different filename to save the subproblem 
    freepeak_param["dataFile1"] = dirname + "/" + funname + ".txt";
    freepeak->inputParameters().input(freepeak_param);
    freepeak->initialize(env.get());
    std::shared_ptr<ofec::Random> pro_rnd(new Random(0.5));
    freepeak->setRandom(pro_rnd);
    // number of variable , number of objective, number of constraint
    // multi party multi modal benchmark , for each freepeak , there is only one objective without contraints
    freepeak->setSizes(numDim, numObj, numCon);

    // 
  // Set KD-tree subspace structure:
  // The argument is a vector of pairs: each pair is <tree_name, vector of (subtree_name, space_size)>.
  // 
  // 
// Example (corresponding to the code block below):
// - The `root` node contains two child subtrees: `subtree1` and `subtree2`, each assigned a weight of 0.5.
// - `subtree1` contains two subspaces: `<funname>_s1_1` and `<funname>_s1_2`, each with a weight of 0.5.
// - `subtree2` contains two subspaces: `<funname>_s2_1` (weight 0.3) and `<funname>_s2_2` (weight 0.1).

    freepeak->setKDtree({

        {"root", { { "subtree1", 0.5}, {"subtree2", 0.5} } },
        {"subtree1", { {funname + "_s1_1", 0.1}, {funname + "_s1_2", 0.4} } },
        {"subtree2", { {funname + "_s2_1", 0.3}, {funname + "_s2_2", 0.2} } } });

    // design each subproblem in each subspace 
    for (int idsubtree(1); idsubtree <= 2; ++idsubtree)
    {
        for (int idpeak(1); idpeak <= 2; ++idpeak)
        {
            std::string subspace_name = funname + "_s" + std::to_string(idsubtree) + "_" + std::to_string(idpeak);
            ParameterMap subpro_param;
            subpro_param["subspace"] = subspace_name;
            subpro_param["generation_type"] = std::string("assigned");
            // must set different filename to save the subproblem 
            subpro_param["dataFile1"] = dirname + "/" + subspace_name + ".txt";


            auto subpro(Subproblem::create());
            // refer to the function: addInputParameters of the class to see the parameters and their range
            subpro->initialize(subpro_param, freepeak);
            // generate a distance and send the the distance parameter  to subpro;
            // each subpro has one distance .
            {
                std::cout << "These are the register name of the distance class\t";
                for (auto& it : FactoryFP<DistanceBase>::registerNames()) {
                    std::cout << it << "\t";
                }
                std::cout << std::endl;

                auto dis(FactoryFP<DistanceBase>::produce("Euclidean"));

                // there are no parameters to set for Euclidean distance, but the following step can not be ingored , because the register name will be saved.
                ParameterMap dis_param;
                dis->initialize(freepeak, subspace_name, dis_param);

                // save the parameters and send to the parameter to subpro
                //dis->recordInputParameters();
                subpro->setDistance(dis);

            }

            // each subpro has one function.
            {

                // there are no parameters to set for Euclidean distance, but the following step can not be ingored ,
                // because the register name will be saved.
                ParameterMap fun_param;
                fun_param["generation_type"] = std::string("assigned");
                fun_param["dataFile1"] = dirname + "/" + subspace_name + "_onepeak" + ".txt";

                // there are only one function of FunctionBase class, one_peak
                auto func(FactoryFP<FunctionBase>::produce("one_peak"));
                //func->inputParameters().input(fun_param);
                func->initialize(freepeak, subspace_name, fun_param);
                auto onepeak_func = dynamic_cast<ofec::free_peaks::OnePeakFunction*>(func);

                // each function has one or more onepeaks
                // for multiparty multi modal benchmark, each function has only one peak, 
                // one peak represent one objective function in each subspace
                // generate Onepeak 
                {
                    auto onepeak(FactoryFP<OnePeakBase>::produceRandom(rnd->uniform.next()));
                    ParameterMap onepeak_param;
                    onepeak_param["center_type"] = std::string("assigned");
                    onepeak_param["height"] = rnd->uniform.nextNonStd<Real>(0, 100); // Example height  
                    // in the range [-100, 100]^D
                    std::vector<Real> peak_position(numDim); // Example peak position
                    for (int idx(0); idx < numDim; ++idx) {
                        peak_position[idx] = rnd->uniform.nextNonStd<Real>(OnePeakBase::ms_x_range.first, OnePeakBase::ms_x_range.second);
                    }
                    onepeak_param["center_postion"] = peak_position;
                    onepeak->initialize(freepeak, subspace_name, onepeak_param);
                    //             onepeak->recordInputParameters();
                    onepeak_func->addOnePeaks(onepeak);

                }

                //       func->recordInputParameters();
                subpro->setFunction(func);

            }

            // there are no constraints 
            // subpro->addConstraint(param_con)

            // tranform in X space
            // there can be many transform methods, here we random select two 
            {
                //for (int idx(0); idx < 1; ++idx)
                {
                    ofec::ParameterMap default_param;
                    auto it(FactoryFP<X_TransformBase>::produceRandom(rnd->uniform.next()));

                    it->initialize(freepeak, subspace_name, default_param);
                    // refer to the   function: addInputParameters of the class to see the parameters and their range
                    // use its default parameters   
                    subpro->addVariableTransform(it);
                }
            }

            // transform in Y space
            // subpro->addObjectiveTransform(param_obj)
            {
                ofec::ParameterMap default_param;
                auto it(FactoryFP<Y_TransformBase>::produce("map_objective"));

                it->initialize(freepeak, subspace_name, default_param);
                // refer to the   function: addInputParameters of the class to see the parameters and their range
                // use its default parameters   
                subpro->addObjectiveTransform(it);
            }


            //subpro->recordInputParameters();
            freepeak->setSubproblem(subspace_name, subpro);





        }
    }

    freepeak->bindData();
    // after handmade setting each subproblem, we need to initialize the freepeak again to bind the data
   // freepeak->initialize(env.get());

    //printBenchmarkStructure(freepeak);





//      std::string outputMeshDir = "originMesh.txt";
      //outputProMesh(outputMeshDir, freepeak.get(), 400);
      // test output freepeak


      // output parameters

                  // generate file directory of subproblem

    std::cout << "freepeak path\t" << FreePeaks::directory() + dirname << std::endl;
    std::filesystem::create_directories(FreePeaks::directory() + dirname);
    std::cout << "Subproblem path\t" << Subproblem::directory() + dirname << std::endl;
    std::filesystem::create_directories(Subproblem::directory() + dirname);
    std::cout << "OnePeakFunction path\t" << OnePeakFunction::directory() + dirname << std::endl;
    std::filesystem::create_directories(OnePeakFunction::directory() + dirname);


    // set the generation type to read file to save the parameters  
    freepeak->inputParameters().at("generation_type")->setValue("read_file");
    for (auto& it : freepeak->subspaceTree().name_box_subproblem) {
        if (it.second.second != nullptr) {
            it.second.second->inputParameters().at("generation_type")->setValue("read_file");
            // it.second.second->recordInputParameters();
            it.second.second->function()->inputParameters().at("generation_type")->setValue("read_file");

            //         auto onepeak = CAST_OnePeakFunction(it.second.second->function());
            //         if(onepeak !=nullptr)
            //         {
            //             for (auto& peak : onepeak->onepeaksFun()) {
            //                 peak->inputParameters().at("center_type")->setValue("read_file");
            //                // peak->recordInputParameters();
            //             }
                     //}   
        }
    }

    freepeak->recordInputParameters();
    freepeak->outputTotalFile();


    {       // generate the same freepeak instance by read file
        std::string freepeakName = "free_peaks";

        ParameterMap freepeak_param;
        freepeak_param["generation_type"] = std::string("read_file");
        // must set different filename to save the subproblem 
        freepeak_param["dataFile1"] = dirname + "/" + funname + ".txt";

        std::shared_ptr<Environment> env2(generateEnvironmentByFactory(freepeakName));
        env2->recordInputParameters();
        env2->initialize();

        env2->setProblem(generateProblemByFactory(freepeakName));

        // the parameter of the handmade freepeak
        env2->problem()->inputParameters().input(freepeak_param);
        env2->problem()->recordInputParameters();
        env2->initializeProblem(0.5);



        outputProMesh("freepeak_mesh", env2.get(), 400);


    }



}


namespace ofec {
    void run(int argc, char* argv[]) {
        using namespace ofec;
        using namespace std;


        // Multiparty Multimodal Optimization (MPMMO)
        // Find ALL x* where consensus Q(x*) = min_p F_p(x*) is maximised.
        ofec::g_working_directory = "D:/Research/ThesisProject/freepeaks/Data/ofec_data_new";


        // Run the full MPMMO 2-D benchmark suite
        ofec::runMPMMCoEAExperiment();

    }
}

// outputFreePeaksMultipartyLandscape12
//   Generates high-quality 2D landscape files for all 12 benchmark problems.
void outputFreePeaksMultipartyLandscape12(const std::string& dir) {
    using namespace ofec;
    registerInstance();
    const size_t dimension = 2;
    const int grid_resolution = 400;
    const std::filesystem::path output_root = dir +
        "Visualization/free_peaks_multiparty/landscapes_12/D2";
    std::filesystem::create_directories(output_root);
    for (int suite_id = 1; suite_id <= 12; ++suite_id) {
        std::shared_ptr<Environment> env(generateEnvironmentByFactory("free_peaks_multiparty"));
        env->setProblem(generateProblemByFactory("free_peaks_multiparty"));
        ParameterMap pm;
        pm["suite_id"] = suite_id;
        pm["problem_dimension"] = dimension;
        env->problem()->inputParameters().input(pm);
        const Real instance_seed = static_cast<Real>(suite_id * 1000 + static_cast<int>(dimension));
        env->initializeProblem(instance_seed);
        auto* problem = dynamic_cast<FreePeaksMultiParty*>(env->problem());
        if (!problem) throw std::runtime_error("free_peaks_multiparty initialization failed");
        const auto& spec = problem->currentSpec();
        std::filesystem::path suite_dir = output_root /
            (suitePLabel(suite_id) + "_suite" + std::to_string(suite_id) + "_" + spec.name);
        std::filesystem::create_directories(suite_dir);

        // 400x400 grid (j=x1 outer, i=x0 inner)
        std::vector<std::shared_ptr<SolutionBase>> sols;
        std::vector<std::vector<double>> indiObjs;
        sols.reserve(static_cast<size_t>(grid_resolution) * grid_resolution);
        indiObjs.reserve(static_cast<size_t>(grid_resolution) * grid_resolution);
        for (int j = 0; j < grid_resolution; ++j) {
            for (int i = 0; i < grid_resolution; ++i) {
                const Real x0 = static_cast<Real>(i) / static_cast<Real>(grid_resolution - 1);
                const Real x1 = static_cast<Real>(j) / static_cast<Real>(grid_resolution - 1);
                std::shared_ptr<SolutionBase> sol(problem->createSolution());
                auto& var = dynamic_cast<Solution<VariableVector<Real>>&>(*sol).variable().vector();
                var[0] = x0;
                var[1] = x1;
                sol->evaluate(env.get(), false);
                indiObjs.emplace_back(sol->objective());
                sols.emplace_back(sol);
            }
        }


        const int N = static_cast<int>(indiObjs.size());
        std::vector<double> consensus(N);
        for (int k = 0; k < N; ++k)
            consensus[k] = std::min(indiObjs[k][0], indiObjs[k][1]);

        // Write grid_2d.tsv
        // Column consensus_min = min(F0,F1): the MPMMO quality measure.
        // Shared optima are grid cells where consensus_min ≈ 90.
        std::ofstream grid_f(suite_dir / "grid_2d.tsv");
        grid_f << std::fixed << std::setprecision(10);
        grid_f << "i\tj\tx0_norm\tx1_norm\tx0_plot\tx1_plot\tf_party0\tf_party1\tconsensus_min\n";
        for (int j = 0; j < grid_resolution; ++j) {
            for (int i = 0; i < grid_resolution; ++i) {
                const int idx = j * grid_resolution + i;
                const auto& var = dynamic_cast<const Solution<VariableVector<Real>>&>(*sols[idx]).variable().vector();
                grid_f << i << '\t' << j << '\t'
                    << var[0] << '\t' << var[1] << '\t'
                    << (var[0] * 7.0 - 3.5) << '\t'
                    << (var[1] * 7.0 - 3.5) << '\t'
                    << sols[idx]->objective(0) << '\t'
                    << sols[idx]->objective(1) << '\t'
                    << consensus[idx] << "\n";
            }
        }

        // Write optima_2d.tsv — shared optima where Q(x*) = min(F0,F1) ≈ 90
        std::ofstream optima_f(suite_dir / "optima_2d.tsv");
        optima_f << std::fixed << std::setprecision(10);
        optima_f << "idx\tx0_norm\tx1_norm\tx0_plot\tx1_plot\tf_party0\tf_party1\tconsensus_min\n";
        if (problem->optima()) {
            for (size_t opt_id = 0; opt_id < problem->optima()->numberSolutions(); ++opt_id) {
                const auto& opt_sol = problem->optima()->solution(opt_id);
                const auto& opt_var = opt_sol.variable().vector();
                const double q = std::min(opt_sol.objective(0), opt_sol.objective(1));
                optima_f << opt_id << '\t'
                    << opt_var[0] << '\t' << opt_var[1] << '\t'
                    << (opt_var[0] * 7.0 - 3.5) << '\t'
                    << (opt_var[1] * 7.0 - 3.5) << '\t'
                    << opt_sol.objective(0) << '\t'
                    << opt_sol.objective(1) << '\t'
                    << q << "\n";
            }
        }

        // Write run_info.txt
        std::ofstream info_f(suite_dir / "run_info.txt");
        info_f << "p_id " << suitePLabel(suite_id) << "\n";
        info_f << "suite " << suite_id << "\n";
        info_f << "name " << spec.name << "\n";
        info_f << "feature " << spec.feature << "\n";
        info_f << "dimension " << dimension << "\n";
        info_f << "instance_seed " << static_cast<int>(instance_seed) << "\n";
        info_f << "grid_resolution " << grid_resolution << "\n";
        info_f << "consensus_definition min(F_party0, F_party1)\n";
        info_f << "shared_optima_count " << (problem->optima() ? problem->optima()->numberSolutions() : 0) << "\n";
        std::cout << "[landscape] " << suitePLabel(suite_id)
            << " suite " << suite_id << " " << spec.name
            << " -> " << suite_dir.string() << "\n";
    }
}

#endif // !OFEC_CUSTOM_METHOD_HPP
