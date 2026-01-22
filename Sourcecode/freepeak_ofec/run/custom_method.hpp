/*#ifndef OFEC_CUSTOM_METHOD_HPP
#define OFEC_CUSTOM_METHOD_HPP

#include "../core/global.h"
#include "interface.h"




#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

void  copyFile(std::string &ofec_dir ,
std::string &file_list, // 包含list(APPEND...)的文件
std::string &target_dir // 目标路径
){
    //if (argc < 4) {
    //    std::cerr << "用法: " << argv[0] << " <OFEC_DIR> <file_list> <target_dir>\n";
    //    return 1;
    //}



    std::ifstream fin(file_list);
    if (!fin) {
        std::cout << "无法打开文件: " << file_list << "\n";
        return;
    }

    std::string line;
    while (std::getline(fin, line)) {
        // 找到 ${OFEC_DIR}
        size_t pos = line.find("${OFEC_DIR}");
        if (pos == std::string::npos) continue;  // 没有路径的行跳过

        // 去掉前后空白和括号
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        // 去掉结尾的 ')'
        if (!line.empty() && (line.back() == ')' || line.back() == ',')) {
            line.pop_back();
        }

        // 替换成实际路径
        std::string relative_path = line.substr(pos + std::string("${OFEC_DIR}").size());
        fs::path src = fs::path(ofec_dir) / relative_path.substr(0); // 去掉 '/'
        fs::path dst = fs::path(target_dir) / relative_path.substr(0);

        try {
            fs::create_directories(dst.parent_path());
            fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
            std::cout << "复制: " << src << " -> " << dst << "\n";
        }
        catch (const std::exception& e) {
            std::cout << "复制失败: " << src << " 错误: " << e.what() << "\n";
        }
    }

    //return 0;
}




namespace ofec {
    void run(int argc, char* argv[]) {
        using namespace ofec;
        using namespace std;
        ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
        ofec::g_working_directory = "//172.29.41.69/share/2018/diaoyiya/ofec-data";
        ofec::g_working_directory = "//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data";
        ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
        ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
        //ofec::g_working_directory = "//172.24.24.151/e/DiaoYiya/code/data/ofec-data";
        ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/";
        ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
        ofec::g_working_directory = "D:/code/Data/ofec_data/";
        //ofec::g_working_directory = "/home/dyy/data/ofec_data/ofec_data/";
        //ofec::g_working_directory = "/data/Share/Student/2018/diaoyiya/data/ofec_data/";
      //  ofec::g_working_directory = "/mnt/dataShare/Student/2018/diaoyiya/data/ofec_data/";

        registerInstance();


        std::string ofec_dir = OFEC_DIR;
        std::string savedir = "D:/code/multi_party_multi_modal/freepeak_ofec";
        std::string list_file = ofec_dir + std::string("/SrcFiles.cmake");
        copyFile(ofec_dir, list_file,savedir);
    }
}

#endif // !OFEC_CUSTOM_METHOD_HPP*/

/*#ifndef OFEC_CUSTOM_METHOD_HPP
#define OFEC_CUSTOM_METHOD_HPP

#include "../core/global.h"
#include "interface.h"
#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../core/problem/solution.h"
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include <iomanip>

namespace ofec {

    // PROBLEM WRAPPER (Low-level Interface to FreePeaks)
    class FreePeaksManual : public FreePeaks {
    public:
        void configure(const std::string& filename) {
            this->m_file_name = filename;
            this->m_generation_type = "read_file";
            this->m_number_objectives = 1;
        }
        void initialize_(Environment* env) override {
            std::string full_path = g_working_directory + "/instance/problem/continuous/free_peaks/" + m_file_name;
            FreePeaks::initialize_(env);
        }
        // Helper to evaluate a full D-dimensional solution
        void evaluateSolution(ofec::Solution<>& sol) {
            std::vector<ofec::Real> objs(this->numberObjectives());
            std::vector<ofec::Real> cons(this->numberConstraints());
            this->evaluate(sol.variable(), objs, cons);
            for (size_t i = 0; i < objs.size(); i++) sol.objective(i) = objs[i];
        }
    };

    // MPMMO BENCHMARK DEFINITION
    //    This class transforms the standard FreePeaks into a Multiparty Problem
    class MPMMO_Benchmark {
    private:
        std::shared_ptr<FreePeaksManual> m_problem;
        int m_num_parties;
        int m_total_vars;

        // Stores which variables belong to which party
        // e.g., party_vars[0] = {0, 1, 2}, party_vars[1] = {3, 4}
        std::vector<std::vector<int>> m_party_indices;

    public:
        MPMMO_Benchmark(std::shared_ptr<FreePeaksManual> prob, int num_parties)
            : m_problem(prob), m_num_parties(num_parties) {

            m_total_vars = m_problem->numberVariables();
            decomposeVariables();
        }

        // Logic to split the decision variables among parties
        void decomposeVariables() {
            m_party_indices.resize(m_num_parties);
            int vars_per_party = m_total_vars / m_num_parties;
            int remainder = m_total_vars % m_num_parties;

            int current_idx = 0;
            for (int i = 0; i < m_num_parties; i++) {
                int count = vars_per_party + (i < remainder ? 1 : 0);
                for (int k = 0; k < count; k++) {
                    m_party_indices[i].push_back(current_idx++);
                }
            }
        }

        // --- THE CORE MPMMO EVALUATION ---
        // Evaluation requires:
        // The Party ID doing the evaluation
        // That Party's partial solution (subset of variables)
        // Context Vectors from all OTHER parties to fill the gaps
        double evaluateParty(int party_id, const std::vector<double>& party_vars, const std::vector<std::vector<double>>& context_vectors) {

            // Create a full solution container
            ofec::Solution<> full_sol(m_problem->numberObjectives(), m_problem->numberConstraints());
            full_sol.variable().resize(m_total_vars);

            // Fill in the Active Party's variables
            const auto& my_indices = m_party_indices[party_id];
            for (size_t i = 0; i < my_indices.size(); i++) {
                full_sol.variable()[my_indices[i]] = party_vars[i];
            }

            // Fill in the Context (Collaborators) from other parties
            for (int other_pid = 0; other_pid < m_num_parties; other_pid++) {
                if (other_pid == party_id) continue;

                const auto& other_indices = m_party_indices[other_pid];
                // Safety check
                if (context_vectors[other_pid].size() != other_indices.size()) {
                    std::cerr << "Error: Context vector size mismatch for party " << other_pid << std::endl;
                    return -1e9;
                }

                for (size_t i = 0; i < other_indices.size(); i++) {
                    full_sol.variable()[other_indices[i]] = context_vectors[other_pid][i];
                }
            }

            // Evaluate on the actual FreePeaks landscape
            m_problem->evaluateSolution(full_sol);
            return full_sol.objective(0);
        }

        // Getters
        const std::vector<int>& getPartyVars(int party_id) const { return m_party_indices[party_id]; }
        int getNumParties() const { return m_num_parties; }
    };

    // SIMPLE COOPERATIVE SOLVER (To Test the Benchmark)
    class Solver_MPM_Basic {
    private:
        std::shared_ptr<MPMMO_Benchmark> m_benchmark;
        int m_pop_size;

        // Population structure: [PartyID][IndividualID][VariableIndex]
        std::vector<std::vector<std::vector<double>>> m_populations;
        // Fitness values: [PartyID][IndividualID]
        std::vector<std::vector<double>> m_fitness;

        // "Representatives" - The best solution found by each party so far
        std::vector<std::vector<double>> m_representatives;

    public:
        Solver_MPM_Basic(std::shared_ptr<MPMMO_Benchmark> bench, int pop_size)
            : m_benchmark(bench), m_pop_size(pop_size) {
        }

        void initialize() {
            int num_parties = m_benchmark->getNumParties();
            m_populations.resize(num_parties);
            m_fitness.resize(num_parties);
            m_representatives.resize(num_parties);

            std::cout << ">>> [MPMMO] Initializing " << num_parties << " parties..." << std::endl;

            for (int p = 0; p < num_parties; p++) {
                int num_vars = m_benchmark->getPartyVars(p).size();
                m_populations[p].resize(m_pop_size, std::vector<double>(num_vars));
                m_fitness[p].resize(m_pop_size);

                // Random Init
                for (int i = 0; i < m_pop_size; i++) {
                    for (int j = 0; j < num_vars; j++) {
                        m_populations[p][i][j] = (double)rand() / RAND_MAX;
                    }
                }
                // Initial Representative (Just pick the first one randomly for now)
                m_representatives[p] = m_populations[p][0];
            }
        }

        void run(int generations) {
            int num_parties = m_benchmark->getNumParties();
            std::cout << ">>> [MPMMO] Run Co-Evolution Loop..." << std::endl;

            for (int gen = 0; gen < generations; gen++) {

                // --- Round Robin Cooperative Co-Evolution ---
                for (int active_party = 0; active_party < num_parties; active_party++) {

                    double best_val_in_party = -1e9;
                    int best_idx_in_party = 0;

                    // Evaluate everyone in this party
                    for (int i = 0; i < m_pop_size; i++) {
                        // Evaluate using Representatives of other parties as context
                        double fit = m_benchmark->evaluateParty(
                            active_party,
                            m_populations[active_party][i],
                            m_representatives
                        );
                        m_fitness[active_party][i] = fit;

                        if (fit > best_val_in_party) {
                            best_val_in_party = fit;
                            best_idx_in_party = i;
                        }
                    }

                    // UPDATE REPRESENTATIVE
                    // The best individual of this party becomes the new context for others
                    m_representatives[active_party] = m_populations[active_party][best_idx_in_party];

                    // Simple "Evolution" (Mutation only for baseline)
                    // If you want full DE/PSO, it goes here.
                    for (int i = 0; i < m_pop_size; i++) {
                        if (i == best_idx_in_party) continue; // Elitism
                        int vars_count = m_populations[active_party][i].size();
                        for (int v = 0; v < vars_count; v++) {
                            // Add small random noise
                            double noise = ((double)rand() / RAND_MAX - 0.5) * 0.1;
                            m_populations[active_party][i][v] += noise;
                            // Clamp
                            if (m_populations[active_party][i][v] > 1.0) m_populations[active_party][i][v] = 1.0;
                            if (m_populations[active_party][i][v] < 0.0) m_populations[active_party][i][v] = 0.0;
                        }
                    }
                }

                // Output global progress
                if ((gen + 1) % 10 == 0) {
                    // Evaluate the global composed solution
                    // We treat Party 0 as the 'evaluator' just to call the function, using its rep + others reps
                    double global_fitness = m_benchmark->evaluateParty(0, m_representatives[0], m_representatives);

                    std::cout << "Gen " << std::setw(3) << gen + 1 << " | Global Fitness (Cooperated): " << global_fitness << " | Representatives: ";
                    for (int p = 0; p < num_parties; p++) std::cout << "[P" << p << " Updated] ";
                    std::cout << std::endl;
                }
            }
        }
    };

    // MAIN EXECUTION
    void run(int argc, char* argv[]) {
        using namespace ofec;
        using namespace std;

        // Path Setup
        ofec::g_working_directory = R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Data\ofec_data_new)";
        std::string my_data_file = "sop/2d_s1_to_s8.txt";

        // Create the Physics (FreePeaks)
        auto raw_problem = std::make_shared<FreePeaksManual>();
        raw_problem->configure(my_data_file);

        try {
            raw_problem->initialize(nullptr);

            // Create the Benchmark Structure (MPMMO Wrapper)
            // Splitting the problem into 2 Parties
            int num_parties = 2;
            auto mpmmo_benchmark = std::make_shared<MPMMO_Benchmark>(raw_problem, num_parties);

            cout << "------------------------------------------------" << endl;
            cout << "BENCHMARK INITIALIZED: MPMMO on FreePeaks" << endl;
            cout << "Total Variables: " << raw_problem->numberVariables() << endl;
            cout << "Number of Parties: " << num_parties << endl;
            cout << "------------------------------------------------" << endl;
            for (int i = 0; i < num_parties; i++) {
                cout << "Party " << i << " controls vars: ";
                for (int idx : mpmmo_benchmark->getPartyVars(i)) cout << idx << " ";
                cout << endl;
            }
            cout << "------------------------------------------------" << endl;

            // 3. Run the Solver
            Solver_MPM_Basic solver(mpmmo_benchmark, 50); // 50 pop size per party
            solver.initialize();
            solver.run(50); // Run for 50 generations

            cout << "------------------------------------------------" << endl;
            cout << "MPMMO EXPERIMENT COMPLETE." << endl;
            cout << "------------------------------------------------" << endl;
        }
        catch (const std::exception& e) {
            cerr << "Error: " << e.what() << endl;
        }
        system("pause");
    }
}

#endif // !OFEC_CUSTOM_METHOD_HPP*/

/*#ifndef OFEC_CUSTOM_METHOD_HPP
#define OFEC_CUSTOM_METHOD_HPP

#include "../core/global.h"
#include "interface.h"
#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../core/problem/solution.h"
#include <iostream>
#include <vector>
#include <numeric>

namespace ofec {


    // THE BASIC PROBLEM (FreePeaks)

    // This allows to load the map manually without Factory errors
    class FreePeaksManual : public FreePeaks {
    public:
        void configure(const std::string& filename) {
            this->m_file_name = filename;
            this->m_generation_type = "read_file";
            this->m_number_objectives = 1;
        }
        void initialize_(Environment* env) override {
            std::string full_path = g_working_directory + "/instance/problem/continuous/free_peaks/" + m_file_name;
            FreePeaks::initialize_(env);
        }
        // Helper to expose the protected evaluate function
        void evaluateSolution(ofec::Solution<>& sol) {
            std::vector<ofec::Real> objs(this->numberObjectives());
            std::vector<ofec::Real> cons(this->numberConstraints());
            this->evaluate(sol.variable(), objs, cons);
            for (size_t i = 0; i < objs.size(); i++) sol.objective(i) = objs[i];
        }
    };


    // MPMMO BENCHMARK WRAPPER
    class MPMMO_Benchmark {
    private:
        std::shared_ptr<FreePeaksManual> m_problem;
        int m_num_parties;
        std::vector<std::vector<int>> m_party_indices; // Which vars belong to which party

    public:
        MPMMO_Benchmark(std::shared_ptr<FreePeaksManual> prob, int num_parties)
            : m_problem(prob), m_num_parties(num_parties) {

            decomposeVariables();
        }

        // Split the variables (10 variables mean Party A gets 5, Party B gets 5)
        void decomposeVariables() {
            int total_vars = m_problem->numberVariables();
            m_party_indices.resize(m_num_parties);

            int current_var = 0;
            // Simple equal distribution
            for (int p = 0; p < m_num_parties; p++) {
                int count = total_vars / m_num_parties;
                // Handle remainders (if 5 vars, 2 parties get 3 and 2)
                if (p < total_vars % m_num_parties) count++;

                for (int k = 0; k < count; k++) {
                    m_party_indices[p].push_back(current_var++);
                }
            }
        }

        // Contextual Evaluation
        // Party 'p' wants to test their 'vars', but needs 'context' from others to calculate fitness.
        double evaluateParty(int party_id, const std::vector<double>& my_vars, const std::vector<std::vector<double>>& context) {

            // Create a full solution container
            ofec::Solution<> full_sol(m_problem->numberObjectives(), m_problem->numberConstraints());
            full_sol.variable().resize(m_problem->numberVariables());

            // Fill in MY variables
            const auto& my_indices = m_party_indices[party_id];
            for (size_t i = 0; i < my_indices.size(); i++) {
                full_sol.variable()[my_indices[i]] = my_vars[i];
            }

            // Fill in OTHER parties' variables (Context)
            for (int other = 0; other < m_num_parties; other++) {
                if (other == party_id) continue;

                const auto& other_indices = m_party_indices[other];
                for (size_t i = 0; i < other_indices.size(); i++) {
                    full_sol.variable()[other_indices[i]] = context[other][i];
                }
            }

            // Ask FreePeaks for the score
            m_problem->evaluateSolution(full_sol);
            return full_sol.objective(0);
        }

        const std::vector<int>& getPartyIndices(int p) { return m_party_indices[p]; }
    };

    // MAIN
    void run(int argc, char* argv[]) {
        using namespace ofec;
        using namespace std;

        ofec::g_working_directory = R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Data\ofec_data_new)";
        std::string my_data_file = "sop/2d_s1_to_s8.txt";

        cout << ">>> 1. Loading FreePeaks Map..." << endl;
        auto core_problem = std::make_shared<FreePeaksManual>();
        core_problem->configure(my_data_file);

        try {
            core_problem->initialize(nullptr);

            // Creating MPMMO ENVIRONMENT
            int num_parties = 2; // Biparty Optimization
            cout << ">>> 2. Splitting into " << num_parties << " Parties (MPMMO Setup)..." << endl;

            MPMMO_Benchmark benchmark(core_problem, num_parties);

            // the split
            for (int p = 0; p < num_parties; p++) {
                cout << "    Party " << p << " controls Variables: [ ";
                for (int idx : benchmark.getPartyIndices(p)) cout << idx << " ";
                cout << "]" << endl;
            }

            // Simulating a cooperation
            cout << ">>> 3. Testing Cooperation..." << endl;

            // Every party starts at 0.5
            vector<vector<double>> context_vectors(num_parties);
            for (int p = 0; p < num_parties; p++) {
                context_vectors[p].resize(benchmark.getPartyIndices(p).size(), 0.5);
            }

            // Party 0 tries to improve its part
            cout << "    Party 0 is evaluating a new solution..." << endl;
            vector<double> party0_test = { 0.8 }; // Imagine Party 0 tries 0.8
            if (party0_test.size() != benchmark.getPartyIndices(0).size()) party0_test.resize(benchmark.getPartyIndices(0).size(), 0.8);

            double fit = benchmark.evaluateParty(0, party0_test, context_vectors);

            cout << "    Resulting Fitness: " << fit << endl;
            cout << ">>> MPMMO Benchmark Ready." << endl;

        }
        catch (const std::exception& e) {
            cerr << "Error: " << e.what() << endl;
        }
        system("pause");
    }
}

#endif // !OFEC_CUSTOM_METHOD_HPP*/

#ifndef OFEC_CUSTOM_METHOD_HPP
#define OFEC_CUSTOM_METHOD_HPP


#include "../core/global.h"
#include "interface.h"
#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../core/problem/solution.h"
#include <iostream>
#include <vector>
#include <numeric>   // for std::iota
#include <algorithm> // for std::max_element

namespace ofec {

    class FreePeaksManual : public FreePeaks {
    public:
        void configure(const std::string& filename) {
            this->m_file_name = filename;
            this->m_generation_type = "read_file";
            this->m_number_objectives = 1;
        }

        // Override to ensure the file path is correct
        void initialize_(Environment* env) override {
            std::string full_path = g_working_directory + "/instance/problem/continuous/free_peaks/" + m_file_name;
            FreePeaks::initialize_(env);
        }

        // Helper to run the evaluation logic
        void evaluateSolution(ofec::Solution<>& sol) {
            std::vector<ofec::Real> objs(this->numberObjectives());
            std::vector<ofec::Real> cons(this->numberConstraints());
            this->evaluate(sol.variable(), objs, cons);
            for (size_t i = 0; i < objs.size(); i++) sol.objective(i) = objs[i];
        }
    };

    // MPMMO_Benchmark 
    class MPMMO_Benchmark {
    private:
        std::shared_ptr<FreePeaksManual> m_problem;
        int m_num_parties;
        std::vector<std::vector<int>> m_party_indices; // Stores which variables belong to whom

    public:
        MPMMO_Benchmark(std::shared_ptr<FreePeaksManual> prob, int num_parties)
            : m_problem(prob), m_num_parties(num_parties) {
            decomposeVariables();
        }

        // Decomposition of decision variables
        void decomposeVariables() {
            int total_vars = m_problem->numberVariables();
            m_party_indices.resize(m_num_parties);

            int current_var = 0;
            for (int p = 0; p < m_num_parties; p++) {
                // Distribute variables evenly
                int count = total_vars / m_num_parties;
                if (p < total_vars % m_num_parties) count++; // Handle remainders

                for (int k = 0; k < count; k++) {
                    m_party_indices[p].push_back(current_var++);
                }
            }
        }

        // Context Vector Design
        // To evaluate Party P, need representatives from all other parties.
        double evaluateParty(int party_id, const std::vector<double>& my_vars, const std::vector<std::vector<double>>& context) {

            // Creating a full-length solution
            ofec::Solution<> full_sol(m_problem->numberObjectives(), m_problem->numberConstraints());
            full_sol.variable().resize(m_problem->numberVariables());

            // Inserting the Active Party's variables
            const auto& my_indices = m_party_indices[party_id];
            for (size_t i = 0; i < my_indices.size(); i++) {
                full_sol.variable()[my_indices[i]] = my_vars[i];
            }

            // Inserting Context (Collaborators)
            for (int other = 0; other < m_num_parties; other++) {
                if (other == party_id) continue;
                const auto& other_indices = m_party_indices[other];
                for (size_t i = 0; i < other_indices.size(); i++) {
                    full_sol.variable()[other_indices[i]] = context[other][i];
                }
            }

            // Calculating Fitness
            m_problem->evaluateSolution(full_sol);
            return full_sol.objective(0);
        }

        // Getters
        const std::vector<int>& getPartyIndices(int p) { return m_party_indices[p]; }
        int getNumParties() { return m_num_parties; }
    };

    // The Solver Loop
    class MPMMO_Solver {
        std::shared_ptr<MPMMO_Benchmark> m_bench;
        int m_pop_size;
        // Population: [PartyID][IndividualID][VariableValues]
        std::vector<std::vector<std::vector<double>>> m_populations;
        // Fitness: [PartyID][IndividualID]
        std::vector<std::vector<double>> m_fitness_vals;
        // Representatives: [PartyID][BestVariables]
        std::vector<std::vector<double>> m_representatives;

    public:
        MPMMO_Solver(std::shared_ptr<MPMMO_Benchmark> b, int pop) : m_bench(b), m_pop_size(pop) {}

        void initialize() {
            int n_parties = m_bench->getNumParties();
            m_populations.resize(n_parties);
            m_fitness_vals.resize(n_parties);
            m_representatives.resize(n_parties);

            std::cout << ">>> [Solver] Initializing Populations for " << n_parties << " Parties..." << std::endl;

            for (int p = 0; p < n_parties; p++) {
                int n_vars = m_bench->getPartyIndices(p).size();
                m_populations[p].resize(m_pop_size, std::vector<double>(n_vars));
                m_fitness_vals[p].resize(m_pop_size);

                // Random Initialization
                for (int i = 0; i < m_pop_size; i++) {
                    for (int v = 0; v < n_vars; v++) {
                        m_populations[p][i][v] = (double)rand() / RAND_MAX;
                    }
                }
                // Pick initial random representative
                m_representatives[p] = m_populations[p][0];
            }
        }

        void run_one_generation(int gen) {
            int n_parties = m_bench->getNumParties();

            // Loop through each party (Round Robin)
            for (int active_p = 0; active_p < n_parties; active_p++) {

                double best_val = -1e9;
                int best_idx = 0;

                // Evaluate every individual in this party
                for (int i = 0; i < m_pop_size; i++) {
                    // Asking Benchmark to evaluate using CURRENT representatives as context
                    double fit = m_bench->evaluateParty(active_p, m_populations[active_p][i], m_representatives);
                    m_fitness_vals[active_p][i] = fit;

                    if (fit > best_val) {
                        best_val = fit;
                        best_idx = i;
                    }
                }

                // Updating Representative for this party
                // In full MPMMO, will update multiple niche representatives here
                m_representatives[active_p] = m_populations[active_p][best_idx];

                // Simple Mutation (Evolution)
                for (int i = 0; i < m_pop_size; i++) {
                    if (i == best_idx) continue; // Keep best
                    for (size_t v = 0; v < m_populations[active_p][i].size(); v++) {
                        // Add noise
                        m_populations[active_p][i][v] += ((double)rand() / RAND_MAX - 0.5) * 0.1;
                    }
                }
            }
        }

        double getGlobalBest() {
            // Check global fitness of the combined representatives
            return m_bench->evaluateParty(0, m_representatives[0], m_representatives);
        }
    };

    // MAIN FUNCTION
    void run(int argc, char* argv[]) {
        using namespace ofec;
        using namespace std;

        ofec::g_working_directory = R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Data\ofec_data_new)";
        std::string my_data_file = "sop/1d_1_s1_2_s1_h95_s6_h95.txt";

        // Loading FreePeaks
        auto core_problem = std::make_shared<FreePeaksManual>();
        core_problem->configure(my_data_file);

        try {
            core_problem->initialize(nullptr);

            // Create MPMMO Benchmark (2 Parties)
            MPMMO_Benchmark benchmark(core_problem, 2);

            // Create Solver and Run
            MPMMO_Solver solver(std::make_shared<MPMMO_Benchmark>(benchmark), 50); // 50 Individuals
            solver.initialize();

            cout << ">>> Starting MPMMO Optimization Loop..." << endl;
            for (int g = 0; g < 20; g++) {
                solver.run_one_generation(g);
                cout << "Gen " << g + 1 << " | Current Best Fitness: " << solver.getGlobalBest() << endl;
            }
            cout << ">>> Done." << endl;

        }
        catch (const std::exception& e) {
            cerr << "Error: " << e.what() << endl;
        }
        system("pause");
    }
}

#endif // !OFEC_CUSTOM_METHOD_HPP