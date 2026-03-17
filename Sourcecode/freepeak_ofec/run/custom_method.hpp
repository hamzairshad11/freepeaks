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
std::string &file_list, // °üº¬list(APPEND...)µÄÎÄ¼þ
std::string &target_dir // Ä¿±êÂ·¾¶
){
    //if (argc < 4) {
    //    std::cerr << "ÓÃ·¨: " << argv[0] << " <OFEC_DIR> <file_list> <target_dir>\n";
    //    return 1;
    //}



    std::ifstream fin(file_list);
    if (!fin) {
        std::cout << "ÎÞ·¨´ò¿ªÎÄ¼þ: " << file_list << "\n";
        return;
    }

    std::string line;
    while (std::getline(fin, line)) {
        // ÕÒµ½ ${OFEC_DIR}
        size_t pos = line.find("${OFEC_DIR}");
        if (pos == std::string::npos) continue;  // Ã»ÓÐÂ·¾¶µÄÐÐÌø¹ý

        // È¥µôÇ°ºó¿Õ°×ºÍÀ¨ºÅ
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        // È¥µô½áÎ²µÄ ')'
        if (!line.empty() && (line.back() == ')' || line.back() == ',')) {
            line.pop_back();
        }

        // Ìæ»»³ÉÊµ¼ÊÂ·¾¶
        std::string relative_path = line.substr(pos + std::string("${OFEC_DIR}").size());
        fs::path src = fs::path(ofec_dir) / relative_path.substr(0); // È¥µô '/'
        fs::path dst = fs::path(target_dir) / relative_path.substr(0);

        try {
            fs::create_directories(dst.parent_path());
            fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
            std::cout << "¸´ÖÆ: " << src << " -> " << dst << "\n";
        }
        catch (const std::exception& e) {
            std::cout << "¸´ÖÆÊ§°Ü: " << src << " ´íÎó: " << e.what() << "\n";
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

/*#ifndef OFEC_CUSTOM_METHOD_HPP
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
            //std::string full_path = g_working_directory + "/instance/problem/continuous/free_peaks/" + m_file_name;
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
        std::string my_data_file = "sop/2d_2_s1_1_1.txt";

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

            cout << "Starting MPMMO Optimization" << endl;
            for (int g = 0; g < 20; g++) {
                solver.run_one_generation(g);
                cout << "Gen " << g + 1 << " | Current Best Fitness: " << solver.getGlobalBest() << endl;
            }
            cout << "Done." << endl;

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
#include <numeric>   // for std::iota
#include <algorithm> // for std::max_element
#include <random>

namespace ofec {

    class FreePeaksManual : public FreePeaks {
    public:
        void configure(const std::string& filename) {
            this->m_file_name = filename;
            this->m_generation_type = "read_file";
            this->m_number_objectives = 1;
        }

        // Override to use base class initialization (correctly uses g_working_directory)
        void initialize_(Environment* env) override {
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
        int getTotalVariables() const { return m_problem->numberVariables(); }
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

        // NBC-based initialization 
        void initializeWithNBC() {
            int n_parties = m_bench->getNumParties();
            m_populations.resize(n_parties);
            m_fitness_vals.resize(n_parties);
            m_representatives.resize(n_parties);

            // Step 1: Create a large initial population for NBC
            const int NBC_POP_SIZE = 500;
            std::vector<std::vector<double>> global_pop(NBC_POP_SIZE);
            std::vector<double> global_fitness(NBC_POP_SIZE);

            // Generate random solutions
            for (int i = 0; i < NBC_POP_SIZE; i++) {
                global_pop[i].resize(m_bench->getTotalVariables());
                for (size_t v = 0; v < global_pop[i].size(); v++) {
                    global_pop[i][v] = (double)rand() / RAND_MAX;
                }
                // Evaluate using a neutral context (e.g., all 0.5)
                std::vector<std::vector<double>> neutral_context(n_parties);
                for (int p = 0; p < n_parties; p++) {
                    int n_vars = m_bench->getPartyIndices(p).size();
                    neutral_context[p] = std::vector<double>(n_vars, 0.5);
                }
                global_fitness[i] = m_bench->evaluateParty(0, global_pop[i], neutral_context);
            }

            // Step 2: Apply NBC clustering (simplified version), reference Lin et al. (2021)
            // Here, we just sort and split into niches
            std::vector<int> indices(NBC_POP_SIZE);
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](int a, int b) { return global_fitness[a] > global_fitness[b]; });

            // Distribute top individuals to each party's population
            for (int p = 0; p < n_parties; p++) {
                int n_vars = m_bench->getPartyIndices(p).size();
                m_populations[p].resize(m_pop_size, std::vector<double>(n_vars));
                m_fitness_vals[p].resize(m_pop_size);

                for (int i = 0; i < m_pop_size; i++) {
                    int src_idx = indices[(p * m_pop_size + i) % NBC_POP_SIZE];
                    const auto& my_indices = m_bench->getPartyIndices(p);
                    for (size_t v = 0; v < my_indices.size(); v++) {
                        m_populations[p][i][v] = global_pop[src_idx][my_indices[v]];
                    }
                }
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
                        // Ensure bounds [0, 1]
                        if (m_populations[active_p][i][v] < 0) m_populations[active_p][i][v] = 0;
                        if (m_populations[active_p][i][v] > 1) m_populations[active_p][i][v] = 1;
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

        // Path
        ofec::g_working_directory = R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\freepeak_ofec\freepeak_ofec)";
        std::string my_data_file = "sop/2d_2_s1_1_1.txt"; //  Using a 2D instance

        // Loading FreePeaks
        auto core_problem = std::make_shared<FreePeaksManual>();
        core_problem->configure(my_data_file);

        try {
            core_problem->initialize(nullptr);

            // Create MPMMO Benchmark (2 Parties)
            MPMMO_Benchmark benchmark(core_problem, 2);

            // Create Solver and Run
            MPMMO_Solver solver(std::make_shared<MPMMO_Benchmark>(benchmark), 50); // 50 Individuals

            //  Use NBC-based init as per thesis proposal
            solver.initializeWithNBC();

            cout << "Starting MPMMO Optimization with NBC Initialization" << endl;
            for (int g = 0; g < 20; g++) {
                solver.run_one_generation(g);
                cout << "Gen " << g + 1 << " | Current Best Fitness: " << solver.getGlobalBest() << endl;
            }
            cout << "Done." << endl;

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
#include <algorithm>
#include <cmath>
#include <queue>
#include <map>
#include <iomanip>  // For std::setprecision
#include <cstdlib>  // For rand(), srand()
#include <ctime>    // For time()

namespace ofec {

    // FREEPEAKS WRAPPER (to safely expose protected members via public interface)
    class FreePeaksWrapper : public FreePeaks {
    public:
        void configure(const std::string& filename) {
            this->m_file_name = filename;
            this->m_generation_type = "read_file";
            this->m_number_objectives = 1;
        }

        // Public wrapper for protected evaluate() - uses Solution<> interface
        void evaluateSolution(ofec::Solution<>& sol) {
            std::vector<ofec::Real> objs(this->numberObjectives());
            std::vector<ofec::Real> cons(this->numberConstraints());
            this->evaluate(sol.variable(), objs, cons);
            for (size_t i = 0; i < objs.size(); i++) {
                sol.objective(i) = objs[i];
            }
        }
    };

    // NEAREST-BETTER CLUSTERING (Graph-based)
    struct NBC_Edge {
        int from, to;
        double weight;
    };

    std::vector<std::vector<int>> nearestBetterClustering(
        const std::vector<std::vector<double>>& population,
        const std::vector<double>& fitness,
        double phi = 2.0)
    {
        int n = static_cast<int>(population.size());
        if (n == 0) return {};

        // Step 1: Build nearest-better graph
        std::vector<NBC_Edge> edges;
        for (int i = 0; i < n; i++) {
            double min_dist = std::numeric_limits<double>::max();
            int nb_idx = -1;

            for (int j = 0; j < n; j++) {
                if (fitness[j] > fitness[i] + 1e-9) { // strictly better (maximization)
                    double d = 0.0;
                    for (size_t k = 0; k < population[i].size(); k++) {
                        double diff = population[i][k] - population[j][k];
                        d += diff * diff;
                    }
                    d = std::sqrt(d);

                    if (d < min_dist) {
                        min_dist = d;
                        nb_idx = j;
                    }
                }
            }

            if (nb_idx != -1) {
                edges.push_back({ i, nb_idx, min_dist });
            }
        }

        // Step 2: Compute pruning threshold
        if (edges.empty()) {
            std::vector<std::vector<int>> clusters;
            for (int i = 0; i < n; i++) {
                clusters.push_back({ i });
            }
            return clusters;
        }

        double sum_weights = 0.0;
        for (const auto& e : edges) sum_weights += e.weight;
        double mean_weight = sum_weights / edges.size();
        double threshold = phi * mean_weight;

        // Step 3: Build adjacency list with pruned edges
        std::vector<std::vector<int>> adj(n);
        for (const auto& e : edges) {
            if (e.weight <= threshold) {
                adj[e.from].push_back(e.to);
                adj[e.to].push_back(e.from); // undirected for clustering
            }
        }

        // Step 4: Find connected components (clusters)
        std::vector<std::vector<int>> clusters;
        std::vector<bool> visited(n, false);

        for (int i = 0; i < n; i++) {
            if (!visited[i]) {
                std::vector<int> cluster;
                std::queue<int> q;
                q.push(i);
                visited[i] = true;

                while (!q.empty()) {
                    int u = q.front(); q.pop();
                    cluster.push_back(u);

                    for (int v : adj[u]) {
                        if (!visited[v]) {
                            visited[v] = true;
                            q.push(v);
                        }
                    }
                }

                if (!cluster.empty()) {
                    clusters.push_back(cluster);
                }
            }
        }

        return clusters;
    }

    // MPMMO BENCHMARK (FreePeaks wrapper with MPMMO semantics)
    class MPMMO_Benchmark {
    private:
        std::shared_ptr<FreePeaksWrapper> m_problem;
        int m_num_dms;
        std::vector<std::vector<int>> m_party_indices;

    public:
        MPMMO_Benchmark(std::shared_ptr<FreePeaksWrapper> prob, int num_dms)
            : m_problem(prob), m_num_dms(num_dms) {
            decomposeVariables();
        }

        void decomposeVariables() {
            int total_vars = m_problem->numberVariables();
            m_party_indices.resize(m_num_dms);

            int base_count = total_vars / m_num_dms;
            int remainder = total_vars % m_num_dms;

            int current_var = 0;
            for (int p = 0; p < m_num_dms; p++) {
                int count = base_count + (p < remainder ? 1 : 0);
                for (int k = 0; k < count; k++) {
                    m_party_indices[p].push_back(current_var++);
                }
            }
        }

        // Evaluate solution for a specific DM using context vectors
        // Using public evaluateSolution() instead of protected evaluate()
        double evaluateDM(int dm_id, const std::vector<double>& my_vars,
            const std::vector<std::vector<double>>& context)
        {
            // Build full solution from context vectors
            std::vector<double> full_sol(m_problem->numberVariables(), 0.5);

            // Insert active DM's variables
            const auto& my_indices = m_party_indices[dm_id];
            for (size_t i = 0; i < my_indices.size(); i++) {
                full_sol[my_indices[i]] = my_vars[i];
            }

            // Insert context from other DMs
            for (int other = 0; other < m_num_dms; other++) {
                if (other == dm_id) continue;
                const auto& other_indices = m_party_indices[other];
                for (size_t i = 0; i < other_indices.size(); i++) {
                    full_sol[other_indices[i]] = context[other][i];
                }
            }

            // Evaluate using FreePeaksWrapper's PUBLIC interface
            ofec::Solution<> sol(m_problem->numberObjectives(), m_problem->numberConstraints());
            sol.variable().resize(full_sol.size());
            for (size_t i = 0; i < full_sol.size(); i++) {
                sol.variable()[i] = static_cast<ofec::Real>(full_sol[i]);
            }
            m_problem->evaluateSolution(sol);  // PUBLIC method to avoid protected access error
            return static_cast<double>(sol.objective(0));
        }

        int getNumDMs() const { return m_num_dms; }
        int getNumVariables() const { return m_problem->numberVariables(); }
        const std::vector<int>& getPartyIndices(int p) const { return m_party_indices[p]; }
    };

    // MPM-CoEA SOLVER (Two-tiered co-evolutionary framework)
    class MPM_CoEA {
    private:
        std::shared_ptr<MPMMO_Benchmark> m_bench;
        int m_pop_size;
        int m_med_pop_size;
        int m_K;

        // Competing populations: [DM][subpopulation_id][individual][variable]
        std::vector<std::vector<std::vector<std::vector<double>>>> m_competing_pops;
        std::vector<std::vector<std::vector<double>>> m_competing_fitness;

        // Mediating population: [individual][variable]
        std::vector<std::vector<double>> m_mediating_pop;
        std::vector<double> m_mediating_fitness;

        // Context vectors: [DM][variable]
        std::vector<std::vector<double>> m_context_vectors;

    public:
        MPM_CoEA(std::shared_ptr<MPMMO_Benchmark> bench, int pop_size,
            int med_pop_size, int K)
            : m_bench(bench), m_pop_size(pop_size),
            m_med_pop_size(med_pop_size), m_K(K) {
        }

        void initialize() {
            int num_dms = m_bench->getNumDMs();
            int dim = m_bench->getNumVariables();

            std::cout << " [MPM-CoEA] Initializing with " << num_dms
                << " decision-makers, " << dim << " dimensions" << std::endl;

            // Step 1: Create large initial population for NBC
            const int NBC_POP_SIZE = 500;
            std::vector<std::vector<double>> global_pop(NBC_POP_SIZE, std::vector<double>(dim));
            std::vector<std::vector<double>> global_fitness(num_dms,
                std::vector<double>(NBC_POP_SIZE));

            // Random initialization in [0,1]^D
            for (int i = 0; i < NBC_POP_SIZE; i++) {
                for (int d = 0; d < dim; d++) {
                    global_pop[i][d] = (double)rand() / RAND_MAX;
                }
                // Evaluate for all DMs using neutral context (0.5 for all variables)
                std::vector<std::vector<double>> neutral_context(num_dms, std::vector<double>(1, 0.5));
                for (int dm = 0; dm < num_dms; dm++) {
                    global_fitness[dm][i] = m_bench->evaluateDM(dm, global_pop[i], neutral_context);
                }
            }

            // Step 2: Apply NBC per DM to create subpopulations
            m_competing_pops.resize(num_dms);
            m_competing_fitness.resize(num_dms);

            for (int dm = 0; dm < num_dms; dm++) {
                auto clusters = nearestBetterClustering(global_pop, global_fitness[dm], 2.0);

                m_competing_pops[dm].resize(clusters.size());
                m_competing_fitness[dm].resize(clusters.size());

                for (size_t c = 0; c < clusters.size(); c++) {
                    int take = std::min((int)clusters[c].size(), m_pop_size);
                    m_competing_pops[dm][c].resize(take, std::vector<double>(dim));
                    m_competing_fitness[dm][c].resize(take);

                    for (int i = 0; i < take; i++) {
                        m_competing_pops[dm][c][i] = global_pop[clusters[c][i]];
                        m_competing_fitness[dm][c][i] = global_fitness[dm][clusters[c][i]];
                    }
                }
                std::cout << "  DM " << dm << ": " << clusters.size() << " niches detected" << std::endl;
            }

            // Step 3: Initialize mediating population
            m_mediating_pop.resize(m_med_pop_size, std::vector<double>(dim));
            m_mediating_fitness.resize(m_med_pop_size);

            // Collect top candidates
            std::vector<std::pair<double, std::vector<double>>> candidates;
            for (int dm = 0; dm < num_dms; dm++) {
                for (size_t c = 0; c < m_competing_pops[dm].size(); c++) {
                    for (size_t i = 0; i < m_competing_pops[dm][c].size(); i++) {
                        candidates.emplace_back(m_competing_fitness[dm][c][i],
                            m_competing_pops[dm][c][i]);
                    }
                }
            }

            // Sort by fitness (descending for maximization)
            std::sort(candidates.begin(), candidates.end(),
                [](const auto& a, const auto& b) { return a.first > b.first; });

            // Select diverse top solutions
            int selected = 0;
            for (size_t i = 0; i < candidates.size() && selected < m_med_pop_size; i++) {
                // Simple diversity check
                bool is_duplicate = false;
                for (int j = 0; j < selected; j++) {
                    double dist = 0.0;
                    for (int d = 0; d < dim; d++) {
                        double diff = candidates[i].second[d] - m_mediating_pop[j][d];
                        dist += diff * diff;
                    }
                    if (std::sqrt(dist) < 0.01) {
                        is_duplicate = true;
                        break;
                    }
                }
                if (!is_duplicate) {
                    m_mediating_pop[selected] = candidates[i].second;
                    m_mediating_fitness[selected] = candidates[i].first;
                    selected++;
                }
            }

            // Fill remaining slots if needed
            while (selected < m_med_pop_size) {
                for (int d = 0; d < dim; d++) {
                    m_mediating_pop[selected][d] = (double)rand() / RAND_MAX;
                }
                std::vector<std::vector<double>> neutral_context(num_dms, std::vector<double>(1, 0.5));
                m_mediating_fitness[selected] = m_bench->evaluateDM(0, m_mediating_pop[selected], neutral_context);
                selected++;
            }

            // Step 4: Initialize context vectors
            m_context_vectors.resize(num_dms);
            for (int dm = 0; dm < num_dms; dm++) {
                int best_subpop = 0, best_indiv = 0;
                double best_fit = -1e9;

                for (size_t c = 0; c < m_competing_pops[dm].size(); c++) {
                    for (size_t i = 0; i < m_competing_pops[dm][c].size(); i++) {
                        if (m_competing_fitness[dm][c][i] > best_fit) {
                            best_fit = m_competing_fitness[dm][c][i];
                            best_subpop = static_cast<int>(c);
                            best_indiv = static_cast<int>(i);
                        }
                    }
                }

                const auto& indices = m_bench->getPartyIndices(dm);
                m_context_vectors[dm].resize(indices.size());
                for (size_t i = 0; i < indices.size(); i++) {
                    m_context_vectors[dm][i] = m_competing_pops[dm][best_subpop][best_indiv][indices[i]];
                }
            }

            std::cout << " [MPM-CoEA] Initialization complete" << std::endl;
        }

        void run_one_generation(int gen) {
            int num_dms = m_bench->getNumDMs();
            int dim = m_bench->getNumVariables();

            // Evolve competing populations (simple mutation)
            for (int dm = 0; dm < num_dms; dm++) {
                for (size_t subpop = 0; subpop < m_competing_pops[dm].size(); subpop++) {
                    for (size_t i = 0; i < m_competing_pops[dm][subpop].size(); i++) {
                        if (i == 0) continue; // Keep best individual

                        for (int d = 0; d < dim; d++) {
                            m_competing_pops[dm][subpop][i][d] +=
                                ((double)rand() / RAND_MAX - 0.5) * 0.1;
                            // Boundary handling [0,1]
                            if (m_competing_pops[dm][subpop][i][d] < 0) m_competing_pops[dm][subpop][i][d] = 0;
                            if (m_competing_pops[dm][subpop][i][d] > 1) m_competing_pops[dm][subpop][i][d] = 1;
                        }

                        // Re-evaluate with current context vectors
                        m_competing_fitness[dm][subpop][i] =
                            m_bench->evaluateDM(dm, m_competing_pops[dm][subpop][i], m_context_vectors);
                    }

                    // Sort subpopulation by fitness
                    std::vector<size_t> indices(m_competing_pops[dm][subpop].size());
                    std::iota(indices.begin(), indices.end(), 0);
                    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
                        return m_competing_fitness[dm][subpop][a] > m_competing_fitness[dm][subpop][b];
                        });

                    // Reorder population
                    std::vector<std::vector<double>> new_pop;
                    std::vector<double> new_fit;
                    for (size_t idx : indices) {
                        new_pop.push_back(m_competing_pops[dm][subpop][idx]);
                        new_fit.push_back(m_competing_fitness[dm][subpop][idx]);
                    }
                    m_competing_pops[dm][subpop] = new_pop;
                    m_competing_fitness[dm][subpop] = new_fit;
                }
            }

            // Update context vectors
            for (int dm = 0; dm < num_dms; dm++) {
                double best_fit = -1e9;
                std::vector<double> best_vars;

                for (size_t subpop = 0; subpop < m_competing_pops[dm].size(); subpop++) {
                    if (m_competing_fitness[dm][subpop][0] > best_fit) {
                        best_fit = m_competing_fitness[dm][subpop][0];
                        const auto& indices = m_bench->getPartyIndices(dm);
                        best_vars.resize(indices.size());
                        for (size_t i = 0; i < indices.size(); i++) {
                            best_vars[i] = m_competing_pops[dm][subpop][0][indices[i]];
                        }
                    }
                }
                m_context_vectors[dm] = best_vars;
            }

            // Evolve mediating population with DE
            std::vector<std::vector<double>> new_med_pop = m_mediating_pop;
            std::vector<double> new_med_fit(m_med_pop_size);

            for (int i = 0; i < m_med_pop_size; i++) {
                // DE/rand/1/bin mutation
                int r1 = rand() % m_med_pop_size;
                int r2 = rand() % m_med_pop_size;
                int r3 = rand() % m_med_pop_size;
                while (r2 == r1) r2 = rand() % m_med_pop_size;
                while (r3 == r1 || r3 == r2) r3 = rand() % m_med_pop_size;

                double F = 0.5, CR = 0.9;
                for (int d = 0; d < dim; d++) {
                    if ((rand() / (double)RAND_MAX) < CR || d == dim - 1) {
                        new_med_pop[i][d] = m_mediating_pop[r1][d] +
                            F * (m_mediating_pop[r2][d] - m_mediating_pop[r3][d]);
                        // Boundary handling
                        if (new_med_pop[i][d] < 0) new_med_pop[i][d] = 0;
                        if (new_med_pop[i][d] > 1) new_med_pop[i][d] = 1;
                    }
                }

                // Evaluate consensus fitness
                std::vector<std::vector<double>> neutral_context(num_dms, std::vector<double>(1, 0.5));
                new_med_fit[i] = m_bench->evaluateDM(0, new_med_pop[i], neutral_context);
            }

            // Selection
            for (int i = 0; i < m_med_pop_size; i++) {
                if (new_med_fit[i] > m_mediating_fitness[i]) {
                    m_mediating_pop[i] = new_med_pop[i];
                    m_mediating_fitness[i] = new_med_fit[i];
                }
            }
        }

        double getBestConsensusFitness() const {
            if (m_mediating_fitness.empty()) return -1e9;
            return *std::max_element(m_mediating_fitness.begin(), m_mediating_fitness.end());
        }
    };

    // MAIN FUNCTION
    void run_mpmcoea() {
        srand(static_cast<unsigned int>(time(nullptr)));

        // Working directory
        ofec::g_working_directory = R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Data\ofec_data_new)";

        // Create FreePeaks instance via wrapper
        auto free_peaks = std::make_shared<FreePeaksWrapper>();
        free_peaks->configure("sop/2d_2_s1_1_1.txt");

        // Initialize with dummy environment
        class DummyEnv : public ofec::Environment {
        public:
            DummyEnv() = default;
        };
        DummyEnv env;
        free_peaks->initialize(&env);

        std::cout << " [MPM-CoEA] Loaded FreePeaks instance" << std::endl;
        std::cout << " Dimensions: " << free_peaks->numberVariables() << std::endl;

        // Create benchmark and solver
        auto benchmark = std::make_shared<MPMMO_Benchmark>(free_peaks, 2);
        MPM_CoEA solver(benchmark, 50, 100, 3);  // pop_size=50, med_pop=100, K=3
        solver.initialize();

        std::cout << "\ [MPM-CoEA] Starting optimization (20 generations)" << std::endl;
        for (int gen = 0; gen < 20; gen++) {
            solver.run_one_generation(gen);
            std::cout << "Gen " << gen + 1 << " | Best Consensus Fitness: "
                << std::fixed << std::setprecision(4) << solver.getBestConsensusFitness() << std::endl;
        }
        std::cout << "\n [MPM-CoEA] Optimization complete" << std::endl;
    }
}

#endif // !OFEC_CUSTOM_METHOD_HPP */

#ifndef OFEC_CUSTOM_METHOD_HPP
#define OFEC_CUSTOM_METHOD_HPP

#include "../core/global.h"
#include "interface.h"




#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../instance/problem/continuous/free_peaks/subproblem/subproblem.h"
#include "../instance/problem/continuous/free_peaks/subproblem/function/one_peak_function.h"
#include "../instance/problem/continuous/free_peaks/factory.h"


namespace fs = std::filesystem;


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
        std::cout << parent << ":\n";  // 打印父节点
        for (const auto& [child, value] : children) {
            std::cout << "  " << child << " -> " << value << "\n";  // 打印子节点及值
        }
        std::cout << "\n";  // 每个父节点空一行，便于阅读
    }
}




// 网格采样函数
void gridSampling(ofec::Environment* env, const std::string& filename,
    int grid_size = 1e2) {

    using namespace ofec;


    // 获取问题域
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

    // 打开输出文件
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // 设置输出精度
    outfile << std::fixed << std::setprecision(10);

    // 写入文件头
    outfile << "x1,x2,f(x1,x2)" << std::endl;

    // 获取变量范围
    double x1_min = domain[0].limit.first;
    double x1_max = domain[0].limit.second;
    double x2_min = domain[1].limit.first;
    double x2_max = domain[1].limit.second;

    // 计算步长
    double dx1 = (x1_max - x1_min) / (grid_size - 1);
    double dx2 = (x2_max - x2_min) / (grid_size - 1);

    std::cout << "\n=== Grid Sampling Parameters ===" << std::endl;
    std::cout << "Grid size: " << grid_size << "x" << grid_size
        << " = " << grid_size * grid_size << " points" << std::endl;
    std::cout << "x1 step: " << dx1 << std::endl;
    std::cout << "x2 step: " << dx2 << std::endl;

    // 创建解决方案对象
    std::shared_ptr<ofec::SolutionBase> sol(env->problem()->createSolution());
    auto& solvec = dynamic_cast<ofec::Solution<>&>(*sol).variable().vector();
    int total_points = grid_size * grid_size;
    int current_point = 0;

    // 网格采样
    std::cout << "\n=== Start Grid Sampling ===" << std::endl;

    for (int i = 0; i < grid_size; ++i) {
        double x1 = x1_min + i * dx1;

        for (int j = 0; j < grid_size; ++j) {
            double x2 = x2_min + j * dx2;

            // 设置变量值
            //std::vector<double> vars = { x1, x2 };

            solvec = { x1,x2 };
            // 评估函数
            sol->evaluate(env, false);

            // 写入文件
            outfile << x1 << "," << x2 << "," << sol->objective().front() << std::endl;

            // 进度显示
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

                // there are no parameters to set for Euclidean distance, but the following step can not be ingored , because the register name will be saved.
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
        /*ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
        ofec::g_working_directory = "//172.29.41.69/share/2018/diaoyiya/ofec-data";
        ofec::g_working_directory = "//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data";
        ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
        ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
        //ofec::g_working_directory = "//172.24.24.151/e/DiaoYiya/code/data/ofec-data";
        ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/";
        ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
        ofec::g_working_directory = "D:/code/Data/ofec_data/";*/



        // multiparty multiobjective
        // 
        // 
        ofec::g_working_directory = "E:/HITSZ/Research/Multimodal_Multiparty_Optimization/ThesisProject/Data/ofec_data_new";

        //ofec::g_working_directory = "/home/dyy/data/ofec_data/ofec_data/";
        //ofec::g_working_directory = "/data/Share/Student/2018/diaoyiya/data/ofec_data/";
      //  ofec::g_working_directory = "/mnt/dataShare/Student/2018/diaoyiya/data/ofec_data/";

        registerInstance();

        //generateMPMMOBenchmark("mpmmo_2p2g1l", 2, 2, 1, 2);
        generateHandMadeFreePeak("mpmmo_1");


        //std::string ofec_dir = OFEC_DIR;
        //std::string savedir = "D:/code/multi_party_multi_modal/freepeak_ofec";
        //std::string list_file = ofec_dir + std::string("/SrcFiles.cmake");
        //copyFile(ofec_dir, list_file,savedir);
    }
}

#endif // !OFEC_CUSTOM_METHOD_HPP

/*#ifndef OFEC_CUSTOM_METHOD_HPP
#define OFEC_CUSTOM_METHOD_HPP

#include "../core/global.h"
#include "interface.h"

// --- Includes for Supervisor's Map Generation ---
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../instance/problem/continuous/free_peaks/subproblem/subproblem.h"
#include "../instance/problem/continuous/free_peaks/subproblem/function/one_peak_function.h"
#include "../instance/problem/continuous/free_peaks/factory.h"

// --- Includes for Your Solver ---
#include "../core/problem/solution.h"
#include <vector>
#include <numeric>   
#include <algorithm> 
#include <random>

namespace fs = std::filesystem;

namespace ofec {

    // ==============================================================================
    // PHASE 1: SUPERVISOR'S GENERATION CODE
    // ==============================================================================
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

        ParameterMap freepeak_param;
        freepeak_param["generation_type"] = std::string("assigned");
        freepeak_param["dataFile1"] = dirname + "/" + funname + ".txt";
        freepeak->inputParameters().input(freepeak_param);
        freepeak->initialize(env.get());

        std::shared_ptr<ofec::Random> pro_rnd(new Random(0.5));
        freepeak->setRandom(pro_rnd);
        freepeak->setSizes(numDim, numObj, numCon);

        freepeak->setKDtree({
            {"root", { { "subtree1", 0.5}, {"subtree2", 0.5} } },
            {"subtree1", { {funname + "_s1_1", 0.1}, {funname + "_s1_2", 0.4} } },
            {"subtree2", { {funname + "_s2_1", 0.3}, {funname + "_s2_2", 0.2} } } });

        for (int idsubtree(1); idsubtree <= 2; ++idsubtree) {
            for (int idpeak(1); idpeak <= 2; ++idpeak) {
                std::string subspace_name = funname + "_s" + std::to_string(idsubtree) + "_" + std::to_string(idpeak);
                ParameterMap subpro_param;
                subpro_param["subspace"] = subspace_name;
                subpro_param["generation_type"] = std::string("assigned");
                subpro_param["dataFile1"] = dirname + "/" + subspace_name + ".txt";

                auto subpro(Subproblem::create());
                subpro->initialize(subpro_param, freepeak);

                auto dis(FactoryFP<DistanceBase>::produce("Euclidean"));
                ParameterMap dis_param;
                dis->initialize(freepeak, subspace_name, dis_param);
                subpro->setDistance(dis);

                ParameterMap fun_param;
                fun_param["generation_type"] = std::string("assigned");
                fun_param["dataFile1"] = dirname + "/" + subspace_name + "_onepeak" + ".txt";
                auto func(FactoryFP<FunctionBase>::produce("one_peak"));
                func->initialize(freepeak, subspace_name, fun_param);
                auto onepeak_func = dynamic_cast<ofec::free_peaks::OnePeakFunction*>(func);

                auto onepeak(FactoryFP<OnePeakBase>::produceRandom(rnd->uniform.next()));
                ParameterMap onepeak_param;
                onepeak_param["center_type"] = std::string("assigned");
                onepeak_param["height"] = rnd->uniform.nextNonStd<Real>(0, 100);
                std::vector<Real> peak_position(numDim);
                for (int idx(0); idx < numDim; ++idx) {
                    // Using explicit valid bounds instead of uninitialized static variables
                    //peak_position[idx] = rnd->uniform.nextNonStd<Real>(-100.0, 100.0);
                    peak_position[idx] = rnd->uniform.nextNonStd<Real>(OnePeakBase::ms_x_range.first, OnePeakBase::ms_x_range.second);
                }
                onepeak_param["center_postion"] = peak_position;
                onepeak->initialize(freepeak, subspace_name, onepeak_param);
                onepeak_func->addOnePeaks(onepeak);
                subpro->setFunction(func);

                ofec::ParameterMap default_param;
                auto transformX(FactoryFP<X_TransformBase>::produceRandom(rnd->uniform.next()));
                transformX->initialize(freepeak, subspace_name, default_param);
                subpro->addVariableTransform(transformX);

                auto transformY(FactoryFP<Y_TransformBase>::produce("map_objective"));
                transformY->initialize(freepeak, subspace_name, default_param);
                subpro->addObjectiveTransform(transformY);

                freepeak->setSubproblem(subspace_name, subpro);
            }
        }

        freepeak->bindData();

        std::filesystem::create_directories(FreePeaks::directory() + dirname);
        std::filesystem::create_directories(Subproblem::directory() + dirname);
        std::filesystem::create_directories(OnePeakFunction::directory() + dirname);

        freepeak->inputParameters().at("generation_type")->setValue("read_file");
        for (auto& it : freepeak->subspaceTree().name_box_subproblem) {
            if (it.second.second != nullptr) {
                it.second.second->inputParameters().at("generation_type")->setValue("read_file");
                it.second.second->function()->inputParameters().at("generation_type")->setValue("read_file");
            }
        }

        freepeak->recordInputParameters();
        freepeak->outputTotalFile();
        std::cout << ">>> Landscape [" << funname << "] successfully generated and saved to disk." << std::endl;
    }


    // ==============================================================================
    // PHASE 2: YOUR MPMMO ALGORITHM CODE
    // ==============================================================================
    class FreePeaksManual : public FreePeaks {
    public:
        void configure(const std::string& filename) {
            this->m_file_name = filename;
            this->m_generation_type = "read_file";
            this->m_number_objectives = 1;
        }
        void initialize_(Environment* env) override {
            FreePeaks::initialize_(env);
        }
        void evaluateSolution(ofec::Solution<>& sol) {
            std::vector<ofec::Real> objs(this->numberObjectives());
            std::vector<ofec::Real> cons(this->numberConstraints());
            this->evaluate(sol.variable(), objs, cons);
            for (size_t i = 0; i < objs.size(); i++) sol.objective(i) = objs[i];
        }
    };

    class MPMMO_Benchmark {
    private:
        std::shared_ptr<FreePeaksManual> m_problem;
        int m_num_parties;
        std::vector<std::vector<int>> m_party_indices;

    public:
        MPMMO_Benchmark(std::shared_ptr<FreePeaksManual> prob, int num_parties)
            : m_problem(prob), m_num_parties(num_parties) {
            decomposeVariables();
        }

        void decomposeVariables() {
            int total_vars = static_cast<int>(m_problem->numberVariables()); // Fixed warning
            m_party_indices.resize(m_num_parties);
            int current_var = 0;
            for (int p = 0; p < m_num_parties; p++) {
                int count = total_vars / m_num_parties;
                if (p < total_vars % m_num_parties) count++;
                for (int k = 0; k < count; k++) {
                    m_party_indices[p].push_back(current_var++);
                }
            }
        }

        double evaluateParty(int party_id, const std::vector<double>& my_vars, const std::vector<std::vector<double>>& context) {
            ofec::Solution<> full_sol(m_problem->numberObjectives(), m_problem->numberConstraints());
            full_sol.variable().resize(m_problem->numberVariables());

            const auto& my_indices = m_party_indices[party_id];
            for (size_t i = 0; i < my_indices.size(); i++) full_sol.variable()[my_indices[i]] = my_vars[i];

            for (int other = 0; other < m_num_parties; other++) {
                if (other == party_id) continue;
                const auto& other_indices = m_party_indices[other];
                for (size_t i = 0; i < other_indices.size(); i++) full_sol.variable()[other_indices[i]] = context[other][i];
            }
            m_problem->evaluateSolution(full_sol);
            return full_sol.objective(0);
        }

        const std::vector<int>& getPartyIndices(int p) { return m_party_indices[p]; }
        int getNumParties() { return m_num_parties; }
        int getTotalVariables() const { return static_cast<int>(m_problem->numberVariables()); } // Fixed warning
    };

    class MPMMO_Solver {
        std::shared_ptr<MPMMO_Benchmark> m_bench;
        int m_pop_size;
        std::vector<std::vector<std::vector<double>>> m_populations;
        std::vector<std::vector<double>> m_fitness_vals;
        std::vector<std::vector<double>> m_representatives;

    public:
        MPMMO_Solver(std::shared_ptr<MPMMO_Benchmark> b, int pop) : m_bench(b), m_pop_size(pop) {}

        void initializeWithNBC() {
            int n_parties = m_bench->getNumParties();
            m_populations.resize(n_parties);
            m_fitness_vals.resize(n_parties);
            m_representatives.resize(n_parties);

            const int NBC_POP_SIZE = 500;
            std::vector<std::vector<double>> global_pop(NBC_POP_SIZE);
            std::vector<double> global_fitness(NBC_POP_SIZE);

            for (int i = 0; i < NBC_POP_SIZE; i++) {
                global_pop[i].resize(m_bench->getTotalVariables());
                for (size_t v = 0; v < global_pop[i].size(); v++) global_pop[i][v] = (double)rand() / RAND_MAX;

                std::vector<std::vector<double>> neutral_context(n_parties);
                for (int p = 0; p < n_parties; p++) neutral_context[p] = std::vector<double>(m_bench->getPartyIndices(p).size(), 0.5);
                global_fitness[i] = m_bench->evaluateParty(0, global_pop[i], neutral_context);
            }

            std::vector<int> indices(NBC_POP_SIZE);
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](int a, int b) { return global_fitness[a] > global_fitness[b]; });

            for (int p = 0; p < n_parties; p++) {
                int n_vars = static_cast<int>(m_bench->getPartyIndices(p).size()); // Fixed warning
                m_populations[p].resize(m_pop_size, std::vector<double>(n_vars));
                m_fitness_vals[p].resize(m_pop_size);

                for (int i = 0; i < m_pop_size; i++) {
                    int src_idx = indices[(p * m_pop_size + i) % NBC_POP_SIZE];
                    const auto& my_indices = m_bench->getPartyIndices(p);
                    for (size_t v = 0; v < my_indices.size(); v++) m_populations[p][i][v] = global_pop[src_idx][my_indices[v]];
                }
                m_representatives[p] = m_populations[p][0];
            }
        }

        void run_one_generation(int gen) {
            int n_parties = m_bench->getNumParties();
            for (int active_p = 0; active_p < n_parties; active_p++) {
                double best_val = -1e9;
                int best_idx = 0;

                for (int i = 0; i < m_pop_size; i++) {
                    double fit = m_bench->evaluateParty(active_p, m_populations[active_p][i], m_representatives);
                    m_fitness_vals[active_p][i] = fit;
                    if (fit > best_val) { best_val = fit; best_idx = i; }
                }

                m_representatives[active_p] = m_populations[active_p][best_idx];

                for (int i = 0; i < m_pop_size; i++) {
                    if (i == best_idx) continue;
                    for (size_t v = 0; v < m_populations[active_p][i].size(); v++) {
                        m_populations[active_p][i][v] += ((double)rand() / RAND_MAX - 0.5) * 0.1;
                        if (m_populations[active_p][i][v] < 0) m_populations[active_p][i][v] = 0;
                        if (m_populations[active_p][i][v] > 1) m_populations[active_p][i][v] = 1;
                    }
                }
            }
        }

        double getGlobalBest() { return m_bench->evaluateParty(0, m_representatives[0], m_representatives); }

        std::vector<double> getGlobalBestVariables() {
            std::vector<double> full_sol(m_bench->getTotalVariables());
            for (int p = 0; p < m_bench->getNumParties(); p++) {
                const auto& my_indices = m_bench->getPartyIndices(p);
                for (size_t i = 0; i < my_indices.size(); i++) full_sol[my_indices[i]] = m_representatives[p][i];
            }
            return full_sol;
        }
    };


    // ==============================================================================
    // MAIN RUN FUNCTION
    // ==============================================================================
    void run(int argc, char* argv[]) {
        using namespace ofec;
        using namespace std;

        ofec::g_working_directory = "E:/HITSZ/Research/Multimodal_Multiparty_Optimization/ThesisProject/Data/ofec_data_new/";

        registerInstance();

        cout << "=========================================================" << endl;
        cout << " PHASE 1: GENERATING THE BENCHMARK (Dynamic Map Builder) " << endl;
        cout << "=========================================================" << endl;

        std::string benchmark_name = "mpmmo_test_run";
        generateHandMadeFreePeak(benchmark_name);

        cout << "\n=========================================================" << endl;
        cout << " PHASE 2: RUNNING MPCoEA ALGORITHM (Multimodal Solver)   " << endl;
        cout << "=========================================================" << endl;

        std::string target_file = "multiparty_multimodal/" + benchmark_name + ".txt";

        auto core_problem = std::make_shared<FreePeaksManual>();
        core_problem->configure(target_file);

        try {
            core_problem->initialize(nullptr);

            MPMMO_Benchmark benchmark(core_problem, 2);
            MPMMO_Solver solver(std::make_shared<MPMMO_Benchmark>(benchmark), 50);

            solver.initializeWithNBC();

            for (int g = 0; g < 20; g++) {
                solver.run_one_generation(g);
                std::vector<double> best_vars = solver.getGlobalBestVariables();

                cout << "Gen " << g + 1 << " | Best Fitness: " << solver.getGlobalBest()
                    << " | Coordinates: [ ";
                for (double v : best_vars) cout << v << " ";
                cout << "]" << endl;
            }
            cout << ">>> Full Pipeline Execution Complete." << endl;

        }
        catch (const std::exception& e) {
            cerr << "Error: " << e.what() << endl;
        }
        system("pause");
    }
}

#endif // !OFEC_CUSTOM_METHOD_HPP*/
