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


#endif // !OFEC_CUSTOM_METHOD_HPP
