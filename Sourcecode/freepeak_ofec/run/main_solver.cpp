#include "mpmcoea_solver.hpp"
#include "../instance/problem/continuous/free_peaks/factory.h"
#include "free_peaks_multiparty_sampler.hpp"
#include "multiparty_nbn_cmaes_solver.hpp"

#include <algorithm>
#include <cstdlib>
#include <string>

namespace {

ofec::multiparty_nbn_cmaes::Config parseNbnConfig(int argc, char* argv[], int first) {
    ofec::multiparty_nbn_cmaes::Config config;
    if (argc > first) config.iterations = std::stoi(argv[first]);
    if (argc > first + 1) config.party_pop_size = std::stoi(argv[first + 1]);
    if (argc > first + 2) config.coord_pop_size = std::stoi(argv[first + 2]);
    if (argc > first + 3) config.nbn_source_cap = std::stoi(argv[first + 3]);
    if (argc > first + 4) config.max_evaluations = std::stoi(argv[first + 4]);
    config.coord_children = std::max(20, static_cast<int>(0.55 * config.coord_pop_size));
    config.niche_lambda = std::max(6, static_cast<int>(0.12 * config.party_pop_size));
    return config;
}

}

int main(int argc, char* argv[]) {
    try {
        if (argc > 1 && std::string(argv[1]) == "--sample-multiparty") {
            std::cout << ">>> [FreePeaksMultiParty] Sampling 12 suites\n" << std::endl;
            ofec::sampleFreePeaksMultiPartySuites();
            return 0;
        }
        if (argc > 1 && std::string(argv[1]) == "--run-multiparty-nbn-cmaes") {
            int suite_id = 0;
            size_t dimension = 2;
            std::string output_root = "Visualization/multiparty_nbn_cmaes";
            if (argc > 2) suite_id = std::stoi(argv[2]);
            if (argc > 3) dimension = static_cast<size_t>(std::stoul(argv[3]));
            if (argc > 4) output_root = argv[4];
            auto config = parseNbnConfig(argc, argv, 5);
            std::cout << ">>> [MultiParty-NBN-CMAES] Running suite "
                      << (suite_id == 0 ? std::string("1..12") : std::to_string(suite_id))
                      << ", D=" << dimension
                      << ", iterations=" << config.iterations
                      << ", MaxFE=" << config.max_evaluations
                      << ", partyPop=" << config.party_pop_size
                      << ", coordPop=" << config.coord_pop_size
                      << "\n" << std::endl;
            ofec::multiparty_nbn_cmaes::runMultiPartyNbnCmaes(suite_id, dimension, output_root, config);
            return 0;
        }
        if (argc > 1 && std::string(argv[1]) == "--run-multiparty-nbn-cmaes-all-dims") {
            std::string output_root = "Visualization/multiparty_nbn_cmaes";
            if (argc > 2) output_root = argv[2];
            auto config = parseNbnConfig(argc, argv, 3);
            std::cout << ">>> [MultiParty-NBN-CMAES] Running suites 1..12, D in {2,3,5,10}"
                      << ", iterations=" << config.iterations
                      << ", MaxFE=" << config.max_evaluations
                      << ", partyPop=" << config.party_pop_size
                      << ", coordPop=" << config.coord_pop_size
                      << "\n" << std::endl;
            ofec::multiparty_nbn_cmaes::runMultiPartyNbnCmaesAllDimensions(output_root, config);
            return 0;
        }

        std::cout << ">>> [MPM-CoEA] Starting Experiment with Visualization\n" << std::endl;
        ofec::runMPMCoEAExperiment();
    }
    catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
