#define OFEC_MPMMO_EXPORTS
#include "ofec_mpmmo_api.h"

#include "../run/interface.h"
#include "../core/problem/solution.h"
#include "../instance/problem/continuous/free_peaks/free_peaks_multiparty.h"

#include <algorithm>
#include <exception>
#include <memory>
#include <string>
#include <vector>

namespace {
    thread_local std::string g_last_error;

    struct MPMMOHandle {
        std::shared_ptr<ofec::Environment> env;
        ofec::FreePeaksMultiParty* problem = nullptr;
        int dimension = 2;
    };

    void set_error(const std::string& msg) {
        g_last_error = msg;
    }
}

extern "C" {

OFEC_MPMMO_API void* ofec_mpmmo_create(int suite_id, int dimension) {
    try {
        g_last_error.clear();
        ofec::registerInstance();

        auto* raw_env = ofec::generateEnvironmentByFactory("free_peaks_multiparty");
        if (!raw_env) throw std::runtime_error("generateEnvironmentByFactory failed");

        auto handle = std::make_unique<MPMMOHandle>();
        handle->env.reset(raw_env);
        handle->env->setProblem(ofec::generateProblemByFactory("free_peaks_multiparty"));

        ofec::ParameterMap pm;
        pm["suite_id"] = suite_id;
        pm["problem_dimension"] = static_cast<size_t>(std::max(2, dimension));
        handle->env->problem()->inputParameters().input(pm);
        handle->env->initializeProblem(0.5);

        handle->problem = dynamic_cast<ofec::FreePeaksMultiParty*>(handle->env->problem());
        if (!handle->problem) throw std::runtime_error("free_peaks_multiparty initialization failed");
        handle->dimension = static_cast<int>(handle->problem->numberVariables());
        return handle.release();
    }
    catch (const std::exception& e) {
        set_error(e.what());
        return nullptr;
    }
    catch (...) {
        set_error("unknown exception in ofec_mpmmo_create");
        return nullptr;
    }
}

OFEC_MPMMO_API void ofec_mpmmo_destroy(void* handle) {
    delete static_cast<MPMMOHandle*>(handle);
}

OFEC_MPMMO_API int ofec_mpmmo_evaluate(void* handle, const double* decisions, int n_points, double* objectives) {
    try {
        g_last_error.clear();
        auto* h = static_cast<MPMMOHandle*>(handle);
        if (!h || !h->env || !h->problem) throw std::runtime_error("invalid OFEC MPMMO handle");
        if (!decisions || !objectives || n_points < 0) throw std::runtime_error("invalid evaluate buffers");

        auto sol = std::unique_ptr<ofec::SolutionBase>(h->problem->createSolution());
        auto& x = dynamic_cast<ofec::Solution<ofec::VariableVector<ofec::Real>>&>(*sol).variable().vector();
        const int dim = h->dimension;
        for (int i = 0; i < n_points; ++i) {
            for (int d = 0; d < dim; ++d)
                x[d] = decisions[i * dim + d];
            sol->evaluate(h->env.get(), false);
            objectives[i * 2 + 0] = sol->objective(0);
            objectives[i * 2 + 1] = sol->objective(1);
        }
        return 1;
    }
    catch (const std::exception& e) {
        set_error(e.what());
        return 0;
    }
    catch (...) {
        set_error("unknown exception in ofec_mpmmo_evaluate");
        return 0;
    }
}

OFEC_MPMMO_API int ofec_mpmmo_shared_count(void* handle) {
    auto* h = static_cast<MPMMOHandle*>(handle);
    if (!h || !h->problem) return -1;
    return static_cast<int>(h->problem->currentSpec().shared_optima.size());
}

OFEC_MPMMO_API int ofec_mpmmo_shared_optima(void* handle, double* out_points) {
    try {
        auto* h = static_cast<MPMMOHandle*>(handle);
        if (!h || !h->problem || !out_points) throw std::runtime_error("invalid shared optima buffers");
        const auto& optima = h->problem->currentSpec().shared_optima;
        const int dim = h->dimension;
        for (size_t i = 0; i < optima.size(); ++i) {
            for (int d = 0; d < dim; ++d)
                out_points[i * dim + d] = d < static_cast<int>(optima[i].size()) ? optima[i][d] : 0.5;
        }
        return static_cast<int>(optima.size());
    }
    catch (const std::exception& e) {
        set_error(e.what());
        return -1;
    }
    catch (...) {
        set_error("unknown exception in ofec_mpmmo_shared_optima");
        return -1;
    }
}

OFEC_MPMMO_API const char* ofec_mpmmo_last_error() {
    return g_last_error.c_str();
}

}
