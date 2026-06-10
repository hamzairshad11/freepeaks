#define NOMINMAX
#include "mex.h"
#include <windows.h>
#include <algorithm>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using CreateFn = void* (*)(int, int);
using DestroyFn = void (*)(void*);
using EvaluateFn = int (*)(void*, const double*, int, double*);
using CountFn = int (*)(void*);
using SharedFn = int (*)(void*, double*);
using ErrorFn = const char* (*)();

struct Api {
    HMODULE dll = nullptr;
    CreateFn create = nullptr;
    DestroyFn destroy = nullptr;
    EvaluateFn evaluate = nullptr;
    CountFn count = nullptr;
    SharedFn shared = nullptr;
    ErrorFn last_error = nullptr;
};

struct Key {
    int suite = 1;
    int dim = 2;
    bool operator<(const Key& rhs) const {
        return suite < rhs.suite || (suite == rhs.suite && dim < rhs.dim);
    }
};

static Api g_api;
static std::map<Key, void*> g_handles;

static std::string dll_error() {
    if (g_api.last_error) {
        const char* msg = g_api.last_error();
        if (msg && *msg) return msg;
    }
    return "OFEC MPMMO DLL call failed";
}

static void cleanup() {
    for (auto& kv : g_handles) {
        if (kv.second && g_api.destroy) g_api.destroy(kv.second);
    }
    g_handles.clear();
    if (g_api.dll) {
        FreeLibrary(g_api.dll);
        g_api = Api{};
    }
}

static FARPROC must_get(const char* name) {
    FARPROC ptr = GetProcAddress(g_api.dll, name);
    if (!ptr) mexErrMsgIdAndTxt("OFEC:MPMMO:MissingSymbol", "Missing DLL symbol: %s", name);
    return ptr;
}

static void load_api() {
    if (g_api.dll) return;
    g_api.dll = LoadLibraryA("ofec_mpmmo.dll");
    if (!g_api.dll) {
        mexErrMsgIdAndTxt("OFEC:MPMMO:LoadDLL",
            "Cannot load ofec_mpmmo.dll. Put it next to this MEX file or add its folder to PATH.");
    }
    g_api.create = reinterpret_cast<CreateFn>(must_get("ofec_mpmmo_create"));
    g_api.destroy = reinterpret_cast<DestroyFn>(must_get("ofec_mpmmo_destroy"));
    g_api.evaluate = reinterpret_cast<EvaluateFn>(must_get("ofec_mpmmo_evaluate"));
    g_api.count = reinterpret_cast<CountFn>(must_get("ofec_mpmmo_shared_count"));
    g_api.shared = reinterpret_cast<SharedFn>(must_get("ofec_mpmmo_shared_optima"));
    g_api.last_error = reinterpret_cast<ErrorFn>(must_get("ofec_mpmmo_last_error"));
    mexAtExit(cleanup);
}

static void* get_handle(int suite, int dim) {
    load_api();
    Key key{ suite, dim < 2 ? 2 : dim };
    auto it = g_handles.find(key);
    if (it != g_handles.end()) return it->second;
    void* h = g_api.create(key.suite, key.dim);
    if (!h) mexErrMsgIdAndTxt("OFEC:MPMMO:Create", "%s", dll_error().c_str());
    g_handles[key] = h;
    return h;
}

static std::string get_command(const mxArray* arg) {
    if (!mxIsChar(arg)) mexErrMsgIdAndTxt("OFEC:MPMMO:Input", "First input must be a command string.");
    char* raw = mxArrayToString(arg);
    std::string cmd(raw ? raw : "");
    mxFree(raw);
    return cmd;
}

static int scalar_int(const mxArray* arg, const char* name) {
    if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
        mexErrMsgIdAndTxt("OFEC:MPMMO:Input", "%s must be a real scalar.", name);
    return static_cast<int>(mxGetScalar(arg));
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs < 1) mexErrMsgIdAndTxt("OFEC:MPMMO:Input", "Usage: ofec_mpmmo_mex(command, ...)");
    const std::string cmd = get_command(prhs[0]);

    if (cmd == "clear") {
        cleanup();
        return;
    }

    if (cmd == "evaluate" || cmd == "eval") {
        if (nrhs < 4) mexErrMsgIdAndTxt("OFEC:MPMMO:Input", "evaluate requires suite, dimension, and PopDec.");
        const int suite = scalar_int(prhs[1], "suite");
        const int dim = scalar_int(prhs[2], "dimension");
        const mxArray* X = prhs[3];
        if (!mxIsDouble(X) || mxIsComplex(X) || mxGetN(X) != static_cast<size_t>(dim))
            mexErrMsgIdAndTxt("OFEC:MPMMO:Input", "PopDec must be an N-by-D real double matrix.");
        const mwSize n = mxGetM(X);
        const double* x_col_major = mxGetPr(X);
        std::vector<double> x_row_major(static_cast<size_t>(n) * dim);
        for (mwSize i = 0; i < n; ++i)
            for (int d = 0; d < dim; ++d)
                x_row_major[static_cast<size_t>(i) * dim + d] = x_col_major[i + static_cast<mwSize>(d) * n];

        std::vector<double> obj_row_major(static_cast<size_t>(n) * 2);
        if (!g_api.evaluate(get_handle(suite, dim), x_row_major.data(), static_cast<int>(n), obj_row_major.data()))
            mexErrMsgIdAndTxt("OFEC:MPMMO:Evaluate", "%s", dll_error().c_str());

        plhs[0] = mxCreateDoubleMatrix(n, 2, mxREAL);
        double* out = mxGetPr(plhs[0]);
        for (mwSize i = 0; i < n; ++i) {
            out[i] = obj_row_major[static_cast<size_t>(i) * 2 + 0];
            out[i + n] = obj_row_major[static_cast<size_t>(i) * 2 + 1];
        }
        return;
    }

    if (cmd == "shared") {
        if (nrhs < 3) mexErrMsgIdAndTxt("OFEC:MPMMO:Input", "shared requires suite and dimension.");
        const int suite = scalar_int(prhs[1], "suite");
        const int dim = scalar_int(prhs[2], "dimension");
        void* h = get_handle(suite, dim);
        const int count = g_api.count(h);
        if (count < 0) mexErrMsgIdAndTxt("OFEC:MPMMO:Shared", "%s", dll_error().c_str());
        std::vector<double> pts_row_major(static_cast<size_t>(count) * dim);
        const int got = g_api.shared(h, pts_row_major.data());
        if (got < 0) mexErrMsgIdAndTxt("OFEC:MPMMO:Shared", "%s", dll_error().c_str());
        plhs[0] = mxCreateDoubleMatrix(count, dim, mxREAL);
        double* out = mxGetPr(plhs[0]);
        for (int i = 0; i < count; ++i)
            for (int d = 0; d < dim; ++d)
                out[i + static_cast<mwSize>(d) * count] = pts_row_major[static_cast<size_t>(i) * dim + d];
        return;
    }

    mexErrMsgIdAndTxt("OFEC:MPMMO:Command", "Unknown command: %s", cmd.c_str());
}
