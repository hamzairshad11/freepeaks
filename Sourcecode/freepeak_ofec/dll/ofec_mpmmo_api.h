#ifndef OFEC_MPMMO_API_H
#define OFEC_MPMMO_API_H

#ifdef _WIN32
#  ifdef OFEC_MPMMO_EXPORTS
#    define OFEC_MPMMO_API __declspec(dllexport)
#  else
#    define OFEC_MPMMO_API __declspec(dllimport)
#  endif
#else
#  define OFEC_MPMMO_API
#endif

extern "C" {
    OFEC_MPMMO_API void* ofec_mpmmo_create(int suite_id, int dimension);
    OFEC_MPMMO_API void ofec_mpmmo_destroy(void* handle);
    OFEC_MPMMO_API int ofec_mpmmo_evaluate(void* handle, const double* decisions, int n_points, double* objectives);
    OFEC_MPMMO_API int ofec_mpmmo_shared_count(void* handle);
    OFEC_MPMMO_API int ofec_mpmmo_shared_optima(void* handle, double* out_points);
    OFEC_MPMMO_API const char* ofec_mpmmo_last_error();
}

#endif
