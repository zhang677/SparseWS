#include "../../utils/dataloader.h"
#include "../../utils/lib.h"
#include <math.h> 

void run_mttkrp_splatt(splatt_csf *B, splatt_kruskal *mats) { 
    double *cpd_opts = splatt_default_opts();
    cpd_opts[SPLATT_OPTION_NTHREADS] = 1;
    cpd_opts[SPLATT_OPTION_NITER] = 1;
    cpd_opts[SPLATT_OPTION_CSF_ALLOC] = SPLATT_CSF_ONEMODE;
    cpd_opts[SPLATT_OPTION_TILE] = SPLATT_NOTILE;

    const int mode = B->dim_perm[0];
    splatt_mttkrp(mode, mats->rank, B, mats->factors, mats->factors[mode], cpd_opts);
    splatt_free_opts(cpd_opts);
}

double mttkrp_splatt(const std::string& tensor_path, int J, int warmup, int repeat) {
    double *cpd_opts = splatt_default_opts();
    cpd_opts[SPLATT_OPTION_NTHREADS] = omp_get_num_threads();
    cpd_opts[SPLATT_OPTION_NITER] = 0;
    cpd_opts[SPLATT_OPTION_CSF_ALLOC] = SPLATT_CSF_ONEMODE;
    cpd_opts[SPLATT_OPTION_TILE] = SPLATT_NOTILE;
    cpd_opts[SPLATT_OPTION_VERBOSITY] = SPLATT_VERBOSITY_NONE;

    splatt_idx_t nmodes;
    splatt_csf *tensor;
    splatt_csf_load(tensor_path.c_str(), &nmodes, &tensor, cpd_opts);
    splatt_free_opts(cpd_opts);
    std::cout << "Finish input" << std::endl;
    splatt_kruskal factor_matrices = gen_factor_matrices(J, tensor);
    std::cout << "Finish dense generation" << std::endl;
    for (int i = 0; i < warmup; ++i) {
        run_mttkrp_splatt(tensor, &factor_matrices);
    }
    Timer timer;
    timer.reset();
    for (int i = 0; i < repeat; ++i) {
        run_mttkrp_splatt(tensor, &factor_matrices);
    }
    return timer.elapsed();
}