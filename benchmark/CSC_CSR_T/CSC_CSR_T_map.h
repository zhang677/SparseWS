#include "../../utils/dataloader.h"
#include "../../utils/lib.h"

#define w_order 2

struct index_t {
    int32_t crd[w_order];

    bool operator < (const index_t& other) const {
        for (int i = 0; i < w_order; ++i) {
            if (crd[i] != other.crd[i]) {
                return crd[i] < other.crd[i];
            }
        }
        return false;
    }
};

int32_t TryInsert_map(bool& insertFail, std::map<index_t, float>& accumulator, index_t& point, float val) {
    auto handle = accumulator.insert(std::pair<index_t, float>(point, val));
    if (!handle.second) {
        handle.first->second += val;
    }
    insertFail = false; // Always success
    return accumulator.size();
}

int32_t Merge_map(int32_t* COO1_crd, int32_t* COO2_crd, float* COO_vals, std::map<index_t, float>& accumulator) {
    int32_t pCOO = 0;
    for (auto it = accumulator.begin(); it != accumulator.end(); ++it) {
        COO1_crd[pCOO] = it->first.crd[0];
        COO2_crd[pCOO] = it->first.crd[1];
        COO_vals[pCOO] = it->second;
        pCOO++;
    }
    return pCOO;
}

int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B){
    int C1_dimension = (int)(C->dimensions[0]);
    int* restrict C2_pos = (int*)(C->indices[1][0]);
    int* restrict C2_crd = (int*)(C->indices[1][1]);
    float* restrict C_vals = (float*)(C->vals);
    int* restrict A2_pos = (int*)(A->indices[1][0]);
    int* restrict A2_crd = (int*)(A->indices[1][1]);
    float* restrict A_vals = (float*)(A->vals);
    int B1_dimension = (int)(B->dimensions[0]);
    int B2_dimension = (int)(B->dimensions[1]);
    int* restrict B2_pos = (int*)(B->indices[1][0]);
    int* restrict B2_crd = (int*)(B->indices[1][1]);
    float* restrict B_vals = (float*)(B->vals);

    C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
    C2_pos[0] = 0;
    for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
        C2_pos[pC2] = 0;
    }
    int32_t* restrict w1_pos = 0;
    w1_pos = (int32_t*)malloc(sizeof(int32_t) * 2);
    int32_t* restrict w1_crd = 0;
    int32_t* restrict w2_crd = 0;
    float* restrict w_vals = 0;
    bool w_insertFail = 0;
    index_t w_point;

    std::map<index_t, float> w_accumulator;
    for (int32_t j = 0; j < B1_dimension; j++) {
        for (int32_t iA = A2_pos[j]; iA < A2_pos[(j + 1)]; iA++) {
            int32_t i = A2_crd[iA];
            w_point.crd[1] = i;
            for (int32_t kB = B2_pos[j]; kB < B2_pos[(j + 1)]; kB++) {
                int32_t k = B2_crd[kB];
                w_point.crd[0] = k;
                TryInsert_map(w_insertFail, w_accumulator, w_point, A_vals[iA] * B_vals[kB]);
            }
        }
    }
    int32_t w_all_size = w_accumulator.size();
    w1_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_size);
    w2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_size);
    w_vals = (float*)malloc(sizeof(float) * w_all_size);
    w_all_size = Merge_map(w1_crd, w2_crd, w_vals, w_accumulator);

    w1_pos[0] = 0;
    w1_pos[1] = w_all_size;
    int32_t kw = w1_pos[0];
    int32_t pw1_end = w1_pos[1];

    while (kw < pw1_end) {
        int32_t k = w1_crd[kw];
        int32_t w1_segend = kw + 1;
        while (w1_segend < pw1_end && w1_crd[w1_segend] == k) {
        w1_segend++;
        }
        C2_pos[k + 1] = w1_segend - kw;
        kw = w1_segend;
    }
    w_accumulator.clear();
    free(w1_crd);
    free(w1_pos);
    
    int32_t csC2 = 0;
    for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
        csC2 += C2_pos[pC2];
        C2_pos[pC2] = csC2;
    }

    C->indices[1][0] = (int32_t*)(C2_pos);
    C->indices[1][1] = (int32_t*)(w2_crd);
    C->vals = (float*)w_vals;
    return 0;
}

double CSC_CSR_T_map(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t warmup, int32_t repeat, bool bench=false, bool print = false) {
    for (int i = 0; i < warmup; i++) {
        compute(C,A,B);
    }
    Timer t;
    t.reset();
    for (int i = 0; i < repeat; i++) {
        compute(C,A,B);
        if (bench && i != repeat - 1) {
            free(C->vals);
            free(C->indices[1][0]);
            free(C->indices[1][1]);
        }
    }
    double duration = t.elapsed() / repeat;
    if (print) {
        print_taco_tensor_DC(A);
        print_taco_tensor_DC(B);
        print_taco_tensor_DC(C);
    }
    return duration;
}