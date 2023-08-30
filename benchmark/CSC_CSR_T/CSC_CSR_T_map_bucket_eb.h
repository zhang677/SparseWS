#include "../../utils/dataloader.h"
#include "../../utils/lib.h"

#define w_order 1
#define PThreads 8

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

typedef std::map<index_t, float> map_bucket_t;

int32_t TryInsert_map_bucket(bool& insertFail, int32_t map_size, map_bucket_t* accumulator, bool* accumulator_init, int32_t bucket_id, index_t& point, float val) {
    if (!accumulator_init[bucket_id]) {
        accumulator[bucket_id] = map_bucket_t();
        accumulator_init[bucket_id] = true;
    }
    auto handle = accumulator[bucket_id].insert(std::pair<index_t, float>(point, val));
    if (!handle.second) {
        handle.first->second += val;
    } else {
        map_size++;
    }
    insertFail = false; // Always success
    return map_size;
}

int32_t Merge_map_bucket(int32_t* COO1_crd, int32_t* COO2_crd, float* COO_vals, int32_t bucket_size, map_bucket_t* accumulator, bool* accumulator_init) {
    int32_t pCOO = 0;
    float nz_row = 0;
    #pragma unroll
    for (int32_t bucket_id = 0; bucket_id < bucket_size; ++bucket_id) {
        if (accumulator_init[bucket_id]) {
            for (auto it = accumulator[bucket_id].begin(); it != accumulator[bucket_id].end(); ++it) {
                COO1_crd[pCOO] = bucket_id;
                COO2_crd[pCOO] = it->first.crd[0];
                COO_vals[pCOO] = it->second;
                pCOO++;
            }
            nz_row ++;
            accumulator[bucket_id].clear();
        }
    }
    // std::cout << "NZ Row: " << nz_row << "/" << bucket_size << " : " << nz_row / bucket_size << std::endl; // Mostly 1
    return pCOO;
}

int32_t Merge_map_bucket_parallel(int32_t* COO1_crd, int32_t* COO2_crd, float* COO_vals, int32_t elements_num, int32_t bucket_size, map_bucket_t* accumulator, bool* accumulator_init) {
    int32_t bucket_ids[PThreads];
    int32_t COO_ids[PThreads];
    // (accumulator[bucket_ids[i]], accumulator[bucket_ids[i+1]]) containes Total element/PThreads elements
    int32_t elements_accum = 0;
    int32_t barrier_id = 0;
    bucket_ids[0] = barrier_id;
    COO_ids[0] = 0;
    // Workload balance. Can't handle very imbalanced workload, for example, 1, 1, 1, 1000000
    for (int i = 0; i < bucket_size; i++) {
        if (accumulator_init[i]) {
            elements_accum += accumulator[i].size();
            if (elements_accum >= (barrier_id + 1) * (elements_num / PThreads)) {
                barrier_id++;
                bucket_ids[barrier_id] = i;
                COO_ids[barrier_id] = elements_accum - accumulator[i].size();
            }
            if (barrier_id == PThreads - 1) {
                break;
            }
        }
    }
    // Each thread copies elements from [bucket_ids[i], bucket_ids[i+1]) to COO
    // #pragma omp parallel for schedule(runtime)
    // for (int32_t i = 0; i < PThreads; i++) {
    //     int32_t pCOO = COO_ids[i];
    //     int32_t iterend;
    //     if (i == PThreads - 1) {
    //         iterend = bucket_size;
    //     } else {
    //         iterend = bucket_ids[i + 1];
    //     }
    //     for (int32_t bucket_id = bucket_ids[i]; bucket_id < iterend; ++bucket_id) {
    //         if (accumulator_init[bucket_id]) {
    //             for (auto it = accumulator[bucket_id].begin(); it != accumulator[bucket_id].end(); ++it) {
    //                 COO1_crd[pCOO] = bucket_id;
    //                 COO2_crd[pCOO] = it->first.crd[0];
    //                 COO_vals[pCOO] = it->second;
    //                 pCOO++;
    //             }
    //             accumulator[bucket_id].clear();
    //         }
    //     }
    // }

    // We have PThreads in total. Each thread copies elements from [bucket_ids[i], bucket_ids[i+1]) to COO
    // #pragma omp parallel num_threads(PThreads) 
    #pragma omp parallel for schedule(runtime)
    {
        int i = omp_get_thread_num();
        int32_t pCOO = COO_ids[i];
        int32_t iterend;
        if (i == PThreads - 1) {
            iterend = bucket_size;
        } else {
            iterend = bucket_ids[i + 1];
        }
        #pragma unroll
        for (int32_t bucket_id = bucket_ids[i]; bucket_id < iterend; ++bucket_id) {
            if (accumulator_init[bucket_id]) {
                for (auto it = accumulator[bucket_id].begin(); it != accumulator[bucket_id].end(); ++it) {
                    COO1_crd[pCOO] = bucket_id;
                    COO2_crd[pCOO] = it->first.crd[0];
                    COO_vals[pCOO] = it->second;
                    pCOO++;
                }
                accumulator[bucket_id].clear();
            }
        }
    }

    return elements_num;
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
    int32_t bucket_size = C1_dimension;
    int32_t bucket_id = 0;
    int32_t w_all_size = 0;

    map_bucket_t* w_accumulator = 0;
    w_accumulator = new map_bucket_t[bucket_size];
    bool* w_accumulator_init = 0;
    w_accumulator_init = (bool*)calloc(bucket_size, sizeof(bool));

    for (int32_t j = 0; j < B1_dimension; j++) {
        for (int32_t iA = A2_pos[j]; iA < A2_pos[(j + 1)]; iA++) {
            int32_t i = A2_crd[iA];
            w_point.crd[0] = i;
            for (int32_t kB = B2_pos[j]; kB < B2_pos[(j + 1)]; kB++) {
                int32_t k = B2_crd[kB];
                bucket_id = k;
                w_all_size = TryInsert_map_bucket(w_insertFail, w_all_size, w_accumulator, w_accumulator_init, bucket_id, w_point, A_vals[iA] * B_vals[kB]);
            }
        }
    }

    w1_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_size);
    w2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_size);
    w_vals = (float*)malloc(sizeof(float) * w_all_size);
    w_all_size = Merge_map_bucket_parallel(w1_crd, w2_crd, w_vals, w_all_size, bucket_size, w_accumulator, w_accumulator_init);

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
    delete[] w_accumulator;
    free(w_accumulator_init);
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