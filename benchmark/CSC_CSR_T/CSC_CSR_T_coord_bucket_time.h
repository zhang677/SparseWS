#include "../../utils/dataloader.h"
#include "../../utils/lib.h"

#define w_order 2

Timer timer;
std::ofstream outfile;
int num_tryinsert;
int num_alloc_tryinsert;
int num_free_coord_merge;
int curr_alloc_mem;
int num_merge;

struct wspace_t {
  int32_t* crd[w_order];
  float* val;
};

struct w_mode_t {
    bool* w_index_set;
    int32_t* w_index_list_size;
    int32_t* w_index_list_capacity;
    int32_t** w_index_list;
    int32_t* w_index_set_ids;
    int32_t w_index_num;
};

struct w_last_mode_t {
    bool* w_index_set;
    int32_t* w_index_list;
    int32_t* w_index_set_ids;
    int32_t w_index_num;
};

void InitFirstMode(w_mode_t& first_mode, int32_t first_mode_size) {
    first_mode.w_index_set = (bool*)calloc(first_mode_size, sizeof(bool));
    first_mode.w_index_list = (int32_t**)malloc(sizeof(int32_t*) * first_mode_size);
    first_mode.w_index_list_size = (int32_t*)calloc(first_mode_size, sizeof(int32_t));
    first_mode.w_index_list_capacity = (int32_t*)calloc(first_mode_size, sizeof(int32_t));
    first_mode.w_index_set_ids = (int32_t*)calloc(first_mode_size, sizeof(int32_t));
    first_mode.w_index_num = 0;
}

void FreeFristMode(w_mode_t& first_mode, int32_t first_mode_size) {
    first_mode.w_index_num = 0;
    free(first_mode.w_index_set);
    free(first_mode.w_index_list);
    free(first_mode.w_index_list_size);
    free(first_mode.w_index_list_capacity);
    free(first_mode.w_index_set_ids);
}

// void InitModes(w_mode_t* rest_modes, int32_t* rest_modes_size) {
//     for (int i = 0; i < w_order - 1; i++) {
//         rest_modes[i].w_index_set = (bool*)malloc(sizeof(bool) * rest_modes_size[i]);
//         rest_modes[i].w_index_list = (int32_t**)malloc(sizeof(int32_t*) * rest_modes_size[i]);
//         rest_modes[i].w_index_list_size = (int32_t*)calloc(sizeof(int32_t) * rest_modes_size[i]);
//         rest_modes[i].w_index_list_capacity = (int32_t*)calloc(sizeof(int32_t) * rest_modes_size[i]);
//     }
// }

void InitLastMode(w_last_mode_t& last_mode, int32_t last_mode_size) {
    last_mode.w_index_set = (bool*)calloc(last_mode_size, sizeof(bool));
    last_mode.w_index_list = (int32_t*)malloc(sizeof(int32_t) * last_mode_size);
    last_mode.w_index_set_ids = (int32_t*)calloc(last_mode_size, sizeof(int32_t));
    last_mode.w_index_num = 0;
}

void FreeLastMode(w_last_mode_t& last_mode, int32_t last_mode_size) {
    last_mode.w_index_num = 0;
    free(last_mode.w_index_set);
    free(last_mode.w_index_list);
    free(last_mode.w_index_set_ids);
}

int compare(wspace_t& acc, int32_t l, int32_t r, int32_t* __restrict__ crd0, int32_t* __restrict__  crd1) {
  if (acc.crd[0][l] == crd0[r]) {
    if (acc.crd[1][l] == crd1[r]) {
      return 0;
    } else if (acc.crd[1][l] < crd1[r]) {
      return -1;
    } else {
      return 1;
    }
  } else if (acc.crd[0][l] < crd0[r]) {
    return -1;
  } else {
    return 1;
  }

}

int compare_int(const void* numA, const void* numB)
{
    const int* num1 = (const int*)numA;
    const int* num2 = (const int*)numB;

    if (*num1 > *num2) {
        return 1;
    }
    else {
        if (*num1 == *num2)
            return 0;
        else
            return -1;
    }
}

void TryInsert_coord_bucket(bool& insertFail, int32_t& acc_size, int32_t acc_capacity, w_mode_t& first_mode, int32_t* point, float val, wspace_t& acc){
    if (acc_size == acc_capacity) {
        insertFail = true;
        /// TODO: Increase the acc_capacity (need to realloc the acc)
        return;
    }

    num_tryinsert++;
    bool do_alloc=false;

    insertFail = false;
    int32_t point_id = acc_size;
    acc.crd[0][point_id] = point[0];
    acc.crd[1][point_id] = point[1];
    acc.val[point_id] = val;
    acc_size += 1;

    int32_t bucket_id = point[0];
    int32_t tsize = first_mode.w_index_list_size[bucket_id];
    int32_t tcap;
    if (first_mode.w_index_set[bucket_id]) {
        tcap = first_mode.w_index_list_capacity[bucket_id];
        if (tsize == tcap) {
            tcap *= 2;
            first_mode.w_index_list[bucket_id] = (int32_t*)realloc(first_mode.w_index_list[bucket_id], sizeof(int32_t) * tcap); // error
            first_mode.w_index_list_capacity[bucket_id] = tcap;
            curr_alloc_mem += (tcap / 2) * 4;
            num_alloc_tryinsert ++;
            do_alloc = true;
        }
        first_mode.w_index_list[bucket_id][tsize] = point_id;
        tsize += 1;
        first_mode.w_index_list_size[bucket_id] = tsize;
    } 
    else {
        tcap = 1;
        first_mode.w_index_list[bucket_id] = (int32_t*)malloc(sizeof(int32_t) * tcap);
        first_mode.w_index_list_capacity[bucket_id] = tcap;
        first_mode.w_index_list[bucket_id][tsize] = point_id;
        tsize += 1;
        first_mode.w_index_list_size[bucket_id] = tsize;
        first_mode.w_index_set[bucket_id] = true;
        first_mode.w_index_set_ids[first_mode.w_index_num] = bucket_id;
        first_mode.w_index_num++;
        curr_alloc_mem += 1 * 4;
        num_alloc_tryinsert ++;
        do_alloc = true;
    }
    if (do_alloc) {
        outfile << timer.elapsed() << "," << curr_alloc_mem << "," << "T" << num_tryinsert << "+A" << num_alloc_tryinsert << std::endl; // T for TryInsert; A for Alloc
    } else {
        outfile << timer.elapsed() << "," << curr_alloc_mem << "," << "T" << num_tryinsert << std::endl;
    }
    return;
}

void SortMerge_coord_bucket(int32_t* __restrict__ * COO1_crd_ref, int32_t* __restrict__ * COO2_crd_ref, float* __restrict__ * COO_vals_ref, int32_t& COO_size, int32_t& COO_cap, wspace_t& acc, w_mode_t& first_mode, w_last_mode_t& last_mode, int32_t* acc_index) {
    int acc_index_id = 0; 
    int id;
    int point_id;
    int bucket_id;
    int list_id;
    int list_size;

    bool lastmode_set;
    int lastmode_id;
    int lastmode_bucket_id;
    int lastmode_point_id;
    
    num_merge++;
    // Seperate sort and merge now. Fuse them together to improve ILP in the future.
    // Sort
    qsort(first_mode.w_index_set_ids, first_mode.w_index_num, sizeof(int32_t), compare_int);
    for (id = 0; id < first_mode.w_index_num; id++) {
        bucket_id = first_mode.w_index_set_ids[id];
        list_size = first_mode.w_index_list_size[bucket_id];
        // Bucket sort the last mode
        for (list_id = 0; list_id < list_size; list_id++) {
            point_id = first_mode.w_index_list[bucket_id][list_id];
            lastmode_bucket_id = acc.crd[1][point_id];
            lastmode_set = last_mode.w_index_set[lastmode_bucket_id]; 
            if (lastmode_set) {
                lastmode_point_id = last_mode.w_index_list[lastmode_bucket_id]; 
                acc.val[lastmode_point_id] += acc.val[point_id]; // Reduce the value with the same position.
            }
            else {
                last_mode.w_index_list[lastmode_bucket_id] = point_id; 
                last_mode.w_index_set[lastmode_bucket_id] = true; 
                last_mode.w_index_set_ids[last_mode.w_index_num] = lastmode_bucket_id;
                last_mode.w_index_num++;
            }
        }
        // Clear first mode bucket
        curr_alloc_mem -= first_mode.w_index_list_capacity[bucket_id] * 4;
        outfile << timer.elapsed() << "," << curr_alloc_mem << "," << "FC_" << id << "_" << bucket_id << "_" << last_mode.w_index_num << "/" << first_mode.w_index_list_capacity[bucket_id] << std::endl; // FC is short for freeing acc array list for the first order; the first number means the id; the second number means the coord; the final number shows the reduction proportion. 
        first_mode.w_index_set[bucket_id] = false;
        first_mode.w_index_list_size[bucket_id] = 0;
        free(first_mode.w_index_list[bucket_id]);
        first_mode.w_index_list[bucket_id] = NULL;
        first_mode.w_index_list_capacity[bucket_id] = 0;
        // Fill the acc_index
        qsort(last_mode.w_index_set_ids, last_mode.w_index_num, sizeof(int32_t), compare_int);
        for (lastmode_id = 0; lastmode_id < last_mode.w_index_num; lastmode_id++) {
            lastmode_bucket_id = last_mode.w_index_set_ids[lastmode_id];
            point_id = last_mode.w_index_list[lastmode_bucket_id];
            acc_index[acc_index_id] = point_id;
            acc_index_id++;
            last_mode.w_index_set[lastmode_bucket_id] = false; // Clear last mode bucket
        }
        last_mode.w_index_num = 0;
    }
    // Clear first mode bucket
    first_mode.w_index_num = 0;
    // Merge
    int acc_size = acc_index_id;
    if (COO_size == 0) {
        if (acc_size > COO_cap) {
            curr_alloc_mem += (acc_size - COO_cap) * 3 * 4;
            COO_cap = acc_size;
            (*COO1_crd_ref) = (int32_t*)realloc((*COO1_crd_ref), sizeof(int32_t) * COO_cap);
            (*COO2_crd_ref) = (int32_t*)realloc((*COO2_crd_ref), sizeof(int32_t) * COO_cap);
            (*COO_vals_ref) = (float*)realloc((*COO_vals_ref), sizeof(float) * COO_cap);
        }
        for (int i = 0; i < acc_size; i++) {
            (*COO1_crd_ref)[i] = acc.crd[0][acc_index[i]];
            (*COO2_crd_ref)[i] = acc.crd[1][acc_index[i]];
            (*COO_vals_ref)[i] = acc.val[acc_index[i]];
        }
        COO_size = acc_size;
        outfile << timer.elapsed() << "," << curr_alloc_mem << "," << "FM" << std::endl;
        return;
    }
    int tmp_cap = COO_size + acc_size;
    int32_t * tmp_COO_crd[2];
    float* tmp_COO_vals;
    tmp_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * tmp_cap);
    tmp_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * tmp_cap);
    tmp_COO_vals = (float*)malloc(sizeof(float) * tmp_cap);
    curr_alloc_mem += tmp_cap * 3 * 4;
    outfile << timer.elapsed() << "," << curr_alloc_mem << "," << "MA" << num_merge << std::endl; // "M" for Merge; "A" for Alloc
    int accumulator_pointer = 0;
    int content_pointer = 0;
    int target_pointer = 0;
    while(accumulator_pointer < acc_size && content_pointer < COO_size) {
        int cmp = compare(acc, acc_index[accumulator_pointer], content_pointer, (*COO1_crd_ref), (*COO2_crd_ref));
        if (cmp == 0) {
        tmp_COO_crd[0][target_pointer] = acc.crd[0][acc_index[accumulator_pointer]];
        tmp_COO_crd[1][target_pointer] = acc.crd[1][acc_index[accumulator_pointer]];
        tmp_COO_vals[target_pointer] = acc.val[acc_index[accumulator_pointer]] + (*COO_vals_ref)[content_pointer];
        accumulator_pointer ++;
        content_pointer ++;
        target_pointer ++;
        } else if (cmp < 0) {
        tmp_COO_crd[0][target_pointer] = acc.crd[0][acc_index[accumulator_pointer]];
        tmp_COO_crd[1][target_pointer] = acc.crd[1][acc_index[accumulator_pointer]];
        tmp_COO_vals[target_pointer] = acc.val[acc_index[accumulator_pointer]];
        accumulator_pointer ++;
        target_pointer ++;
        } else {
        tmp_COO_crd[0][target_pointer] = (*COO1_crd_ref)[content_pointer];
        tmp_COO_crd[1][target_pointer] = (*COO2_crd_ref)[content_pointer];
        tmp_COO_vals[target_pointer] = (*COO_vals_ref)[content_pointer];
        content_pointer ++;
        target_pointer ++;
        }
    }
    while(accumulator_pointer < acc_size) {
        tmp_COO_crd[0][target_pointer] = acc.crd[0][acc_index[accumulator_pointer]];
        tmp_COO_crd[1][target_pointer] = acc.crd[1][acc_index[accumulator_pointer]];
        tmp_COO_vals[target_pointer] = acc.val[acc_index[accumulator_pointer]];
        accumulator_pointer ++;
        target_pointer ++;
    }
    while(content_pointer < COO_size) {
        tmp_COO_crd[0][target_pointer] = (*COO1_crd_ref)[content_pointer];
        tmp_COO_crd[1][target_pointer] = (*COO2_crd_ref)[content_pointer];
        tmp_COO_vals[target_pointer] = (*COO_vals_ref)[content_pointer];
        content_pointer ++;
        target_pointer ++;
    }
    // Realloc and Copy
    if (target_pointer > COO_cap) {
        curr_alloc_mem += (target_pointer - COO_cap) * 3 * 4;
        outfile << timer.elapsed() << "," << curr_alloc_mem << "," << "MR" << num_merge << std::endl; // "M" for Merge; "R" for Realloc
        COO_cap = target_pointer;
        (*COO1_crd_ref) = (int32_t*)realloc((*COO1_crd_ref), sizeof(int32_t) * COO_cap);
        (*COO2_crd_ref) = (int32_t*)realloc((*COO2_crd_ref), sizeof(int32_t) * COO_cap);
        (*COO_vals_ref) = (float*)realloc((*COO_vals_ref), sizeof(float) * COO_cap);
    }
    for (int i = 0; i < target_pointer; i++) {
        (*COO1_crd_ref)[i] = tmp_COO_crd[0][i];
        (*COO2_crd_ref)[i] = tmp_COO_crd[1][i];
        (*COO_vals_ref)[i] = tmp_COO_vals[i];
    }
    free(tmp_COO_crd[0]);
    free(tmp_COO_crd[1]);
    free(tmp_COO_vals);
    curr_alloc_mem -= tmp_cap * 3 * 4;
    outfile << timer.elapsed() << "," << curr_alloc_mem << "," << "M" << num_merge << std::endl; // "M" for Merge
    COO_size = target_pointer;
}


int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B, int w_accumulator_capacity){
    int C1_dimension = (int)(C->dimensions[0]);
    int C2_dimension = (int)(C->dimensions[1]);
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

    timer.reset();
    
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
    int32_t w_all_capacity = w_accumulator_capacity; 
    int32_t w_all_size = 0;
    w1_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
    w2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity); 
    w_vals = (float*)malloc(sizeof(float) * w_all_capacity);
    bool w_insertFail = 0;
    int32_t w_point[2];
    wspace_t w_acc;
    w_acc.crd[0] = (int32_t*)malloc(sizeof(int32_t) * w_accumulator_capacity); 
    w_acc.crd[1] = (int32_t*)malloc(sizeof(int32_t) * w_accumulator_capacity); 
    w_acc.val = (float*)malloc(sizeof(float) * w_accumulator_capacity); 
    int32_t* w_acc_index;
    w_acc_index = (int32_t*)malloc(sizeof(int32_t) * w_accumulator_capacity);
    w_mode_t w_first_mode;
    w_last_mode_t w_last_mode;
    // w_mode_t rest_modes[w_order - 1]; // assert w_order > 1
    // int32_t rest_modes_size[w_order - 1];
    int32_t w_firstmode_size = C1_dimension;
    int32_t w_lastmode_size = C2_dimension;
    InitFirstMode(w_first_mode, w_firstmode_size);
    InitLastMode(w_last_mode, w_lastmode_size);
    int32_t w_acc_size = 0;
    float val;

    curr_alloc_mem = ((C1_dimension+1) + (w_all_capacity * 3) + (w_accumulator_capacity * 4) + (w_firstmode_size * 5) + (w_lastmode_size * 3)) * 4;
    outfile << timer.elapsed() << "," << curr_alloc_mem << "," << "Init" << std::endl;
    for (int32_t j = 0; j < B1_dimension; j++) {
        for (int32_t iA = A2_pos[j]; iA < A2_pos[j+1]; iA++) {
            int32_t i = A2_crd[iA];
            w_point[1] = i;
            for (int32_t kB = B2_pos[j]; kB < B2_pos[j+1]; kB++) {
                int32_t k = B2_crd[kB];
                w_point[0] = k;
                val = A_vals[iA] * B_vals[kB];

                TryInsert_coord_bucket(w_insertFail, w_acc_size, w_accumulator_capacity, w_first_mode, w_point, val, w_acc);
                if (w_insertFail) {
                    // Enlarge in Merge function to avoid unnecessary memory allocation
                    SortMerge_coord_bucket(&w1_crd, &w2_crd, &w_vals, w_all_size, w_all_capacity, w_acc, w_first_mode, w_last_mode, w_acc_index);
                    w_acc_size = 0;
                    TryInsert_coord_bucket(w_insertFail, w_acc_size, w_accumulator_capacity, w_first_mode, w_point, val, w_acc);
                }
            }
        }
    }
    if (w_acc_size > 0) {
        SortMerge_coord_bucket(&w1_crd, &w2_crd, &w_vals, w_all_size, w_all_capacity, w_acc, w_first_mode, w_last_mode, w_acc_index);
        w_acc_size = 0;
    }

    free(w_acc.crd[0]);
    free(w_acc.crd[1]);
    free(w_acc.val);
    free(w_acc_index);
    FreeFristMode(w_first_mode, w_firstmode_size);
    FreeLastMode(w_last_mode, w_lastmode_size);

    curr_alloc_mem -= ((w_accumulator_capacity * 4) + (w_firstmode_size * 5) + (w_lastmode_size * 3)) * 4;
    outfile << timer.elapsed() << "," << curr_alloc_mem << "," << "FF" << std::endl; // "FF" for final free

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

    int32_t csC2 = 0;
    for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
        csC2 += C2_pos[pC2];
        C2_pos[pC2] = csC2;
    }

    C->indices[1][0] = (int32_t*)(C2_pos);
    C->indices[1][1] = (int32_t*)(w2_crd);
    C->vals = (float*)w_vals;
    curr_alloc_mem -= pw1_end * 4;
    free(w1_crd);
    free(w1_pos);
    outfile << timer.elapsed() << "," << curr_alloc_mem << "," << "OR" << std::endl; // "OR" for COO to CSR
    return 0;
}
// Use the same calling convension with hash methods. But it doesn't use hash to reduce points with identical coordiantes.
double CSC_CSR_T_hash(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t w_cap, int32_t repeat, const string& filename, bool print = false) {
  double start = clock();
  for (int i = 0; i < repeat; i++) {
    outfile.open(filename + "_" + std::to_string(i) + ".csv", std::ios_base::app);
    compute(C,A,B,w_cap);
    outfile.close();
  }
  double end = clock();
  double duration = (double)(end - start) / (CLOCKS_PER_SEC * repeat);
  if (print) {
    print_taco_tensor_DC(A);
    print_taco_tensor_DC(B);
    print_taco_tensor_DC(C);
  }
  return duration;
}