#include "../../utils/dataloader.h"
#include "../../utils/lib.h"

#define w_order 2

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
        /// TODO: Increase the acc_capacity
        return;
    }

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
            lastmode_set = last_mode.w_index_set[lastmode_bucket_id]; // error
            if (lastmode_set) {
                lastmode_point_id = last_mode.w_index_list[lastmode_bucket_id]; // error
                acc.val[lastmode_point_id] += acc.val[point_id]; // Reduce the value with the same position.
            }
            else {
                last_mode.w_index_list[lastmode_bucket_id] = point_id; // error
                last_mode.w_index_set[lastmode_bucket_id] = true; // error
                last_mode.w_index_set_ids[last_mode.w_index_num] = lastmode_bucket_id;
                last_mode.w_index_num++;
            }
        }
        // Clear first mode bucket
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
            last_mode.w_index_set[lastmode_bucket_id] = false; // Clear last mode bucket // error
        }
        last_mode.w_index_num = 0;

    }
    // Clear first mode bucket
    first_mode.w_index_num = 0;
    // Merge
    int acc_size = acc_index_id;
    if (COO_size == 0) {
        if (acc_size > COO_cap) {
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
        return;
    }
    int tmp_cap = COO_size + acc_size;
    int32_t * tmp_COO_crd[2];
    float* tmp_COO_vals;
    tmp_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * tmp_cap);
    tmp_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * tmp_cap);
    tmp_COO_vals = (float*)malloc(sizeof(float) * tmp_cap);
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
    COO_size = target_pointer;
}

int compute(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C, taco_tensor_t *D, int32_t w_accumulator_capacity) {
  int A1_dimension = (int)(A->dimensions[0]);
  int A2_dimension = (int)(A->dimensions[1]);
  int* restrict B1_pos = (int*)(B->indices[0][0]);
  int* restrict B1_crd = (int*)(B->indices[0][1]);
  int* restrict B2_pos = (int*)(B->indices[1][0]);
  int* restrict B2_crd = (int*)(B->indices[1][1]);
  int* restrict B3_pos = (int*)(B->indices[2][0]);
  int* restrict B3_crd = (int*)(B->indices[2][1]);
  float* restrict B_vals = (float*)(B->vals);
  int* restrict C1_pos = (int*)(C->indices[0][0]);
  int* restrict C1_crd = (int*)(C->indices[0][1]);
  int* restrict C2_pos = (int*)(C->indices[1][0]);
  int* restrict C2_crd = (int*)(C->indices[1][1]);
  float* restrict C_vals = (float*)(C->vals);
  int* restrict D1_pos = (int*)(D->indices[0][0]);
  int* restrict D1_crd = (int*)(D->indices[0][1]);
  int* restrict D2_pos = (int*)(D->indices[1][0]);
  int* restrict D2_crd = (int*)(D->indices[1][1]);
  float* restrict D_vals = (float*)(D->vals);

  int32_t w_accumulator_size = 0;
  int32_t w_all_capacity = w_accumulator_capacity; 
  int32_t w_all_size = 0;
  int32_t* w1_crd = 0;
  int32_t* w2_crd = 0;
  float* w_vals = 0;
  w1_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  w2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  w_vals = (float*)malloc(sizeof(float) * w_all_capacity);
  bool w_insertFail = false; 
  int32_t w_point[w_order];
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
  int32_t w_firstmode_size = A1_dimension;
  int32_t w_lastmode_size = A2_dimension;
  InitFirstMode(w_first_mode, w_firstmode_size);
  InitLastMode(w_last_mode, w_lastmode_size);
  int32_t w_acc_size = 0;
  float val;

  int32_t kB = B1_pos[0];
  int32_t pB1_end = B1_pos[1];
  int32_t kC = C1_pos[0];
  int32_t pC1_end = C1_pos[1];
  //printf("kB: %d, pB1_end: %d, kC: %d, pC1_end: %d\n", kB, pB1_end, kC, pC1_end);

  while (kB < pB1_end && kC < pC1_end) {
    int32_t kB0 = B1_crd[kB];
    int32_t kC0 = C1_crd[kC];
    int32_t k = TACO_MIN(kB0,kC0);
    if (kB0 == k && kC0 == k) {
      int32_t lB = B2_pos[kB];
      int32_t pB2_end = B2_pos[(kB + 1)];
      int32_t lD = D1_pos[0];
      int32_t pD1_end = D1_pos[1];

      while (lB < pB2_end && lD < pD1_end) {
        int32_t lB0 = B2_crd[lB];
        int32_t lD0 = D1_crd[lD];
        int32_t l = TACO_MIN(lB0,lD0);
        if (lB0 == l && lD0 == l) {
          for (int32_t iB = B3_pos[lB]; iB < B3_pos[(lB + 1)]; iB++) {
            int32_t i = B3_crd[iB];
            w_point[0] = i;
            int32_t jC = C2_pos[kC];
            int32_t pC2_end = C2_pos[(kC + 1)];
            int32_t jD = D2_pos[lD];
            int32_t pD2_end = D2_pos[(lD + 1)];

            while (jC < pC2_end && jD < pD2_end) {
              int32_t jC0 = C2_crd[jC];
              int32_t jD0 = D2_crd[jD];
              int32_t j = TACO_MIN(jC0,jD0);
              
              if (jC0 == j && jD0 == j) {
                w_point[1] = j;
                val = B_vals[iB] * C_vals[jC] * D_vals[jD];
                // std::cout << "TryInsert: " << w_point[0] << " , " << w_point[1] << std::endl;
                // TryInsert_hash(w_insertFail, &w_accumulator, w_point, (B_vals[iB] * C_vals[jC] * D_vals[jD]), A2_dimension);
                TryInsert_coord_bucket(w_insertFail, w_acc_size, w_accumulator_capacity, w_first_mode, w_point, val, w_acc);

                if (w_insertFail) {
                  // Enlarge in Merge function to avoid unnecessary memory allocation
                  SortMerge_coord_bucket(&w1_crd, &w2_crd, &w_vals, w_all_size, w_all_capacity, w_acc, w_first_mode, w_last_mode, w_acc_index);
                  w_acc_size = 0;
                  TryInsert_coord_bucket(w_insertFail, w_acc_size, w_accumulator_capacity, w_first_mode, w_point, val, w_acc);
                }
              }
              jC += (int32_t)(jC0 == j);
              jD += (int32_t)(jD0 == j);
            }
          }
        }
        lB += (int32_t)(lB0 == l);
        lD += (int32_t)(lD0 == l);
      }
    }
    kB += (int32_t)(kB0 == k);
    kC += (int32_t)(kC0 == k);
  }
  if (w_acc_size > 0) {
      SortMerge_coord_bucket(&w1_crd, &w2_crd, &w_vals, w_all_size, w_all_capacity, w_acc, w_first_mode, w_last_mode, w_acc_index);
      w_acc_size = 0;
  }

  A->indices[1][0] = (int32_t*)(w1_crd); 
  A->indices[2][0] = (int32_t*)(w2_crd); // It should be A->indices[2][0]. However, it works because we initial the CSR and CSR has indices[1][1]. Therefore, don't change the init code!!!
  A->vals = (float*)(w_vals);
  A->vals_size = w_all_size;

  free(w_acc.crd[0]);
  free(w_acc.crd[1]);
  free(w_acc.val);
  free(w_acc_index);
  FreeFristMode(w_first_mode, w_firstmode_size);
  FreeLastMode(w_last_mode, w_lastmode_size);

  return 0;
}

double COO_CSF_DCSR_DCSR_hash(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, taco_tensor_t* D, int32_t w_cap, int32_t warmup, int32_t repeat, bool bench = false, bool print = false) {
  for (int i = 0; i < warmup; i++) {
    compute(A,B,C,D,w_cap);
    if (bench) {
      free(A->vals);
      free(A->indices[1][0]);
      free(A->indices[2][0]);
    }
  }
  double start = clock();
  for (int i = 0; i < repeat; i++) {
    compute(A,B,C,D,w_cap);
    if (bench && i != repeat - 1) {
      free(A->vals);
      free(A->indices[1][0]);
      free(A->indices[2][0]);
    }
  }
  double end = clock();
  double duration = (double)(end - start) / (CLOCKS_PER_SEC * repeat);
  if (print) {
    std::cout << "A: " << std::endl;
    print_taco_tensor_COO(A);
  }
  return duration;
}