#include "../../utils/dataloader.h"
#include <math.h> 

#define w_order 3

int acc_capacity;
struct wspace_t {
  int32_t* crd[w_order];
  float* val;
};
wspace_t w_accumulator;
int32_t* w_accumulator_index;

int esc_cmp_t(const void* a, const void* b) {
  int l = *(int32_t *)a;
  int r = *(int32_t *)b;
  for (int i = 0; i < w_order; i++) {
    if (w_accumulator.crd[i][l] == w_accumulator.crd[i][r]) continue;
    return w_accumulator.crd[i][l] - w_accumulator.crd[i][r];
  }
  return w_accumulator.crd[w_order - 1][l] - w_accumulator.crd[w_order - 1][r];
}

int esc_cmp_t_rev(const void* a, const void* b) {
  int l = *(int32_t *)a;
  int r = *(int32_t *)b;
  for (int i = 0; i < w_order; i++) {
    if (w_accumulator.crd[i][r] == w_accumulator.crd[i][l]) continue;
    return w_accumulator.crd[i][r] - w_accumulator.crd[i][l];
  }
  return w_accumulator.crd[w_order - 1][r] - w_accumulator.crd[w_order - 1][l];
}

int Sort_t(size_t size, bool rev) {
  if (rev) {
    qsort(w_accumulator_index, size, sizeof(int32_t), esc_cmp_t_rev);
    for (int i = size - 1; i >= 0; i--) {
      if(w_accumulator.crd[0][w_accumulator_index[i]] != -1) {
        return i + 1;
      }
    }
  } 
  qsort(w_accumulator_index, size, sizeof(int32_t), esc_cmp_t);
  return size;
}

int compare(int32_t l, int32_t r, int32_t* crd0, int32_t* crd1, int32_t* crd2) {
  if (w_accumulator.crd[0][l] == crd0[r]) {
    if (w_accumulator.crd[1][l] == crd1[r]) {
      if (w_accumulator.crd[2][l] == crd2[r]) {
        return 0;
      } else if (w_accumulator.crd[2][l] < crd2[r]) {
        return -1;
      } else {
        return 1;
      }
    } else if (w_accumulator.crd[1][l] < crd1[r]) {
      return -1;
    } else {
      return 1;
    }
  } else if (w_accumulator.crd[0][l] < crd0[r]) {
    return -1;
  } else {
    return 1;
  }
}

int32_t TryInsert_coord_t(bool* insertFail, int32_t accumulator_size, int32_t* crds, float val) {
  if (accumulator_size == acc_capacity) {
    *insertFail = true;
    return Sort_t(acc_capacity, false);
  } else {
    w_accumulator_index[accumulator_size] = accumulator_size;
    w_accumulator.crd[0][accumulator_size] = crds[0];
    w_accumulator.crd[1][accumulator_size] = crds[1];
    w_accumulator.crd[2][accumulator_size] = crds[2];
    w_accumulator.val[accumulator_size] = val;
    *insertFail = false;
    accumulator_size ++;
    return accumulator_size;
  }
}

int32_t Merge_coord_t(int32_t* COO1_crd_0, int32_t* COO2_crd_0, int32_t* COO3_crd_0,  float* COO_vals_0,int32_t* COO1_crd_1, int32_t* COO2_crd_1, int32_t* COO3_crd_1, float* COO_vals_1, int32_t COO_size, int32_t accumulator_size) {
  int32_t* COO1_crd;
  int32_t* COO2_crd;
  int32_t* COO3_crd;
  float* COO_vals;
  COO1_crd = COO1_crd_0;
  COO2_crd = COO2_crd_0;
  COO3_crd = COO3_crd_0;
  COO_vals = COO_vals_0;
  int32_t* tmp_COO_crd[w_order];
  float* tmp_COO_vals;
  tmp_COO_crd[0] = COO1_crd_1;
  tmp_COO_crd[1] = COO2_crd_1;
  tmp_COO_crd[2] = COO3_crd_1;
  tmp_COO_vals = COO_vals_1;

  if (COO_size == 0) {
    tmp_COO_crd[0][0] = w_accumulator.crd[0][w_accumulator_index[0]];
    tmp_COO_crd[1][0] = w_accumulator.crd[1][w_accumulator_index[0]];
    tmp_COO_crd[2][0] = w_accumulator.crd[2][w_accumulator_index[0]];
    tmp_COO_vals[0] = w_accumulator.val[w_accumulator_index[0]];
    int target_pointer = 0;
    for (int i = 1; i < accumulator_size; i++) {
      if (compare(w_accumulator_index[i], w_accumulator_index[target_pointer], COO1_crd, COO2_crd, COO3_crd) == 0) {
        tmp_COO_vals[target_pointer] += w_accumulator.val[w_accumulator_index[i]];
      } else {
        target_pointer++;
        tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[i]];
        tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[i]];
        tmp_COO_crd[2][target_pointer] = w_accumulator.crd[2][w_accumulator_index[i]];
        tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[i]];
      }
    }
    return target_pointer + 1;
  }
  int accumulator_pointer = 0;
  int content_pointer = 0;
  int target_pointer = 0;
  int cmp = compare(w_accumulator_index[accumulator_pointer], content_pointer, COO1_crd, COO2_crd, COO3_crd);
  if (cmp == 0) {
    tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
    tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
    tmp_COO_crd[2][target_pointer] = w_accumulator.crd[2][w_accumulator_index[accumulator_pointer]];
    tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]] + COO_vals[content_pointer];
    accumulator_pointer ++;
    content_pointer ++;
  } else if (cmp < 0) {
    tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
    tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
    tmp_COO_crd[2][target_pointer] = w_accumulator.crd[2][w_accumulator_index[accumulator_pointer]];
    tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]];
    accumulator_pointer ++;
  } else {
    tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
    tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
    tmp_COO_crd[2][target_pointer] = COO3_crd[content_pointer];
    tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
    content_pointer ++;
  }

  while(accumulator_pointer < accumulator_size && content_pointer < COO_size) {
    cmp = compare(w_accumulator_index[accumulator_pointer], content_pointer, COO1_crd, COO2_crd, COO3_crd);
    if (cmp == 0) {
      if (compare(w_accumulator_index[accumulator_pointer], target_pointer, tmp_COO_crd[0], tmp_COO_crd[1], tmp_COO_crd[2]) == 0) {
        tmp_COO_vals[target_pointer] += w_accumulator.val[w_accumulator_index[accumulator_pointer]] + COO_vals[content_pointer];
      } else {
        target_pointer ++;
        tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
        tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
        tmp_COO_crd[2][target_pointer] = w_accumulator.crd[2][w_accumulator_index[accumulator_pointer]];
        tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]] + COO_vals[content_pointer];
      }
      accumulator_pointer ++;
      content_pointer ++;
    } else if (cmp < 0) {
      if (compare(w_accumulator_index[accumulator_pointer], target_pointer, tmp_COO_crd[0], tmp_COO_crd[1], tmp_COO_crd[2]) == 0) {
        tmp_COO_vals[target_pointer] += w_accumulator.val[w_accumulator_index[accumulator_pointer]];
      } else {
        target_pointer ++;
        tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
        tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
        tmp_COO_crd[2][target_pointer] = w_accumulator.crd[2][w_accumulator_index[accumulator_pointer]];
        tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]];
      }
      accumulator_pointer ++;
    } else {
      if(tmp_COO_crd[0][target_pointer] == COO1_crd[content_pointer] && tmp_COO_crd[1][target_pointer] == COO2_crd[content_pointer] && tmp_COO_crd[2][target_pointer] == COO3_crd[content_pointer]) {
        tmp_COO_vals[target_pointer] += COO_vals[content_pointer];
      } else {
        target_pointer ++;
        tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
        tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
        tmp_COO_crd[2][target_pointer] = COO3_crd[content_pointer];
        tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      }
      content_pointer ++;
    }
  }
  while(accumulator_pointer<accumulator_size) {
    if (compare(w_accumulator_index[accumulator_pointer], target_pointer, tmp_COO_crd[0], tmp_COO_crd[1], tmp_COO_crd[2]) == 0) {
      tmp_COO_vals[target_pointer] += w_accumulator.val[w_accumulator_index[accumulator_pointer]];
    } else {
      target_pointer ++;
      tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
      tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
      tmp_COO_crd[2][target_pointer] = w_accumulator.crd[2][w_accumulator_index[accumulator_pointer]];
      tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]];
    }
    accumulator_pointer ++;
  }
  while(content_pointer<COO_size) {
    if(tmp_COO_crd[0][target_pointer] == COO1_crd[content_pointer] && tmp_COO_crd[1][target_pointer] == COO2_crd[content_pointer] && tmp_COO_crd[2][target_pointer] == COO3_crd[content_pointer]) {
      tmp_COO_vals[target_pointer] += COO_vals[content_pointer];
    } else {
      target_pointer ++;
      tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
      tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
      tmp_COO_crd[2][target_pointer] = COO3_crd[content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
    }
    content_pointer ++;
  }
  return target_pointer + 1;
}

int compute_coord(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C, int32_t w_accumulator_capacity) {
  int A2_dimension = (int)(A->dimensions[1]);
  int A3_dimension = (int)(A->dimensions[2]);
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

  int32_t w_accumulator_size = 0;
  int32_t w_all_capacity_0 = acc_capacity; 
  int32_t w_all_capacity_1 = acc_capacity; 
  int32_t w_all_size = 0;
  int32_t* w1_crd_0 = 0;
  int32_t* w2_crd_0 = 0;
  int32_t* w3_crd_0 = 0;
  float* w_vals_0 = 0;
  w1_crd_0 = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity_0);
  w2_crd_0 = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity_0);
  w3_crd_0 = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity_0);
  w_vals_0 = (float*)malloc(sizeof(float) * w_all_capacity_0);
  int32_t* w1_crd_1 = 0;
  int32_t* w2_crd_1 = 0;
  int32_t* w3_crd_1 = 0;
  float* w_vals_1 = 0;
  w1_crd_1 = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity_1);
  w2_crd_1 = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity_1);
  w3_crd_1 = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity_1);
  w_vals_1 = (float*)malloc(sizeof(float) * w_all_capacity_1);
  int32_t* w1_crd = 0;
  int32_t* w2_crd = 0;
  int32_t* w3_crd = 0;
  float* w_vals = 0;
  int32_t all_array_id = 0;
  bool w_insertFail = false; 
  int32_t w_point[w_order];

  // int iterations = 0;
  // std::vector<int> w_all_caps;
  // w_all_caps.push_back(w_all_capacity);
  for (int32_t kB = B1_pos[0]; kB < B1_pos[1]; kB++) {
    int32_t k = B1_crd[kB];
    int32_t iB = B2_pos[kB];
    int32_t pB2_end = B2_pos[(kB + 1)];
    int32_t iC = C1_pos[0];
    int32_t pC1_end = C1_pos[1];
    w_point[0] = k;
    while (iB < pB2_end && iC < pC1_end) {
        int32_t iB0 = B2_crd[iB];
        int32_t iC0 = C1_crd[iC];
        int32_t i = TACO_MIN(iB0,iC0);
        if (iB0 == i && iC0 == i) {
            for (int32_t jB = B3_pos[iB]; jB < B3_pos[(iB + 1)]; jB++) {
                int32_t j = B3_crd[jB];
                int32_t jA = k * A2_dimension + j;
                w_point[1] = j;
                for (int32_t lC = C2_pos[iC]; lC < C2_pos[(iC + 1)]; lC++) {
                    int32_t l = C2_crd[lC];
                    int32_t lA = jA * A3_dimension + l;
                    w_point[2] = l;

                    w_accumulator_size = TryInsert_coord_t(&w_insertFail, w_accumulator_size, w_point, (B_vals[jB] * C_vals[lC]));
                    if (w_insertFail) { 
                      //counter += 1; 
                      if ((all_array_id == 0 && w_accumulator_size + w_all_size > w_all_capacity_1) ||
                          (all_array_id == 1 && w_accumulator_size + w_all_size > w_all_capacity_0)) {
                        //w_all_capacity = w_all_capacity * 2;
                        
                        if (all_array_id == 1) {
                          w_all_capacity_0 = w_accumulator_size + w_all_size;
                          w1_crd_0 = (int32_t*)realloc(w1_crd_0, sizeof(int32_t) * w_all_capacity_0);
                          w2_crd_0 = (int32_t*)realloc(w2_crd_0, sizeof(int32_t) * w_all_capacity_0);
                          w3_crd_0 = (int32_t*)realloc(w3_crd_0, sizeof(int32_t) * w_all_capacity_0);
                          w_vals_0 = (float*)realloc(w_vals_0, sizeof(float) * w_all_capacity_0);
                        } else {
                          w_all_capacity_1 = w_accumulator_size + w_all_size;
                          w1_crd_1 = (int32_t*)realloc(w1_crd_1, sizeof(int32_t) * w_all_capacity_1);
                          w2_crd_1 = (int32_t*)realloc(w2_crd_1, sizeof(int32_t) * w_all_capacity_1);
                          w3_crd_1 = (int32_t*)realloc(w3_crd_1, sizeof(int32_t) * w_all_capacity_1);
                          w_vals_1 = (float*)realloc(w_vals_1, sizeof(float) * w_all_capacity_1);
                        }
                      }
                      if (all_array_id == 0) {
                        w_all_size = Merge_coord_t(w1_crd_0, w2_crd_0, w3_crd_0, w_vals_0, w1_crd_1, w2_crd_1, w3_crd_1, w_vals_1, w_all_size, w_accumulator_size);
                      } else {
                        w_all_size = Merge_coord_t(w1_crd_1, w2_crd_1, w3_crd_1, w_vals_1, w1_crd_0, w2_crd_0, w3_crd_0, w_vals_0, w_all_size, w_accumulator_size);
                      }
                      all_array_id ^= 1;
                      // Heuristic. (Overflow is not considered)
                      if (acc_capacity < 4194304) {
                        acc_capacity *= 2;
                      } else  {
                        acc_capacity *= 1.5;
                      }
                      free(w_accumulator.crd[0]);
                      free(w_accumulator.crd[1]);
                      free(w_accumulator.crd[2]);
                      free(w_accumulator.val);
                      free(w_accumulator_index);
                      w_accumulator.crd[0] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
                      w_accumulator.crd[1] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
                      w_accumulator.crd[2] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
                      w_accumulator.val = (float*)malloc(sizeof(float) * acc_capacity);
                      w_accumulator_index = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
                      w_accumulator_size = 0;
                      w_accumulator_index[0] = 0;
                      w_accumulator_size = TryInsert_coord_t(&w_insertFail, w_accumulator_size, w_point, (B_vals[jB] * C_vals[lC]));
                    }
                }
            }
        }
        iB += (int32_t)(iB0 == i);
        iC += (int32_t)(iC0 == i);
    }
  }
  
  if (w_accumulator_size > 0) {
    // Sort
    w_accumulator_size = Sort_t(w_accumulator_size, false);
    // Merge
    if ((all_array_id == 0 && w_accumulator_size + w_all_size > w_all_capacity_1) ||
        (all_array_id == 1 && w_accumulator_size + w_all_size > w_all_capacity_0)) {
      //w_all_capacity = w_all_capacity * 2;
      
      if (all_array_id == 1) {
        w_all_capacity_0 = w_accumulator_size + w_all_size;
        w1_crd_0 = (int32_t*)realloc(w1_crd_0, sizeof(int32_t) * w_all_capacity_0);
        w2_crd_0 = (int32_t*)realloc(w2_crd_0, sizeof(int32_t) * w_all_capacity_0);
        w3_crd_0 = (int32_t*)realloc(w3_crd_0, sizeof(int32_t) * w_all_capacity_0);
        w_vals_0 = (float*)realloc(w_vals_0, sizeof(float) * w_all_capacity_0);
      } else {
        w_all_capacity_1 = w_accumulator_size + w_all_size;
        w1_crd_1 = (int32_t*)realloc(w1_crd_1, sizeof(int32_t) * w_all_capacity_1);
        w2_crd_1 = (int32_t*)realloc(w2_crd_1, sizeof(int32_t) * w_all_capacity_1);
        w3_crd_1 = (int32_t*)realloc(w3_crd_1, sizeof(int32_t) * w_all_capacity_1);
        w_vals_1 = (float*)realloc(w_vals_1, sizeof(float) * w_all_capacity_1);
      }
    }
 
    if (all_array_id == 0) {
      w_all_size = Merge_coord_t(w1_crd_0, w2_crd_0, w3_crd_0, w_vals_0, w1_crd_1, w2_crd_1, w3_crd_1, w_vals_1, w_all_size, w_accumulator_size);
    } else {
      w_all_size = Merge_coord_t(w1_crd_1, w2_crd_1, w3_crd_1, w_vals_1, w1_crd_0, w2_crd_0, w3_crd_0, w_vals_0, w_all_size, w_accumulator_size);
    }   
    all_array_id ^= 1;
    // Clear
    w_accumulator_size = 0;
  }

  if (all_array_id == 0) {
    w1_crd = w1_crd_0;
    w2_crd = w2_crd_0;
    w3_crd = w3_crd_0;
    w_vals = w_vals_0;
    free(w_vals_1);
    free(w1_crd_1);
    free(w2_crd_1);
    free(w3_crd_1);
  } else {
    w1_crd = w1_crd_1;
    w2_crd = w2_crd_1;
    w3_crd = w3_crd_1;
    w_vals = w_vals_1;
    free(w_vals_0);
    free(w1_crd_0);
    free(w2_crd_0);
    free(w3_crd_0);
  }

  A->indices[1][0] = (int32_t*)(w1_crd);
  A->indices[2][0] = (int32_t*)(w2_crd);
  A->indices[3][0] = (int32_t*)(w3_crd);
  A->vals = (float*)w_vals;
  A->vals_size = w_all_size;

  return 0;
}

void COO_CSF_DCSR_coord(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t w_cap, bool print = false) {
  acc_capacity = w_cap;
  w_accumulator.crd[0] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  w_accumulator.crd[1] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  w_accumulator.crd[2] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  w_accumulator.val = (float*)malloc(sizeof(float) * acc_capacity);
  w_accumulator_index = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  compute_coord(A,B,C,w_cap);
  free(w_accumulator.crd[0]);
  free(w_accumulator.crd[1]);
  free(w_accumulator.crd[2]);
  free(w_accumulator.val);
  free(w_accumulator_index);
  if (print) {
    print_taco_tensor_DC(C);
    print_taco_tensor_COO(A);
  }
  // free(A->vals);
  // free(A->indices[1][0]);
  // free(A->indices[2][0]);
  // free(A->indices[3][0]);
}