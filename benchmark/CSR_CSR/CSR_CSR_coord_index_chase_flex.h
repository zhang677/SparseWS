#include "../../utils/dataloader.h"
#include "../../utils/lib.h"

#define w_order 2

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

int compare(int32_t l, int32_t r, int32_t* crd0, int32_t* crd1) {
  if (w_accumulator.crd[0][l] == crd0[r]) {
    if (w_accumulator.crd[1][l] == crd1[r]) {
      return 0;
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
    w_accumulator.val[accumulator_size] = val;
    *insertFail = false;
    accumulator_size ++;
    return accumulator_size;
  }
}

int32_t Merge_coord_t(int32_t* COO1_crd_0, int32_t* COO2_crd_0, float* COO_vals_0,int32_t* COO1_crd_1, int32_t* COO2_crd_1, float* COO_vals_1, int32_t COO_size, int32_t accumulator_size) {
  int32_t* COO1_crd;
  int32_t* COO2_crd;
  float* COO_vals;
  COO1_crd = COO1_crd_0;
  COO2_crd = COO2_crd_0;
  COO_vals = COO_vals_0;
  int32_t* tmp_COO_crd[2];
  float* tmp_COO_vals;
  tmp_COO_crd[0] = COO1_crd_1;
  tmp_COO_crd[1] = COO2_crd_1;
  tmp_COO_vals = COO_vals_1;

  if (COO_size == 0) {
    tmp_COO_crd[0][0] = w_accumulator.crd[0][w_accumulator_index[0]];
    tmp_COO_crd[1][0] = w_accumulator.crd[1][w_accumulator_index[0]];
    tmp_COO_vals[0] = w_accumulator.val[w_accumulator_index[0]];
    int target_pointer = 0;
    for (int i = 1; i < accumulator_size; i++) {
      if (compare(w_accumulator_index[i], w_accumulator_index[target_pointer], COO1_crd, COO2_crd) == 0) {
        tmp_COO_vals[target_pointer] += w_accumulator.val[w_accumulator_index[i]];
      } else {
        target_pointer++;
        tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[i]];
        tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[i]];
        tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[i]];
      }
    }
    return target_pointer + 1;
  }
  int accumulator_pointer = 0;
  int content_pointer = 0;
  int target_pointer = 0;
  int cmp = compare(w_accumulator_index[accumulator_pointer], content_pointer, COO1_crd, COO2_crd);
  if (cmp == 0) {
    tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
    tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
    tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]] + COO_vals[content_pointer];
    accumulator_pointer ++;
    content_pointer ++;
  } else if (cmp < 0) {
    tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
    tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
    tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]];
    accumulator_pointer ++;
  } else {
    tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
    tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
    tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
    content_pointer ++;
  }

  while(accumulator_pointer < accumulator_size && content_pointer < COO_size) {
    cmp = compare(w_accumulator_index[accumulator_pointer], content_pointer, COO1_crd, COO2_crd);
    if (cmp == 0) {
      if (compare(w_accumulator_index[accumulator_pointer], target_pointer, tmp_COO_crd[0], tmp_COO_crd[1]) == 0) {
        tmp_COO_vals[target_pointer] += w_accumulator.val[w_accumulator_index[accumulator_pointer]] + COO_vals[content_pointer];
      } else {
        target_pointer ++;
        tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
        tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
        tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]] + COO_vals[content_pointer];
      }
      accumulator_pointer ++;
      content_pointer ++;
    } else if (cmp < 0) {
      if (compare(w_accumulator_index[accumulator_pointer], target_pointer, tmp_COO_crd[0], tmp_COO_crd[1]) == 0) {
        tmp_COO_vals[target_pointer] += w_accumulator.val[w_accumulator_index[accumulator_pointer]];
      } else {
        target_pointer ++;
        tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
        tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
        tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]];
      }
      accumulator_pointer ++;
    } else {
      if(tmp_COO_crd[0][target_pointer] == COO1_crd[content_pointer] && tmp_COO_crd[1][target_pointer] == COO2_crd[content_pointer]) {
        tmp_COO_vals[target_pointer] += COO_vals[content_pointer];
      } else {
        target_pointer ++;
        tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
        tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
        tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      }
      content_pointer ++;
    }
  }
  while(accumulator_pointer<accumulator_size) {
    if (compare(w_accumulator_index[accumulator_pointer], target_pointer, tmp_COO_crd[0], tmp_COO_crd[1]) == 0) {
      tmp_COO_vals[target_pointer] += w_accumulator.val[w_accumulator_index[accumulator_pointer]];
    } else {
      target_pointer ++;
      tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
      tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
      tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]];
    }
    accumulator_pointer ++;
  }
  while(content_pointer<COO_size) {
    if(tmp_COO_crd[0][target_pointer] == COO1_crd[content_pointer] && tmp_COO_crd[1][target_pointer] == COO2_crd[content_pointer]) {
      tmp_COO_vals[target_pointer] += COO_vals[content_pointer];
    } else {
      target_pointer ++;
      tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
      tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
    }
    content_pointer ++;
  }
  return target_pointer + 1;
}

int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
  int C1_dimension = (int)(C->dimensions[0]);
  int* restrict C2_pos = (int*)(C->indices[1][0]);
  int* restrict C2_crd = (int*)(C->indices[1][1]);
  float* restrict C_vals = (float*)(C->vals);
  int A1_dimension = (int)(A->dimensions[0]);
  int* restrict A2_pos = (int*)(A->indices[1][0]);
  int* restrict A2_crd = (int*)(A->indices[1][1]);
  float* restrict A_vals = (float*)(A->vals);
  int B2_dimension = (int)(B->dimensions[1]);
  int* restrict B2_pos = (int*)(B->indices[1][0]);
  int* restrict B2_crd = (int*)(B->indices[1][1]);
  float* restrict B_vals = (float*)(B->vals);

  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    C2_pos[pC2] = 0;
  }

  int32_t w_accumulator_size = 0;
  int32_t w_all_capacity_0 = acc_capacity; 
  int32_t w_all_capacity_1 = acc_capacity; 
  int32_t w_all_size = 0;
  int32_t* w1_crd_0 = 0;
  int32_t* w2_crd_0 = 0;
  float* w_vals_0 = 0;
  w1_crd_0 = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity_0);
  w2_crd_0 = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity_0);
  w_vals_0 = (float*)malloc(sizeof(float) * w_all_capacity_0);
  int32_t* w1_crd_1 = 0;
  int32_t* w2_crd_1 = 0;
  float* w_vals_1 = 0;
  w1_crd_1 = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity_1);
  w2_crd_1 = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity_1);
  w_vals_1 = (float*)malloc(sizeof(float) * w_all_capacity_1);
  int32_t* w1_crd = 0;
  int32_t* w2_crd = 0;
  float* w_vals = 0;
  int32_t all_array_id = 0;

  bool w_insertFail = false; 
  int32_t w_point[w_order];
  int32_t w1_pos[w_order];

  //int counter = 0;

  for (int32_t i = 0; i < A1_dimension; i++) {
    for (int32_t jA = A2_pos[i]; jA < A2_pos[i+1]; jA++) {
      int32_t j = A2_crd[jA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[j+1]; kB++) {
        int32_t k = B2_crd[kB];
        w_point[0] = i;
        w_point[1] = k;

        w_accumulator_size = TryInsert_coord_t(&w_insertFail, w_accumulator_size, w_point, A_vals[jA] * B_vals[kB]);
        if (w_insertFail) { 
          //counter += 1; 
          if ((all_array_id == 0 && w_accumulator_size + w_all_size > w_all_capacity_1) ||
              (all_array_id == 1 && w_accumulator_size + w_all_size > w_all_capacity_0)) {
            //w_all_capacity = w_all_capacity * 2;
            
            if (all_array_id == 1) {
              w_all_capacity_0 = w_accumulator_size + w_all_size;
              w1_crd_0 = (int32_t*)realloc(w1_crd_0, sizeof(int32_t) * w_all_capacity_0);
              w2_crd_0 = (int32_t*)realloc(w2_crd_0, sizeof(int32_t) * w_all_capacity_0);
              w_vals_0 = (float*)realloc(w_vals_0, sizeof(float) * w_all_capacity_0);
            } else {
              w_all_capacity_1 = w_accumulator_size + w_all_size;
              w1_crd_1 = (int32_t*)realloc(w1_crd_1, sizeof(int32_t) * w_all_capacity_1);
              w2_crd_1 = (int32_t*)realloc(w2_crd_1, sizeof(int32_t) * w_all_capacity_1);
              w_vals_1 = (float*)realloc(w_vals_1, sizeof(float) * w_all_capacity_1);
            }
          }
          if (all_array_id == 0) {
            w_all_size = Merge_coord_t(w1_crd_0, w2_crd_0, w_vals_0, w1_crd_1, w2_crd_1, w_vals_1, w_all_size, w_accumulator_size);
          } else {
            w_all_size = Merge_coord_t(w1_crd_1, w2_crd_1, w_vals_1, w1_crd_0, w2_crd_0, w_vals_0, w_all_size, w_accumulator_size);
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
          free(w_accumulator.val);
          free(w_accumulator_index);
          w_accumulator.crd[0] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
          w_accumulator.crd[1] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
          w_accumulator.val = (float*)malloc(sizeof(float) * acc_capacity);
          w_accumulator_index = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
          w_accumulator_size = 0;
          w_accumulator_index[0] = 0;
          w_accumulator_size = TryInsert_coord_t(&w_insertFail, w_accumulator_size, w_point, A_vals[jA] * B_vals[kB]);
        }
      }
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
        w_vals_0 = (float*)realloc(w_vals_0, sizeof(float) * w_all_capacity_0);
      } else {
        w_all_capacity_1 = w_accumulator_size + w_all_size;
        w1_crd_1 = (int32_t*)realloc(w1_crd_1, sizeof(int32_t) * w_all_capacity_1);
        w2_crd_1 = (int32_t*)realloc(w2_crd_1, sizeof(int32_t) * w_all_capacity_1);
        w_vals_1 = (float*)realloc(w_vals_1, sizeof(float) * w_all_capacity_1);
      }
    }
 
    if (all_array_id == 0) {
      w_all_size = Merge_coord_t(w1_crd_0, w2_crd_0, w_vals_0, w1_crd_1, w2_crd_1, w_vals_1, w_all_size, w_accumulator_size);
    } else {
      w_all_size = Merge_coord_t(w1_crd_1, w2_crd_1, w_vals_1, w1_crd_0, w2_crd_0, w_vals_0, w_all_size, w_accumulator_size);
    }   
    all_array_id ^= 1;
    // Clear
    w_accumulator_size = 0;
  }
  if (all_array_id == 0) {
    w1_crd = w1_crd_0;
    w2_crd = w2_crd_0;
    w_vals = w_vals_0;
    free(w_vals_1);
    free(w1_crd_1);
    free(w2_crd_1);
  } else {
    w1_crd = w1_crd_1;
    w2_crd = w2_crd_1;
    w_vals = w_vals_1;
    free(w_vals_0);
    free(w1_crd_0);
    free(w2_crd_0);
  }
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

  free(w1_crd);

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

void CSR_CSR_coord(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t w_cap, bool print = false) {
  // C(i,k) = A(i,j) * B(j,k); C: CSR, A: CSR, B: CSR
  // std::cout << "Capacity: " << w_cap << std::endl;
  acc_capacity = w_cap;
  w_accumulator.crd[0] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  w_accumulator.crd[1] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  w_accumulator.val = (float*)malloc(sizeof(float) * acc_capacity);
  w_accumulator_index = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  compute(C,A,B);
  free(w_accumulator.crd[0]);
  free(w_accumulator.crd[1]);
  free(w_accumulator.val);
  free(w_accumulator_index);
  if (print) {
    print_taco_tensor_DC(A);
    print_taco_tensor_DC(B);
    print_taco_tensor_DC(C);
  }
}