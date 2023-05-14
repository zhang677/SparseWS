#include "../utils/dataloader.h"
#include "../utils/lib.h"

#define w_order 2
const int32_t INTMAX = 2147483647;

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
  return w_accumulator.crd[1][l] - w_accumulator.crd[1][r];
}

int esc_cmp_t_rev(const void* a, const void* b) {
  int l = *(int32_t *)a;
  int r = *(int32_t *)b;
  for (int i = 0; i < w_order; i++) {
    if (w_accumulator.crd[i][r] == w_accumulator.crd[i][l]) continue;
    return w_accumulator.crd[i][r] - w_accumulator.crd[i][l];
  }
  return w_accumulator.crd[1][r] - w_accumulator.crd[1][l];
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

int32_t Merge_coord_t(int32_t* COO1_crd, int32_t* COO2_crd, float* COO_vals, int32_t COO_size, int32_t accumulator_size) {
  int p1 = 0;
  int p2 = 1;
  while (p2 < accumulator_size) {
    if (esc_cmp_t(&w_accumulator_index[p2], &w_accumulator_index[p1]) == 0) {
      w_accumulator.val[w_accumulator_index[p1]] += w_accumulator.val[w_accumulator_index[p2]];
    } else {
      if (p2 - p1 > 1) {
        p1++;
        w_accumulator.crd[0][w_accumulator_index[p1]] = w_accumulator.crd[0][w_accumulator_index[p2]];
        w_accumulator.crd[1][w_accumulator_index[p1]] = w_accumulator.crd[1][w_accumulator_index[p2]];
        w_accumulator.val[w_accumulator_index[p1]] = w_accumulator.val[w_accumulator_index[p2]];
      } else {
        p1++;
      }
    }
    p2++;
  }
  accumulator_size = p1 + 1;
  if (COO_size == 0) {
    for (int i = 0; i < accumulator_size; i++) {
      COO1_crd[i] = w_accumulator.crd[0][w_accumulator_index[i]];
      COO2_crd[i] = w_accumulator.crd[1][w_accumulator_index[i]];
      COO_vals[i] = w_accumulator.val[w_accumulator_index[i]];
    }
    return accumulator_size;
  }
  int32_t* tmp_COO_crd[2];
  float* tmp_COO_vals;
  tmp_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_vals = (float*)malloc(sizeof(float) * (accumulator_size + COO_size));
  int accumulator_pointer = 0;
  int content_pointer = 0;
  int target_pointer = 0;
  while(accumulator_pointer < accumulator_size && content_pointer < COO_size) {
    int cmp = compare(w_accumulator_index[accumulator_pointer], content_pointer, COO1_crd, COO2_crd);
    if (cmp == 0) {
      tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
      tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
      tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]] + COO_vals[content_pointer];
      accumulator_pointer ++;
      content_pointer ++;
      target_pointer ++;
    } else if (cmp < 0) {
      tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
      tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
      tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]];
      accumulator_pointer ++;
      target_pointer ++;
    } else {
      tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
      tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      content_pointer ++;
      target_pointer ++;
    }
  }
  while(accumulator_pointer<accumulator_size) {
    tmp_COO_crd[0][target_pointer] = w_accumulator.crd[0][w_accumulator_index[accumulator_pointer]];
    tmp_COO_crd[1][target_pointer] = w_accumulator.crd[1][w_accumulator_index[accumulator_pointer]];
    tmp_COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[accumulator_pointer]];
    accumulator_pointer ++;
    target_pointer ++;
  }
  while(content_pointer<COO_size) {
    tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
    tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
    tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
    content_pointer ++;
    target_pointer ++;
  }
  for (int i = 0; i < target_pointer; i++) {
    COO1_crd[i] = tmp_COO_crd[0][i];
    COO2_crd[i] = tmp_COO_crd[1][i];
    COO_vals[i] = tmp_COO_vals[i];
  }
  free(tmp_COO_crd[0]);
  free(tmp_COO_crd[1]);
  free(tmp_COO_vals);
  return target_pointer;
}

int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B, int32_t w_cap) {
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

  int32_t w_accumulator_size = 0;
  int32_t w_all_capacity = acc_capacity; 
  int32_t w_all_size = 0;
  int32_t* w1_crd = 0;
  int32_t* w2_crd = 0;
  float* w_vals = 0;
  w1_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  w2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  w_vals = (float*)malloc(sizeof(float) * w_all_capacity);
  bool w_insertFail = false; 
  int32_t w_point[w_order];
  int32_t w1_pos[w_order];
  for (int32_t j = 0; j < B1_dimension; j++) {
    for (int32_t iA = A2_pos[j]; iA < A2_pos[j+1]; iA++) {
      int32_t i = A2_crd[iA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[j+1]; kB++) {
        int32_t k = B2_crd[kB];
        w_point[0] = k;
        w_point[1] = i; // Transpose

        w_accumulator_size = TryInsert_coord_t(&w_insertFail, w_accumulator_size, w_point, A_vals[iA] * B_vals[kB]);
        if (w_insertFail) {
          if(w_accumulator_size + w_all_size > w_all_capacity) {
            w_all_capacity = w_accumulator_size + w_all_size;
            //std::cout << w_all_capacity << std::endl;
            w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
            w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
            w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
          }
          w_all_size = Merge_coord_t(w1_crd, w2_crd, w_vals, w_all_size, w_accumulator_size);
          w_accumulator_index[0] = 0;
          w_accumulator_size = 0;
          w_accumulator_size = TryInsert_coord_t(&w_insertFail, w_accumulator_size, w_point, A_vals[iA] * B_vals[kB]);
        }
      }
    }
  }
  if (w_accumulator_size > 0) {
    // Sort
    w_accumulator_size = Sort_t(w_accumulator_size, false);
    // Enlarge
    if(w_accumulator_size + w_all_size > w_all_capacity) {
        w_all_capacity = w_accumulator_size + w_all_size;
        w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
        w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
        w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
    }
    // Merge
    w_all_size = Merge_coord_t(w1_crd, w2_crd, w_vals, w_all_size, w_accumulator_size);
    // Clear
    w_accumulator_size = 0;
    w_accumulator_index[0] = 0;
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

void CSC_CSR_T_coord(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t w_cap, bool print = false) {
  // C(k,i) = A(i,j) * B(j,k); C: CSR, A: CSR, B: CSR
  // std::cout << "Capacity: " << w_cap << std::endl;
  acc_capacity = w_cap;
  w_accumulator.crd[0] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  w_accumulator.crd[1] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  w_accumulator.val = (float*)malloc(sizeof(float) * acc_capacity);
  w_accumulator_index = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  compute(C,A,B,w_cap);
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