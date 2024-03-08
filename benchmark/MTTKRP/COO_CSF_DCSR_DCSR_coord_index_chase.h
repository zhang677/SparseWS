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

int32_t Merge_coord_t(int32_t* COO1_crd, int32_t* COO2_crd, float* COO_vals, int32_t COO_size, int32_t accumulator_size) {
  if (COO_size == 0) {
    COO1_crd[0] = w_accumulator.crd[0][w_accumulator_index[0]];
    COO2_crd[0] = w_accumulator.crd[1][w_accumulator_index[0]];
    COO_vals[0] = w_accumulator.val[w_accumulator_index[0]];
    int target_pointer = 0;
    for (int i = 1; i < accumulator_size; i++) {
      if (compare(w_accumulator_index[i], w_accumulator_index[target_pointer], COO1_crd, COO2_crd) == 0) {
        COO_vals[target_pointer] += w_accumulator.val[w_accumulator_index[i]];
      } else {
        target_pointer++;
        COO1_crd[target_pointer] = w_accumulator.crd[0][w_accumulator_index[i]];
        COO2_crd[target_pointer] = w_accumulator.crd[1][w_accumulator_index[i]];
        COO_vals[target_pointer] = w_accumulator.val[w_accumulator_index[i]];
      }
    }
    return target_pointer + 1;
  }
  int32_t* tmp_COO_crd[2];
  float* tmp_COO_vals;
  tmp_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_vals = (float*)malloc(sizeof(float) * (accumulator_size + COO_size));
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
  for (int i = 0; i <= target_pointer; i++) {
    COO1_crd[i] = tmp_COO_crd[0][i];
    COO2_crd[i] = tmp_COO_crd[1][i];
    COO_vals[i] = tmp_COO_vals[i];
  }
  free(tmp_COO_crd[0]);
  free(tmp_COO_crd[1]);
  free(tmp_COO_vals);
  return target_pointer + 1;
}

int compute(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C, taco_tensor_t *D, int32_t w_accumulator_capacity) {
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
                // std::cout << "TryInsert: " << w_point[0] << " , " << w_point[1] << std::endl;
                // TryInsert_hash(w_insertFail, &w_accumulator, w_point, (B_vals[iB] * C_vals[jC] * D_vals[jD]), A2_dimension);
                w_accumulator_size = TryInsert_coord_t(&w_insertFail, w_accumulator_size, w_point, (B_vals[iB] * C_vals[jC] * D_vals[jD]));

                if (w_insertFail) {
                    if (w_accumulator_size + w_all_size > w_all_capacity) {
                        w_all_capacity = w_accumulator_size + w_all_size;
                        w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
                        w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
                        w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
                    }
                    // w_all_size = Merge_hash(w1_crd, w2_crd, w_vals, w_all_size, &w_accumulator);
                    // refresh_wspace(&w_accumulator);
                    // TryInsert_hash(w_insertFail, &w_accumulator, w_point, (B_vals[iB] * C_vals[jC] * D_vals[jD]), A2_dimension);
                    w_all_size = Merge_coord_t(w1_crd, w2_crd, w_vals, w_all_size, w_accumulator_size);
                    w_accumulator_index[0] = 0;
                    w_accumulator_size = 0;
                    w_accumulator_size = TryInsert_coord_t(&w_insertFail, w_accumulator_size, w_point, (B_vals[iB] * C_vals[jC] * D_vals[jD]));
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
  if (w_accumulator_size > 0) {
    // Sort
    w_accumulator_size = Sort_t(w_accumulator_size, false);
    // Enlarge
    if (w_accumulator_size + w_all_size > w_all_capacity) {
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

  A->indices[1][0] = (int32_t*)(w1_crd); 
  A->indices[2][0] = (int32_t*)(w2_crd); // It should be A->indices[2][0]. However, it works because we initial the CSR and CSR has indices[1][1]. Therefore, don't change the init code!!!
  A->vals = (float*)(w_vals);
  A->vals_size = w_all_size;

  return 0;
}

void COO_CSF_DCSR_DCSR_coord(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, taco_tensor_t* D, int32_t w_cap, bool print = false) {
  acc_capacity = w_cap;
  w_accumulator.crd[0] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  w_accumulator.crd[1] = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  w_accumulator.val = (float*)malloc(sizeof(float) * acc_capacity);
  w_accumulator_index = (int32_t*)malloc(sizeof(int32_t) * acc_capacity);
  compute(A,B,C,D,w_cap);
  free(w_accumulator.crd[0]);
  free(w_accumulator.crd[1]);
  free(w_accumulator.val);
  free(w_accumulator_index);
  if (print) {
    print_taco_tensor_DC(C);
    print_taco_tensor_DC(D);
    print_taco_tensor_COO(A);
  }
  free(A->vals);
  free(A->indices[1][0]);
  free(A->indices[2][0]);
}