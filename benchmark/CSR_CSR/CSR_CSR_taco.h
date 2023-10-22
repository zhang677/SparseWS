#include "../../utils/dataloader.h"
#include "../../utils/lib.h"

int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
  int C1_dimension = (int)(C->dimensions[0]);
  int* restrict C2_pos = (int*)(C->indices[1][0]);
  int* restrict C2_crd = (int*)(C->indices[1][1]);
  float* restrict C_vals = (float*)(C->vals);
  int A1_dimension = (int)(A->dimensions[0]);
  int* restrict A2_pos = (int*)(A->indices[1][0]);
  int* restrict A2_crd = (int*)(A->indices[1][1]);
  float* restrict A_vals = (float*)(A->vals);
  int B1_dimension = (int)(B->dimensions[0]);
  int B2_dimension = (int)(B->dimensions[1]);
  int* restrict B2_pos = (int*)(B->indices[1][0]);
  int* restrict B2_crd = (int*)(B->indices[1][1]);
  float* restrict B_vals = (float*)(B->vals);

  int32_t* restrict C2_nnz = (int32_t*)calloc(A1_dimension, sizeof(int32_t));

  int32_t* restrict qw_index_list = 0;
  qw_index_list = (int32_t*)malloc(sizeof(int32_t) * B2_dimension);
  bool* restrict qw_already_set = (bool*)calloc(B2_dimension, sizeof(bool));

  for (int32_t qi = 0; qi < A1_dimension; qi++) {
    int32_t qw_index_list_size = 0;
    for (int32_t qjA = A2_pos[qi]; qjA < A2_pos[(qi + 1)]; qjA++) {
      int32_t qj = A2_crd[qjA];
      for (int32_t qkB = B2_pos[qj]; qkB < B2_pos[(qj + 1)]; qkB++) {
        int32_t qk = B2_crd[qkB];
        if (!qw_already_set[qk]) {
          qw_index_list[qw_index_list_size] = qk;
          qw_already_set[qk] = 1;
          qw_index_list_size++;
        }
      }
    }
    for (int32_t qw_index_locator = 0; qw_index_locator < qw_index_list_size; qw_index_locator++) {
      int32_t qk = qw_index_list[qw_index_locator];
      C2_nnz[qi] = C2_nnz[qi] + 1;
      qw_already_set[qk] = 0;
    }
  }

  free(qw_index_list);
  free(qw_already_set);

  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t i = 0; i < C1_dimension; i++) {
    C2_pos[i + 1] = C2_pos[i] + C2_nnz[i];
  }
  C2_crd = (int32_t*)malloc(sizeof(int32_t) * C2_pos[C1_dimension]);
  C_vals = (float*)malloc(sizeof(float) * C2_pos[C1_dimension]);

  float* restrict w = 0;
  int32_t* restrict w_index_list = 0;
  w_index_list = (int32_t*)malloc(sizeof(int32_t) * B2_dimension);
  bool* restrict w_already_set = (bool*)calloc(B2_dimension, sizeof(bool));
  w = (float*)malloc(sizeof(float) * B2_dimension);

  for (int32_t i = 0; i < A1_dimension; i++) {
    int32_t w_index_list_size = 0;
    for (int32_t jA = A2_pos[i]; jA < A2_pos[(i + 1)]; jA++) {
      int32_t j = A2_crd[jA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[(j + 1)]; kB++) {
        int32_t k = B2_crd[kB];
        if (!w_already_set[k]) {
          w[k] = A_vals[jA] * B_vals[kB];
          w_index_list[w_index_list_size] = k;
          w_already_set[k] = 1;
          w_index_list_size++;
        }
        else {
          w[k] = w[k] + A_vals[jA] * B_vals[kB];
        }
      }
    }
    qsort(w_index_list, w_index_list_size, sizeof(int32_t), cmp);

    for (int32_t w_index_locator = 0; w_index_locator < w_index_list_size; w_index_locator++) {
      int32_t k = w_index_list[w_index_locator];
      int32_t pC2 = C2_pos[i];
      C2_pos[i] = C2_pos[i] + 1;
      C2_crd[pC2] = k;
      C_vals[pC2] = w[k];
      w_already_set[k] = 0;
    }
  }

  free(w_index_list);
  free(w_already_set);
  free(w);

  for (int32_t p = 0; p < C1_dimension; p++) {
    C2_pos[C1_dimension - p] = C2_pos[((C1_dimension - p) - 1)];
  }
  C2_pos[0] = 0;

  free(C2_nnz);

  C->indices[1][0] = (int32_t*)(C2_pos);
  C->indices[1][1] = (int32_t*)(C2_crd);
  C->vals = (float*)C_vals;
  return 0;
}

double CSR_CSR_taco(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t warmup, int32_t repeat, bool bench = false, bool print = false) {
  // C(i,k) = A(i,j) * B(j,k); C: CSR, A: CSR, B: CSR
  for (int i = 0; i < warmup; i++) {
    compute(C,A,B);
    if (bench) {
      free(C->vals);
      free(C->indices[1][0]);
      free(C->indices[1][1]);
    }
  }
  Timer timer;
  timer.reset();
  for (int i = 0; i < repeat; i++) {
    compute(C,A,B);
    if (bench && i != repeat - 1) {
      free(C->vals);
      free(C->indices[1][0]);
      free(C->indices[1][1]);
    }
  }
  double duration = timer.elapsed() / repeat;
  if (print) {
    print_taco_tensor_DC(A);
    print_taco_tensor_DC(B);
    print_taco_tensor_DC(C);
  }
  return duration;
}