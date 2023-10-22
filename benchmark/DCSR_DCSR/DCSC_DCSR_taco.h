#include "../../utils/dataloader.h"
#include "../../utils/lib.h"

int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
  int C1_dimension = (int)(C->dimensions[0]);
  int A1_dimension = (int)(A->dimensions[0]);
  int B2_dimension = (int)(B->dimensions[1]);
  int* restrict C2_pos = (int*)(C->indices[1][0]);
  int* restrict C2_crd = (int*)(C->indices[1][1]);
  float* restrict C_vals = (float*)(C->vals);
  int* restrict A1_pos = (int*)(A->indices[0][0]);
  int* restrict A1_crd = (int*)(A->indices[0][1]);
  int* restrict A2_pos = (int*)(A->indices[1][0]);
  int* restrict A2_crd = (int*)(A->indices[1][1]);
  float* restrict A_vals = (float*)(A->vals);
  int* restrict B1_pos = (int*)(B->indices[0][0]);
  int* restrict B1_crd = (int*)(B->indices[0][1]);
  int* restrict B2_pos = (int*)(B->indices[1][0]);
  int* restrict B2_crd = (int*)(B->indices[1][1]);
  float* restrict B_vals = (float*)(B->vals);

  int32_t* restrict C2_nnz = (int32_t*)calloc(A1_dimension, sizeof(int32_t));

  int32_t* restrict qw_index_list = 0;
  qw_index_list = (int32_t*)malloc(sizeof(int32_t) * A1_dimension * B2_dimension);
  bool* restrict qw_already_set = (bool*)calloc(A1_dimension * B2_dimension, sizeof(bool));
  int qw_index_list_size = 0;


  int32_t qjA = A1_pos[0];
  int32_t qpA1_end = A1_pos[1];
  int32_t qjB = B1_pos[0];
  int32_t qpB1_end = B1_pos[1];

  while (qjA < qpA1_end && qjB < qpB1_end) {
    int32_t qjA0 = A1_crd[qjA];
    int32_t qjB0 = B1_crd[qjB];
    int32_t qj = TACO_MIN(qjA0, qjB0);
    if (qjA0 == qj && qjB0 == qj) {
      for (int32_t qiA = A2_pos[qjA]; qiA < A2_pos[qjA + 1]; qiA++) {
        int32_t qi = A2_crd[qiA];
        for (int32_t qkB = B2_pos[qjB]; qkB < B2_pos[qjB + 1]; qkB++) {
          int32_t qk = B2_crd[qkB];
          if (!qw_already_set[qi * B2_dimension + qk]) {
            qw_already_set[qi * B2_dimension + qk] = true;
            qw_index_list[qw_index_list_size++] = qi * B2_dimension + qk;
          }
        }
      }
    }
    qjA += (int32_t)(qjA0 == qj);
    qjB += (int32_t)(qjB0 == qj);
  }

  for (int32_t qw_index_locator = 0; qw_index_locator < qw_index_list_size; qw_index_locator++) {
    int32_t qk = qw_index_list[qw_index_locator] % B2_dimension;
    int32_t qi = qw_index_list[qw_index_locator] / B2_dimension;
    C2_nnz[qi] = C2_nnz[qi] + 1;
    qw_already_set[qi * B2_dimension + qk] = false;
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
  w_index_list = (int32_t*)malloc(sizeof(int32_t) * A1_dimension * B2_dimension);
  bool* restrict w_already_set = (bool*)calloc(A1_dimension * B2_dimension, sizeof(bool));
  w = (float*)malloc(sizeof(float) * A1_dimension * B2_dimension);

  int32_t w_index_list_size = 0;

  int32_t jA = A1_pos[0];
  int32_t pA1_end = A1_pos[1];
  int32_t jB = B1_pos[0];
  int32_t pB1_end = B1_pos[1];

  while (jA < pA1_end && jB < pB1_end) {
    int32_t jA0 = A1_crd[jA];
    int32_t jB0 = B1_crd[jB];
    int32_t j = TACO_MIN(jA0, jB0);
    if (jA0 == j && jB0 == j) {
      for (int32_t iA = A2_pos[jA]; iA < A2_pos[jA + 1]; iA++) {
        int32_t i = A2_crd[iA];
        for (int32_t kB = B2_pos[jB]; kB < B2_pos[jB + 1]; kB++) {
          int32_t k = B2_crd[kB];
          if (!w_already_set[i * B2_dimension + k]) {
            w_already_set[i * B2_dimension + k] = true;
            w_index_list[w_index_list_size++] = i * B2_dimension + k;
            w[i * B2_dimension + k] = A_vals[iA] * B_vals[kB];
          }
          else {
            w[i * B2_dimension + k] = w[i * B2_dimension + k] + A_vals[iA] * B_vals[kB];
          }
        }
      }
    }
    jA += (int32_t)(jA0 == j);
    jB += (int32_t)(jB0 == j);
  }
  qsort(w_index_list, w_index_list_size, sizeof(int32_t), cmp);
  for (int32_t w_index_locator = 0; w_index_locator < w_index_list_size; w_index_locator++) {
    int32_t k = w_index_list[w_index_locator] % B2_dimension;
    int32_t i = w_index_list[w_index_locator] / B2_dimension;
    int32_t pC2 = C2_pos[i];
    C2_pos[i] = C2_pos[i] + 1;
    C2_crd[pC2] = k;
    C_vals[pC2] = w[i * B2_dimension + k];
    w_already_set[i * B2_dimension + k] = false;
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

double DCSC_DCSR_taco(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t warmup, int32_t repeat, bool bench = false, bool print = false) {
  // std::cout << "Capacity: " << w_cap << std::endl;
  for (int i = 0; i < warmup; i++) {
    compute(C,A,B);
    if (bench) {
      free(C->vals);
      free(C->indices[1][0]);
      free(C->indices[1][1]);
    }
  }
  double start = clock();
  for (int i = 0; i < repeat; i++) {
    compute(C,A,B);
    if (bench && i != repeat - 1) {
      free(C->vals);
      free(C->indices[1][0]);
      free(C->indices[1][1]);
    }
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