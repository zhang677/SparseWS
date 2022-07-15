#ifndef DENSE_DOT_H
#define DENSE_DOT_H
#include "../utils/dataloader.h"
#include "../utils/lib.h"
#include <iostream>

int benchmark(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C) {
  float* restrict A_vals = (float*)(A->vals);
  int B1_dimension = (int)(B->dimensions[0]);
  float* restrict B_vals = (float*)(B->vals);
  int C1_dimension = (int)(C->dimensions[0]);
  float* restrict C_vals = (float*)(C->vals);

  float A_val = 0.0;
  A_vals = (float*)malloc(sizeof(float) * 1);
  for(int i = 0 ; i<C1_dimension ; i++) {
      A_val += B_vals[i] * C_vals[i];
  }

  A_vals[0] = A_val;
  A->vals = (float*)A_vals;
  return 0;
}


int compute(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C) {
  float* restrict A_vals = (float*)(A->vals);
  int B1_dimension = (int)(B->dimensions[0]);
  float* restrict B_vals = (float*)(B->vals);
  int C1_dimension = (int)(C->dimensions[0]);
  float* restrict C_vals = (float*)(C->vals);

  float A_val = 0.0;
  A_vals = (float*)malloc(sizeof(float) * 1);

  int32_t ws_index_list_size = 0;
  float* restrict ws = 0;
  int32_t* restrict ws_index_list = 0;
  ws_index_list = (int32_t*)malloc(sizeof(int32_t) * C1_dimension);
  bool* restrict ws_already_set = calloc(C1_dimension, sizeof(bool));
  ws = (float*)malloc(sizeof(float) * C1_dimension);
  float* restrict B_new = 0;
  B_new = (float*)malloc(sizeof(float) * C1_dimension);
  for (int32_t pB_new = 0; pB_new < C1_dimension; pB_new++) {
    B_new[pB_new] = 0.0;
  }
  for (int32_t i0 = 0; i0 < 32; i0++) {
    for (int32_t i1 = 0; i1 < 32; i1++) {
      int32_t i_bounded = i0 * 32 + i1;
      B_new[i_bounded] = B_vals[i_bounded];
    }
  }
  float* restrict C_new = 0;
  C_new = (float*)malloc(sizeof(float) * C1_dimension);
  for (int32_t pC_new = 0; pC_new < C1_dimension; pC_new++) {
    C_new[pC_new] = 0.0;
  }
  for (int32_t i0 = 0; i0 < 32; i0++) {
    for (int32_t i1 = 0; i1 < 32; i1++) {
      int32_t i_bounded = i0 * 32 + i1;
      C_new[i_bounded] = C_vals[i_bounded];
    }
  }
  for (int32_t i0 = 0; i0 < 32; i0++) {
    for (int32_t i1 = 0; i1 < 32; i1++) {
      int32_t i_bounded = i0 * 32 + i1;
      if (!ws_already_set[i_bounded]) {
        ws[i_bounded] = B_new[i_bounded] * C_new[i_bounded];
        ws_index_list[ws_index_list_size] = i_bounded;
        ws_already_set[i_bounded] = 1;
        ws_index_list_size++;
      }
      else {
        ws[i_bounded] = B_new[i_bounded] * C_new[i_bounded];
      }
    }
  }
  free(C_new);
  free(B_new);
  for (int32_t i0 = 0; i0 < 32; i0++) {
    for (int32_t ws_index_locator = 0; ws_index_locator < ws_index_list_size; ws_index_locator++) {
      int32_t i1 = ws_index_list[ws_index_locator];
      int32_t i_bounded = i0 * 32 + i1;
      A_val += ws[i_bounded];
      ws_already_set[i1] = 0;
    }
  }
  free(ws_index_list);
  free(ws_already_set);
  free(ws);

  A_vals[0] = A_val;

  A->vals = (float*)A_vals;
  return 0;
}

void test_dense_dot() {
    int N = 1024;
    EigenVector B_eigen = gen_one_vector(N);
    taco_tensor_t B = to_taco_tensor(B_eigen);
    EigenVector C_eigen = gen_one_vector(N);
    taco_tensor_t C = to_taco_tensor(C_eigen);
    taco_tensor_t A;
    compute(&A,&B,&C);
    std::cout<<A.vals[0]<<std::endl;
    taco_tensor_t A_ref;
    benchmark(&A_ref,&B,&C);
    std::cout<<A_ref.vals[0]<<std::endl;
}
#endif