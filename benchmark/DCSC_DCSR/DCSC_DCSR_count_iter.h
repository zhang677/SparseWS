#include "../../utils/dataloader.h"
#include <math.h> 

int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
  int C1_dimension = (int)(C->dimensions[0]);
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

  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    C2_pos[pC2] = 0;
  }
  
  int32_t* restrict w_point = 0;
  w_point = (int32_t*)malloc(sizeof(int32_t) * 2);

  int32_t jA = A1_pos[0];
  int32_t pA1_end = A1_pos[1];
  int32_t jB = B1_pos[0];
  int32_t pB1_end = B1_pos[1];
  float v;
  int iterations = 0;

  while (jA < pA1_end && jB < pB1_end) {
    int32_t jA0 = A1_crd[jA];
    int32_t jB0 = B1_crd[jB];
    int32_t j = TACO_MIN(jA0, jB0);
    if (jA0 == j && jB0 == j) {
      for (int32_t iA = A2_pos[jA]; iA < A2_pos[jA + 1]; iA++) {
        int32_t i = A2_crd[iA];
        w_point[0] = i;
        for (int32_t kB = B2_pos[jB]; kB < B2_pos[jB + 1]; kB++) {
          int32_t k = B2_crd[kB];
          w_point[1] = k;
          iterations++;
        }
      }
    }
    jA += (int32_t)(jA0 == j);
    jB += (int32_t)(jB0 == j);
  }
  return iterations;
}

int DCSC_DCSR_iteration(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C) {
  return compute(C, A, B);
}