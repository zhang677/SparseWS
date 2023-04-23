#include "../utils/dataloader.h"
#include "../utils/lib.h"

/*
  IndexExpr precomputedExpr = A(i, j) * B(j, k);
  C(i, k) = precomputedExpr;
  TensorVar W("w", Type(Float32,{(size_t)NUM_I, (size_t)NUM_K}), {Dense, Dense});
  IndexStmt stmt = C.getAssignment().concretize();
  stmt = stmt.reorder({j,i,k}).precompute(precomputedExpr, {i,k}, {i,k}, W);
*/

int assemble(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
  int C1_dimension = (int)(C->dimensions[0]);
  int* restrict C2_pos = (int*)(C->indices[1][0]);
  int* restrict C2_crd = (int*)(C->indices[1][1]);
  float* restrict C_vals = (float*)(C->vals);
  int A1_dimension = (int)(A->dimensions[0]);
  int B2_dimension = (int)(B->dimensions[1]);

  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    C2_pos[pC2] = 0;
  }
  int32_t C2_crd_size = 1048576;
  C2_crd = (int32_t*)malloc(sizeof(int32_t) * C2_crd_size);
  int32_t kC = 0;

  for (int32_t i = 0; i < A1_dimension; i++) {
    int32_t pC2_begin = kC;

    for (int32_t k = 0; k < B2_dimension; k++) {
      if (C2_crd_size <= kC) {
        C2_crd = (int32_t*)realloc(C2_crd, sizeof(int32_t) * (C2_crd_size * 2));
        C2_crd_size *= 2;
      }
      C2_crd[kC] = k;
      kC++;
    }

    C2_pos[i + 1] = kC - pC2_begin;
  }

  int32_t csC2 = 0;
  for (int32_t pC20 = 1; pC20 < (C1_dimension + 1); pC20++) {
    csC2 += C2_pos[pC20];
    C2_pos[pC20] = csC2;
  }

  C_vals = (float*)malloc(sizeof(float) * kC);

  C->indices[1][0] = (uint8_t*)(C2_pos);
  C->indices[1][1] = (uint8_t*)(C2_crd);
  C->vals = (uint8_t*)C_vals;
  return 0;
}

int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
  int C1_dimension = (int)(C->dimensions[0]);
  float* restrict C_vals = (float*)(C->vals);
  int A1_dimension = (int)(A->dimensions[0]);
  int A2_dimension = (int)(A->dimensions[1]);
  int* restrict A2_pos = (int*)(A->indices[1][0]);
  int* restrict A2_crd = (int*)(A->indices[1][1]);
  float* restrict A_vals = (float*)(A->vals);
  int B1_dimension = (int)(B->dimensions[0]);
  int B2_dimension = (int)(B->dimensions[1]);
  int* restrict B2_pos = (int*)(B->indices[1][0]);
  int* restrict B2_crd = (int*)(B->indices[1][1]);
  float* restrict B_vals = (float*)(B->vals);

  int32_t kC = 0;

  float* restrict w = 0;
  w = (float*)malloc(sizeof(float) * (A1_dimension * B2_dimension));
  for (int32_t pw = 0; pw < (A1_dimension * B2_dimension); pw++) {
    w[pw] = 0.0;
  }
  for (int32_t j = 0; j < B1_dimension; j++) {
    for (int32_t iA = A2_pos[j]; iA < A2_pos[(j + 1)]; iA++) {
      int32_t i = A2_crd[iA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[(j + 1)]; kB++) {
        int32_t k = B2_crd[kB];
        int32_t kw = i * B2_dimension + k;
        w[kw] = w[kw] + A_vals[iA] * B_vals[kB];
      }
    }
  }
  for (int32_t i = 0; i < A1_dimension; i++) {
    for (int32_t k = 0; k < B2_dimension; k++) {
      int32_t kw = i * B2_dimension + k;
      C_vals[kC] = w[kw];
      kC++;
    }
  }
  free(w);

  C->vals = (uint8_t*)C_vals;
  return 0;
}