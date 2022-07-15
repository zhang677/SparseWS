#ifndef SPARSE_DOT_H
#define SPARSE_DOT_H
#include "../utils/dataloader.h"
#include "../utils/lib.h"


void benchmark (taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C) {
  int* restrict A1_pos = (int*)(A->indices[0][0]);
  int* restrict A1_crd = (int*)(A->indices[0][1]);
  float* restrict A_vals = (float*)(A->vals);
  int* restrict B1_pos = (int*)(B->indices[0][0]);
  int* restrict B1_crd = (int*)(B->indices[0][1]);
  float* restrict B_vals = (float*)(B->vals);
  int* restrict C1_pos = (int*)(C->indices[0][0]);
  int* restrict C1_crd = (int*)(C->indices[0][1]);
  float* restrict C_vals = (float*)(C->vals);

  A1_pos = (int32_t*)malloc(sizeof(int32_t) * 2);
  A1_pos[0] = 0;
  int32_t A1_crd_size = 1048576;
  A1_crd = (int32_t*)malloc(sizeof(int32_t) * A1_crd_size);
  int32_t iA = 0;
  int32_t A_capacity = 1048576;
  A_vals = (float*)malloc(sizeof(float) * A_capacity);


  int32_t iB = B1_pos[0];
  int32_t pB1_end = B1_pos[1];
  int32_t iC = C1_pos[0];
  int32_t pC1_end = C1_pos[1];

  while (iB < pB1_end && iC < pC1_end) {
    int32_t iB0 = B1_crd[iB];
    int32_t iC0 = C1_crd[iC];
    int32_t i = TACO_MIN(iB0,iC0);
    if (iB0 == i && iC0 == i) {
      if (A_capacity <= iA) {
        A_vals = (float*)realloc(A_vals, sizeof(float) * (A_capacity * 2));
        A_capacity *= 2;
      }
      A_vals[iA] = B_vals[iB] * C_vals[iC];
      if (A1_crd_size <= iA) {
        A1_crd = (int32_t*)realloc(A1_crd, sizeof(int32_t) * (A1_crd_size * 2));
        A1_crd_size *= 2;
      }
      A1_crd[iA] = i;
      iA++;
    }
    iB += (int32_t)(iB0 == i);
    iC += (int32_t)(iC0 == i);
  }

  A1_pos[1] = iA;

  A->indices[0][0] = (int32_t*)(A1_pos);
  A->indices[0][1] = (int32_t*)(A1_crd);
  A->vals = (float*)A_vals;
}

void split(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C){
  int* restrict A1_pos = (int*)(A->indices[0][0]);
  int* restrict A1_crd = (int*)(A->indices[0][1]);
  float* restrict A_vals = (float*)(A->vals);
  int* restrict B1_pos = (int*)(B->indices[0][0]);
  int* restrict B1_crd = (int*)(B->indices[0][1]);
  float* restrict B_vals = (float*)(B->vals);
  int C1_dimension = (int)(C->dimensions[1]);
  int* restrict C1_pos = (int*)(C->indices[0][0]);
  int* restrict C1_crd = (int*)(C->indices[0][1]);
  float* restrict C_vals = (float*)(C->vals);
  A1_pos = (int32_t*)malloc(sizeof(int32_t) * 2);
  A1_pos[0] = 0;
  int32_t A1_crd_size = 1048576;
  A1_crd = (int32_t*)malloc(sizeof(int32_t) * A1_crd_size);
  int32_t iA = 0;
  int32_t A_capacity = 1048576;
  A_vals = (float*)malloc(sizeof(float) * A_capacity);

  for (int32_t io = 0; io < ((C1_dimension + 3) / 4); io++) {
    int32_t pB1_begin = io * 4;
    int32_t iB = taco_binarySearchAfter(B1_crd, B1_pos[0], B1_pos[1], pB1_begin);
    int32_t pB1_end = B1_pos[1];
    int32_t pC1_begin = io * 4;
    int32_t iC = taco_binarySearchAfter(C1_crd, C1_pos[0], C1_pos[1], pC1_begin);
    int32_t pC1_end = C1_pos[1];
    int32_t iB0 = B1_crd[iB];
    int32_t iC0 = C1_crd[iC];
    int32_t i = TACO_MIN(iB0,iC0);
    int32_t ii = i - io * 4;
    int32_t ii_end = 4;

    while ((iB < pB1_end && ii < ii_end) && iC < pC1_end) {
      iB0 = B1_crd[iB];
      iC0 = C1_crd[iC];
      i = TACO_MIN(iB0,iC0);
      if (iB0 == i && iC0 == i) {
        if (A_capacity <= iA) {
          A_vals = (float*)realloc(A_vals, sizeof(float) * (A_capacity * 2));
          A_capacity *= 2;
        }
        A_vals[iA] = B_vals[iB] * C_vals[iC];
        if (A1_crd_size <= iA) {
          A1_crd = (int32_t*)realloc(A1_crd, sizeof(int32_t) * (A1_crd_size * 2));
          A1_crd_size *= 2;
        }
        A1_crd[iA] = i;
        iA++;
      }
      iB += (int32_t)(iB0 == i);
      iC += (int32_t)(iC0 == i);
      iB0 = B1_crd[iB];
      iC0 = C1_crd[iC];
      i = TACO_MIN(iB0,iC0);
      ii = i - io * 4;
    }

    A1_pos[1] = iA;
  }
  A->indices[0][0] = (int32_t*)(A1_pos);
  A->indices[0][1] = (int32_t*)(A1_crd);
  A->vals = (float*)A_vals;
}

void compute(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C) {
  int* restrict A1_pos = (int*)(A->indices[0][0]);
  int* restrict A1_crd = (int*)(A->indices[0][1]);
  float* restrict A_vals = (float*)(A->vals);
  int* restrict B1_pos = (int*)(B->indices[0][0]);
  int* restrict B1_crd = (int*)(B->indices[0][1]);
  float* restrict B_vals = (float*)(B->vals);
  int C1_dimension = (int)(C->dimensions[1]);
  int* restrict C1_pos = (int*)(C->indices[0][0]);
  int* restrict C1_crd = (int*)(C->indices[0][1]);
  float* restrict C_vals = (float*)(C->vals);

  A1_pos = (int32_t*)malloc(sizeof(int32_t) * 2);
  A1_pos[0] = 0;
  int32_t A1_crd_size = 1048576;
  A1_crd = (int32_t*)malloc(sizeof(int32_t) * A1_crd_size);
  int32_t iA = 0;
  int32_t A_capacity = 1048576;
  A_vals = (float*)malloc(sizeof(float) * A_capacity);

  int32_t ws_index_list_size = 0;
  float* restrict ws = 0;
  int32_t* restrict ws_index_list = 0;
  ws_index_list = (int32_t*)malloc(sizeof(int32_t) * C1_dimension);
  bool* restrict ws_already_set = calloc(C1_dimension, sizeof(bool));
  ws = (float*)malloc(sizeof(float) * C1_dimension);
  for (int32_t io = 0; io < ((C1_dimension + 3) / 4); io++) {
    int32_t pB1_begin = io * 4;
    int32_t iB = taco_binarySearchAfter(B1_crd, B1_pos[0], B1_pos[1], pB1_begin);
    int32_t pB1_end = B1_pos[1];
    int32_t pC1_begin = io * 4;
    int32_t iC = taco_binarySearchAfter(C1_crd, C1_pos[0], C1_pos[1], pC1_begin);
    int32_t pC1_end = C1_pos[1];
    int32_t iB0 = B1_crd[iB];
    int32_t iC0 = C1_crd[iC];
    int32_t i = TACO_MIN(iB0,iC0);
    int32_t ii = i - io * 4;
    int32_t ii_end = 4;

    while ((iB < pB1_end && ii < ii_end) && iC < pC1_end) {
      iB0 = B1_crd[iB];
      iC0 = C1_crd[iC];
      i = TACO_MIN(iB0,iC0);
      if (iB0 == i && iC0 == i) {
        if (!ws_already_set[i]) {
          ws[i] = B_vals[iB] * C_vals[iC];
          ws_index_list[ws_index_list_size] = i;
          ws_already_set[i] = 1;
          ws_index_list_size++;
        }
        else {
          ws[i] = B_vals[iB] * C_vals[iC];
        }
      }
      iB += (int32_t)(iB0 == i);
      iC += (int32_t)(iC0 == i);
      iB0 = B1_crd[iB];
      iC0 = C1_crd[iC];
      i = TACO_MIN(iB0,iC0);
      ii = i - io * 4;
    }
  }
  qsort(ws_index_list, ws_index_list_size, sizeof(int32_t), cmp);
  for (int32_t io = 0; io < ((C1_dimension + 3) / 4); io++) {
    for (int32_t ws_index_locator = 0; ws_index_locator < ws_index_list_size; ws_index_locator++) {
      int32_t ii = ws_index_list[ws_index_locator];
      int32_t i = io * 4 + ii;
      if (i >= C1_dimension)
        continue;

      if (A_capacity <= iA) {
        A_vals = (float*)realloc(A_vals, sizeof(float) * (A_capacity * 2));
        A_capacity *= 2;
      }
      A_vals[iA] = ws[i];
      if (A1_crd_size <= iA) {
        A1_crd = (int32_t*)realloc(A1_crd, sizeof(int32_t) * (A1_crd_size * 2));
        A1_crd_size *= 2;
      }
      A1_crd[iA] = ii;
      iA++;
      ws_already_set[ii] = 0;
    }

    A1_pos[1] = iA;
  }
  free(ws_index_list);
  free(ws_already_set);
  free(ws);
  A->indices[0][0] = (int32_t*)(A1_pos);
  A->indices[0][1] = (int32_t*)(A1_crd);
  A->vals = (float*)A_vals;
}

void modified_v1(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C) {
  int* restrict A1_pos = (int*)(A->indices[0][0]);
  int* restrict A1_crd = (int*)(A->indices[0][1]);
  float* restrict A_vals = (float*)(A->vals);
  int* restrict B1_pos = (int*)(B->indices[0][0]);
  int* restrict B1_crd = (int*)(B->indices[0][1]);
  float* restrict B_vals = (float*)(B->vals);
  int C1_dimension = (int)(C->dimensions[1]);
  int* restrict C1_pos = (int*)(C->indices[0][0]);
  int* restrict C1_crd = (int*)(C->indices[0][1]);
  float* restrict C_vals = (float*)(C->vals);

  A1_pos = (int32_t*)malloc(sizeof(int32_t) * 2);
  A1_pos[0] = 0;
  int32_t A1_crd_size = 1048576;
  A1_crd = (int32_t*)malloc(sizeof(int32_t) * A1_crd_size);
  int32_t iA = 0;
  int32_t A_capacity = 1048576;
  A_vals = (float*)malloc(sizeof(float) * A_capacity);

  int32_t ws_index_list_size = 0;
  float* restrict ws = 0;
  int32_t* restrict ws_index_list = 0;
  ws_index_list = (int32_t*)malloc(sizeof(int32_t) * C1_dimension);
  bool* restrict ws_already_set = calloc(C1_dimension, sizeof(bool));
  ws = (float*)malloc(sizeof(float) * C1_dimension);
  for (int32_t io = 0; io < ((C1_dimension + 3) / 4); io++) {
    int32_t pB1_begin = io * 4;
    int32_t iB = taco_binarySearchAfter(B1_crd, B1_pos[0], B1_pos[1], pB1_begin);
    int32_t pB1_end = B1_pos[1];
    int32_t pC1_begin = io * 4;
    int32_t iC = taco_binarySearchAfter(C1_crd, C1_pos[0], C1_pos[1], pC1_begin);
    int32_t pC1_end = C1_pos[1];
    int32_t iB0 = B1_crd[iB];
    int32_t iC0 = C1_crd[iC];
    int32_t i = TACO_MIN(iB0,iC0);
    int32_t ii = i - io * 4;
    int32_t ii_end = 4;

    while ((iB < pB1_end && ii < ii_end) && iC < pC1_end) {
      iB0 = B1_crd[iB];
      iC0 = C1_crd[iC];
      i = TACO_MIN(iB0,iC0);
      if (iB0 == i && iC0 == i) {
        if (!ws_already_set[i]) {
          ws[i] = B_vals[iB] * C_vals[iC];
          ws_index_list[ws_index_list_size] = i;
          ws_already_set[i] = 1;
          ws_index_list_size++;
        }
        else {
          ws[i] = B_vals[iB] * C_vals[iC];
        }
      }
      iB += (int32_t)(iB0 == i);
      iC += (int32_t)(iC0 == i);
      iB0 = B1_crd[iB];
      iC0 = C1_crd[iC];
      i = TACO_MIN(iB0,iC0);
      ii = i - io * 4;
    }
  }
  qsort(ws_index_list, ws_index_list_size, sizeof(int32_t), cmp);
  for (int32_t io = 0; io < ((C1_dimension + 3) / 4); io++) {
    int32_t p_begin = io * 4;
    int32_t p_end = (io + 1) * 4;
    int32_t ws_index_locator = taco_binarySearchAfter(ws_index_list, 0, ws_index_list_size, p_begin);
    int32_t i =  ws_index_list[ws_index_locator];

    while(i < p_end && ws_index_locator < ws_index_list_size) {
      i =  ws_index_list[ws_index_locator];
      ws_index_locator = ws_index_locator + 1;
      if (i >= p_end)
        break;
      if (A_capacity <= iA) {
        A_vals = (float*)realloc(A_vals, sizeof(float) * (A_capacity * 2));
        A_capacity *= 2;
      }
      A_vals[iA] = ws[i];
      if (A1_crd_size <= iA) {
        A1_crd = (int32_t*)realloc(A1_crd, sizeof(int32_t) * (A1_crd_size * 2));
        A1_crd_size *= 2;
      }
      A1_crd[iA] = i;
      iA++;
      A1_pos[1] = iA;
    }
  }
  free(ws_index_list);
  free(ws_already_set);
  free(ws);
  A->indices[0][0] = (int32_t*)(A1_pos);
  A->indices[0][1] = (int32_t*)(A1_crd);
  A->vals = (float*)A_vals;
}

void modified_v2(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C) {
  int* restrict A1_pos = (int*)(A->indices[0][0]);
  int* restrict A1_crd = (int*)(A->indices[0][1]);
  float* restrict A_vals = (float*)(A->vals);
  int* restrict B1_pos = (int*)(B->indices[0][0]);
  int* restrict B1_crd = (int*)(B->indices[0][1]);
  float* restrict B_vals = (float*)(B->vals);
  int C1_dimension = (int)(C->dimensions[1]);
  int* restrict C1_pos = (int*)(C->indices[0][0]);
  int* restrict C1_crd = (int*)(C->indices[0][1]);
  float* restrict C_vals = (float*)(C->vals);

  A1_pos = (int32_t*)malloc(sizeof(int32_t) * 2);
  A1_pos[0] = 0;
  int32_t A1_crd_size = 1048576;
  A1_crd = (int32_t*)malloc(sizeof(int32_t) * A1_crd_size);
  int32_t iA = 0;
  int32_t A_capacity = 1048576;
  A_vals = (float*)malloc(sizeof(float) * A_capacity);

  int32_t ws_index_list_size = 0;
  float* restrict ws = 0;
  int32_t* restrict ws_index_list = 0;
  ws_index_list = (int32_t*)malloc(sizeof(int32_t) * C1_dimension);
  bool* restrict ws_already_set = calloc(C1_dimension, sizeof(bool));
  ws = (float*)malloc(sizeof(float) * C1_dimension);
  for (int32_t io = 0; io < ((C1_dimension + 3) / 4); io++) {
    int32_t pB1_begin = io * 4;
    int32_t iB = taco_binarySearchAfter(B1_crd, B1_pos[0], B1_pos[1], pB1_begin);
    int32_t pB1_end = B1_pos[1];
    int32_t pC1_begin = io * 4;
    int32_t iC = taco_binarySearchAfter(C1_crd, C1_pos[0], C1_pos[1], pC1_begin);
    int32_t pC1_end = C1_pos[1];
    int32_t iB0 = B1_crd[iB];
    int32_t iC0 = C1_crd[iC];
    int32_t i = TACO_MIN(iB0,iC0);
    int32_t ii = i - io * 4;
    int32_t ii_end = 4;

    while ((iB < pB1_end && ii < ii_end) && iC < pC1_end) {
      iB0 = B1_crd[iB];
      iC0 = C1_crd[iC];
      i = TACO_MIN(iB0,iC0);
      if (iB0 == i && iC0 == i) {
        if (!ws_already_set[i]) {
          ws[i] = B_vals[iB] * C_vals[iC];
          ws_index_list[ws_index_list_size] = i;
          ws_already_set[i] = 1;
          ws_index_list_size++;
        }
        else {
          ws[i] = B_vals[iB] * C_vals[iC];
        }
      }
      iB += (int32_t)(iB0 == i);
      iC += (int32_t)(iC0 == i);
      iB0 = B1_crd[iB];
      iC0 = C1_crd[iC];
      i = TACO_MIN(iB0,iC0);
      ii = i - io * 4;
    }
  }
  qsort(ws_index_list, ws_index_list_size, sizeof(int32_t), cmp);
  for (int32_t io = 0; io < ((C1_dimension + 3) / 4); io++) {
    for (int32_t ii = 0; ii < 4; ii++){
      int32_t i = ii + io * 4;
      if (i >= C1_dimension) {
        continue;
      }
      if (!ws_already_set[i]) {
        continue;
      }
      if (A_capacity <= iA) {
        A_vals = (float*)realloc(A_vals, sizeof(float) * (A_capacity * 2));
        A_capacity *= 2;
      }
      A_vals[iA] = ws[i];
      if (A1_crd_size <= iA) {
        A1_crd = (int32_t*)realloc(A1_crd, sizeof(int32_t) * (A1_crd_size * 2));
        A1_crd_size *= 2;
      }
      A1_crd[iA] = i;
      iA++;
      A1_pos[1] = iA;
    }
  }
  free(ws_index_list);
  free(ws_already_set);
  free(ws);
  A->indices[0][0] = (int32_t*)(A1_pos);
  A->indices[0][1] = (int32_t*)(A1_crd);
  A->vals = (float*)A_vals;
}

void test_sparse_dot() {
    int N = 16;
    const COO B_coo = gen_sparse_vector(N);
    EigenCSR B_eigen = to_eigen_csr(B_coo);
    taco_tensor_t B = to_taco_tensor(B_eigen);
    const COO C_coo = gen_sparse_vector(N);
    EigenCSR C_eigen = to_eigen_csr(C_coo);
    taco_tensor_t C = to_taco_tensor(C_eigen);
    taco_tensor_t A = get_csr_taco_tensor(1,N);
    modified_v2(&A,&B,&C);
    taco_tensor_t A_ref = get_csr_taco_tensor(1,N);
    benchmark(&A_ref,&B,&C);
    print_csr_vector(&A);
    print_csr_vector(&A_ref);
}
#endif