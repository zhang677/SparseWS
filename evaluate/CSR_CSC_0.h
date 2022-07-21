#include "../utils/dataloader.h"
#include "../utils/lib.h"

void compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
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

  int32_t* restrict C2_nnz = (int32_t* )calloc(A1_dimension, sizeof(int32_t));

  #pragma omp parallel for schedule(runtime)
  for (int32_t qi = 0; qi < A1_dimension; qi++) {
    for (int32_t qk = 0; qk < B2_dimension; qk++) {
      bool qtqjC_val = 0;
      int32_t qjA = A2_pos[qi];
      int32_t pA2_end = A2_pos[(qi + 1)];
      int32_t qjB = B2_pos[qk];
      int32_t pB2_end = B2_pos[(qk + 1)];

      while (qjA < pA2_end && qjB < pB2_end) {
        int32_t qjA0 = A2_crd[qjA];
        int32_t qjB0 = B2_crd[qjB];
        int32_t qj = TACO_MIN(qjA0,qjB0);
        if (qjA0 == qj && qjB0 == qj) {
          qtqjC_val = 1;
        }
        qjA += (int32_t)(qjA0 == qj);
        qjB += (int32_t)(qjB0 == qj);
      }
      C2_nnz[qi] = C2_nnz[qi] + (int32_t)qtqjC_val;
    }
  }
  
  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t i = 0; i < C1_dimension; i++) {
    C2_pos[i + 1] = C2_pos[i] + C2_nnz[i];
  }

  C2_crd = (int32_t*)malloc(sizeof(int32_t) * C2_pos[C1_dimension]);
  C_vals = (float* )calloc(C2_pos[C1_dimension], sizeof(float));

  #pragma omp parallel for schedule(runtime)
  for (int32_t i = 0; i < A1_dimension; i++) {
    int32_t pC2 = C2_pos[i];
    for (int32_t k = 0; k < B2_dimension; k++) {
      bool tqjC_val = 0;
      int32_t jA = A2_pos[i];
      int32_t pA2_end0 = A2_pos[(i + 1)];
      int32_t jB = B2_pos[k];
      int32_t pB2_end0 = B2_pos[(k + 1)];

      while (jA < pA2_end0 && jB < pB2_end0) {
        int32_t jA0 = A2_crd[jA];
        int32_t jB0 = B2_crd[jB];
        int32_t j = TACO_MIN(jA0,jB0);
        if (jA0 == j && jB0 == j) {
          //int32_t pC2 = C2_pos[i];
          //C2_pos[i] = C2_pos[i] + 1;
          C2_crd[pC2] = k;
          C_vals[pC2] = C_vals[pC2] + A_vals[jA] * B_vals[jB];
          tqjC_val = 1;
        }
        jA += (int32_t)(jA0 == j);
        jB += (int32_t)(jB0 == j);
      }
      pC2 += (int32_t)tqjC_val;
    }
  }
  

  free(C2_nnz);

  C->indices[1][0] = (int32_t*)(C2_pos);
  C->indices[1][1] = (int32_t*)(C2_crd);
  C->vals = (float*)C_vals;
}

void CSR_CSC_0(const string& A_name, const string& B_name) {
    vector<int> indptr;
    vector<int> indices;
    vector<int> id_buffer;
    vector<float> value;
    int nrow;
    int ncol;
    int nnz;
    read_mtx_csr(A_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t A = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    read_mtx_csc(B_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{1,0});
    taco_tensor_t C;
    init_taco_tensor_DC(&C, nrow, ncol, {0,1});
    compute(&C,&A,&B);
    print_taco_tensor_DC(&A);
    print_taco_tensor_DC(&B);
    print_taco_tensor_DC(&C);
}