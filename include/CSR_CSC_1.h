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

  int32_t* restrict C2_nnz = (int32_t*)calloc(A1_dimension, sizeof(int32_t));

  bool* restrict qw_all = 0;
  int32_t* restrict qw_index_list_all = 0;
  qw_index_list_all = (int32_t*)malloc(sizeof(int32_t) * (B2_dimension * omp_get_max_threads()));
  bool* restrict qw_already_set_all = (bool*)calloc((B2_dimension * omp_get_max_threads()), sizeof(bool));
  qw_all = (bool*)malloc(sizeof(bool) * (B2_dimension * omp_get_max_threads()));

  #pragma omp parallel for schedule(runtime)
  for (int32_t qi = 0; qi < A1_dimension; qi++) {
    int32_t qw_index_list_size = 0;
    bool* restrict qw = qw_all + B2_dimension * omp_get_thread_num();
    int32_t* restrict qw_index_list = qw_index_list_all + B2_dimension * omp_get_thread_num();
    bool* restrict qw_already_set = qw_already_set_all + B2_dimension * omp_get_thread_num();
    for (int32_t qk = 0; qk < B2_dimension; qk++) {
      bool qtqjw_val = 0;
      int32_t qjA = A2_pos[qi];
      int32_t pA2_end = A2_pos[(qi + 1)];
      int32_t qjB = B2_pos[qk];
      int32_t pB2_end = B2_pos[(qk + 1)];

      while (qjA < pA2_end && qjB < pB2_end) {
        int32_t qjA0 = A2_crd[qjA];
        int32_t qjB0 = B2_crd[qjB];
        int32_t qj = TACO_MIN(qjA0,qjB0);
        if (qjA0 == qj && qjB0 == qj) {
          qtqjw_val = 1;
        }
        qjA += (int32_t)(qjA0 == qj);
        qjB += (int32_t)(qjB0 == qj);
      }
      if (!qw_already_set[qk]) {
        qw[qk] = qtqjw_val;
        qw_index_list[qw_index_list_size] = qk;
        qw_already_set[qk] = 1;
        qw_index_list_size++;
      }
      else {
        qw[qk] = qw[qk] || qtqjw_val;
      }
    }
    for (int32_t qw_index_locator = 0; qw_index_locator < qw_index_list_size; qw_index_locator++) {
      int32_t qk = qw_index_list[qw_index_locator];
      C2_nnz[qi] = C2_nnz[qi] + (int32_t)qw[qk];
      qw_already_set[qk] = 0;
    }
  }

  free(qw_index_list_all);
  free(qw_already_set_all);
  free(qw_all);

  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t i = 0; i < C1_dimension; i++) {
    C2_pos[i + 1] = C2_pos[i] + C2_nnz[i];
  }
  C2_crd = (int32_t*)malloc(sizeof(int32_t) * C2_pos[C1_dimension]);
  C_vals = (float*)malloc(sizeof(float) * C2_pos[C1_dimension]);

  float* restrict w_all = 0;
  int32_t* restrict w_index_list_all = 0;
  w_index_list_all = (int32_t*)malloc(sizeof(int32_t) * (B2_dimension * omp_get_max_threads()));
  bool* restrict w_already_set_all = (bool*)calloc((B2_dimension * omp_get_max_threads()), sizeof(bool));
  w_all = (float*)malloc(sizeof(float) * (B2_dimension * omp_get_max_threads()));

  #pragma omp parallel for schedule(runtime)
  for (int32_t i = 0; i < A1_dimension; i++) {
    int32_t w_index_list_size = 0;
    float* restrict w = w_all + B2_dimension * omp_get_thread_num();
    int32_t* restrict w_index_list = w_index_list_all + B2_dimension * omp_get_thread_num();
    bool* restrict w_already_set = w_already_set_all + B2_dimension * omp_get_thread_num();
    for (int32_t k = 0; k < B2_dimension; k++) {
      int32_t jA = A2_pos[i];
      int32_t pA2_end0 = A2_pos[(i + 1)];
      int32_t jB = B2_pos[k];
      int32_t pB2_end0 = B2_pos[(k + 1)];

      while (jA < pA2_end0 && jB < pB2_end0) {
        int32_t jA0 = A2_crd[jA];
        int32_t jB0 = B2_crd[jB];
        int32_t j = TACO_MIN(jA0,jB0);
        if (jA0 == j && jB0 == j) {
          if (!w_already_set[k]) {
            w[k] = A_vals[jA] * B_vals[jB];
            w_index_list[w_index_list_size] = k;
            w_already_set[k] = 1;
            w_index_list_size++;
          }
          else {
            w[k] = w[k] + A_vals[jA] * B_vals[jB];
          }
        }
        jA += (int32_t)(jA0 == j);
        jB += (int32_t)(jB0 == j);
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

  free(w_index_list_all);
  free(w_already_set_all);
  free(w_all);

  for (int32_t p = 0; p < C1_dimension; p++) {
    C2_pos[C1_dimension - p] = C2_pos[((C1_dimension - p) - 1)];
  }
  C2_pos[0] = 0;

  free(C2_nnz);

  C->indices[1][0] = (int32_t*)(C2_pos);
  C->indices[1][1] = (int32_t*)(C2_crd);
  C->vals = (float*)C_vals;
}

void CSR_CSC_1(const string& A_name, const string& B_name) {
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