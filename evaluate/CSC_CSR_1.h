#include "../utils/dataloader.h"
#include "../utils/lib.h"

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

  int32_t restrict C2_nnz = (int32_t*)calloc(A1_dimension, sizeof(int32_t));

  int32_t* restrict qw_index_list_all = 0;
  qw_index_list_all = (int32_t*)malloc(sizeof(int32_t) * (B2_dimension * B1_dimension * omp_get_max_threads()));
  bool* restrict qw_already_set_all = (bool*)calloc((B2_dimension * B1_dimension * omp_get_max_threads()), sizeof(bool));

  #pragma omp parallel for schedule(runtime)
  int32_t qw_index_list_size = 0;
  int32_t* restrict qw_index_list = qw_index_list_all + (B2_dimension * B1_dimension * omp_get_thread_num());
  bool* restrict qw_already_set = qw_already_set_all + (B2_dimension * B1_dimension * omp_get_thread_num());
  for (int32_t qj = 0; qj < B1_dimension; qj++) {
    for (int32_t qiA = A2_pos[qj]; qiA < A2_pos[(qj + 1)]; qiA++) {
      int32_t qi = A2_crd[qiA];
      for (int32_t qkB = B2_pos[qj]; qkB < B2_pos[(qj + 1)]; qkB++) {
        int32_t qk = B2_crd[qkB];
        if (!qw_already_set[(qi * B2_dimension) + qk]) {
          qw_already_set[(qi * B2_dimension) + qk] = true;
          qw_index_list[qw_index_list_size] = (qi * B2_dimension) + qk;
          qw_index_list_size++;
        }
      }
    }
  }
  for (int32_t qw_index_list_i = 0; qw_index_list_i < qw_index_list_size; qw_index_list_i++) {
    int32_t qw = qw_index_list[qw_index_list_i];
    int32_t qi = qw / B2_dimension;
    int32_t qk = qw % B2_dimension;
    for (int32_t qj = 0; qj < B1_dimension; qj++) {
      for (int32_t qiA = A2_pos[qj]; qiA < A2_pos[(qj + 1)]; qiA++) {
        if (qi == A2_crd[qiA]) {
          for (int32_t qkB = B2_pos[qj]; qkB < B2_pos[(qj + 1)]; qkB++) {
            if (qk == B2_crd[qkB]) {
              C2_nnz[qi]++;
            }
          }
        }
      }
    }
  }
}