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

  int32_t* restrict C2_nnz = (int32_t*)calloc(A1_dimension, sizeof(int32_t));

  int32_t* restrict qw_index_list_all = 0;
  qw_index_list_all = (int32_t*)malloc(sizeof(int32_t) * (B2_dimension * omp_get_max_threads()));
  bool* restrict qw_already_set_all = (bool*)calloc((B2_dimension * omp_get_max_threads()), sizeof(bool));

  #pragma omp parallel for schedule(runtime)
  for (int32_t qi = 0; qi < A1_dimension; qi++) {
    int32_t qw_index_list_size = 0;
    int32_t* restrict qw_index_list = qw_index_list_all + B2_dimension * omp_get_thread_num();
    bool* restrict qw_already_set = qw_already_set_all + B2_dimension * omp_get_thread_num();
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

  free(qw_index_list_all);
  free(qw_already_set_all);

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
  return 0;
}

const int w_psize = 1;
struct w_pair{
  int32_t pos[1];
  float val;
} ;
// cmp function is generated when knowing the dimension of pos in w_pair
int esc_cmp(const void *b, const void *a) {
  const int psize = w_psize;
  for (int i=0; i<psize; i++) {
    if (((w_pair*)b)->pos[i] == ((w_pair*)a)->pos[i]) continue;
    return (((w_pair*)b)->pos[i] - ((w_pair*)a)->pos[i]);
  }
  return (((w_pair*)b)->pos[psize-1] - ((w_pair*)a)->pos[psize-1]);
}
bool pos_equal(const w_pair *a, const w_pair *b, const int32_t size) {
  for (int i=0; i<size; i++) {
    if(a->pos[i]!=b->pos[i]) return false;
  }
  return true;
}

void ESC_assemble_compute_v1(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
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

  // With assemble
  int32_t* restrict C2_nnz_upper_bound = (int32_t*)malloc(sizeof(int32_t) * A1_dimension);
  int32_t* restrict C2_nnz = (int32_t*)malloc(sizeof(int32_t)*A1_dimension);
  int32_t* qw_index_list = (int32_t*)malloc(sizeof(int32_t) * B2_dimension);
  bool* restrict qw_already_set = (bool*)malloc(sizeof(bool) * B2_dimension);

  for (int32_t qi = 0; qi < A1_dimension ; qi++) {
    int32_t qw_index_list_size = 0;
    C2_nnz_upper_bound[qi] = 0;
    for (int32_t qjA = A2_pos[qi]; qjA < A2_pos[(qi+1)]; qjA++) {
      int32_t qj = A2_crd[qjA];
      for (int32_t qkB = B2_pos[qj]; qkB < B2_pos[(qj + 1)]; qkB++) {
        int32_t qk = B2_crd[qkB];
        if (!qw_already_set[qk]) {
          qw_index_list[qw_index_list_size] = qk;
          qw_already_set[qk] = 1;
          qw_index_list_size++;
        }
      }
      C2_nnz_upper_bound[qi] += (B2_pos[(qj+1)] - B2_pos[qj]);
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
  free(C2_nnz);
  C2_crd = (int32_t*)malloc(sizeof(int32_t) * C2_pos[C1_dimension]);
  C_vals = (float*)malloc(sizeof(float) * C2_pos[C1_dimension]);

  w_pair* restrict w = 0;
  w = (w_pair*)malloc(sizeof(w_pair) * C2_nnz_upper_bound[0]);
  int32_t w_size = C2_nnz_upper_bound[0];
  for(int32_t i = 0; i<A1_dimension; i++) {
    if(C2_nnz_upper_bound[i]==0) continue;
    if(C2_nnz_upper_bound[i]>w_size) {
      free(w);
      w = (w_pair*)realloc(w, sizeof(w_pair) * C2_nnz_upper_bound[i]);
    }
    w_size = C2_nnz_upper_bound[i];
    int32_t iw = 0;
    for (int32_t jA = A2_pos[i]; jA < A2_pos[(i+1)]; jA++) {
      int32_t j = A2_crd[jA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[(j+1)]; kB++) {
        w[iw].pos[0] = B2_crd[kB]; 
        w[iw].val = A_vals[jA] * B_vals[kB];
        iw ++;
      }
    }
    qsort(w, w_size, sizeof(w_pair), esc_cmp);
    w_pair current = w[0];
    int32_t pC2 = C2_pos[i];
    for (int32_t w_index_locator = 0; w_index_locator < w_size; w_index_locator++) {
      if(!pos_equal(&w[w_index_locator],&current,1)){
        current = w[w_index_locator];
        pC2 ++;
      }
      C2_crd[pC2] = w[w_index_locator].pos[0];
      C_vals[pC2] += w[w_index_locator].val;
    }
  }

  free(w);
  C->indices[1][0] = (int32_t*)(C2_pos);
  C->indices[1][1] = (int32_t*)(C2_crd);
  C->vals = (float*)C_vals;

}




void CSR_CSR_1(const string& A_name, const string& B_name) {
    vector<int> indptr;
    vector<int> indices;
    vector<int> id_buffer;
    vector<float> value;
    int nrow;
    int ncol;
    int nnz;
    read_mtx_csr(A_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t A = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    read_mtx_csr(B_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    taco_tensor_t C;
    init_taco_tensor_DC(&C, nrow, ncol, {0,1});
    //ESC_assemble_compute_v1(&C,&A,&B);
    compute(&C,&A,&B);
    print_taco_tensor_DC(&A);
    print_taco_tensor_DC(&B);
    print_taco_tensor_DC(&C);
}