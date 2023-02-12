#include "../utils/dataloader.h"
#include "../utils/lib.h"

const int w_order = 2;
struct wspace {
  int32_t pos[2];
  float val;
};
// cmp function is generated when knowing the dimension of pos in wspace
int esc_cmp(const void *b, const void *a) {
  const int psize = w_order;
  for (int i=0; i<psize; i++) {
    if (((wspace*)b)->pos[i] == ((wspace*)a)->pos[i]) continue;
    return (((wspace*)b)->pos[i] - ((wspace*)a)->pos[i]);
  }
  return (((wspace*)b)->pos[psize-1] - ((wspace*)a)->pos[psize-1]);
}

int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
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

  int32_t w_accumulator_capacity = 3; // Arbitrary
  int32_t w_accumulator_size = 0;
  wspace* restrict w_accumulator = 0;
  w_accumulator = (wspace*)malloc(sizeof(wspace) * w_accumulator_capacity);
  int32_t w_all_capacity = 2; // Arbitrary (larger than 1)
  int32_t w_all_size = 0;
  wspace* restrict w_all = 0;
  w_all = (wspace*)malloc(sizeof(wspace) * w_all_capacity);

  for (int32_t i = 0; i < A1_dimension; i++) {
    for (int32_t k = 0; k < B2_dimension; k++) {
      int32_t jA = A2_pos[i];
      int32_t pA2_end = A2_pos[(i + 1)];
      int32_t jB = B2_pos[k];
      int32_t pB2_end = B2_pos[(k + 1)];

      while (jA < pA2_end && jB < pB2_end) {
        int32_t jA0 = A2_crd[jA];
        int32_t jB0 = B2_crd[jB];
        int32_t j = TACO_MIN(jA0, jB0);
        if (jA0 == j && jB0 == j) {
            w_accumulator[w_accumulator_size].pos[0] = i;
            w_accumulator[w_accumulator_size].pos[1] = k; //pos[2] => coord
            w_accumulator[w_accumulator_size].val = A_vals[jA] * B_vals[jB];
            w_accumulator_size ++;
            if (w_accumulator_size >= w_accumulator_capacity) {
                qsort(w_accumulator, w_accumulator_size, sizeof(wspace), esc_cmp);
                int32_t begin = (w_all_size == 0);
                if (w_all_size == 0) {
                    w_all_size ++;
                    w_all[w_all_size-1] = w_accumulator[0];
                }
                for (int32_t wi = begin; wi < w_accumulator_size; wi++) {
                    if(esc_cmp(&w_accumulator[wi],&w_all[w_all_size-1]) == 0) {
                        w_all[w_all_size-1].val += w_accumulator[wi].val;
                    }
                    else {
                        if(w_all_size >= w_all_capacity) {
                            w_all_capacity = w_all_capacity * 2;
                            w_all = (wspace*)realloc(w_all, sizeof(wspace) * w_all_capacity);
                        }
                        w_all_size ++;
                        w_all[w_all_size-1] = w_accumulator[wi];
                    }
                }
                w_accumulator_size = 0;
            }
        }
        jA += (int32_t)(jA0 == j);
        jB += (int32_t)(jB0 == j);
      }
    }
  }
  if(w_accumulator_size > 0) {
    qsort(w_accumulator, w_accumulator_size, sizeof(wspace), esc_cmp);
    for (int32_t wi = 0; wi < w_accumulator_size; wi++) {
        if(esc_cmp(&w_accumulator[wi],&w_all[w_all_size-1]) == 0) {
            w_all[w_all_size-1].val += w_accumulator[wi].val;
        }
        else {
            if(w_all_size >= w_all_capacity) {
                w_all_capacity = w_all_capacity * 2;
                w_all = (wspace*)realloc(w_all, sizeof(wspace) * w_all_capacity);
            }
            w_all_size ++;
            w_all[w_all_size-1] = w_accumulator[wi];
        }
    }
    w_accumulator_size = 0;
  }
  int  C_COO1_pos[2] = {0,w_all_size};
  int32_t* C_COO_crd[2];
  for (int wi = 0; wi < w_order; wi++) {
    C_COO_crd[wi] = (int32_t*)malloc(sizeof(int32_t) * w_all_size);
    for (int wk = 0; wk < w_all_size; wk++) {
        C_COO_crd[wi][wk] = w_all[wk].pos[wi];
    }
  }
  float*  C_COO_vals = (float*)malloc(sizeof(float) * w_all_size);
  for (int wk = 0; wk < w_all_size; wk++) {
    C_COO_vals[wk] = w_all[wk].val;
  }
  pack_C(C, C_COO1_pos, C_COO_crd[0], C_COO_crd[1], C_COO_vals);
  return 0;
}

void CSR_CSC_4(const string& A_name, const string& B_name, taco_tensor_t* C, bool print = false) {
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
    init_taco_tensor_DC(C, nrow, ncol, {0,1});
    compute(C,&A,&B);
    if (print) {
      print_taco_tensor_DC(&A);
      print_taco_tensor_DC(&B);
      print_taco_tensor_DC(C);
    }
}