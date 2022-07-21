#include "../utils/dataloader.h"
#include "../utils/lib.h"

/*
 * The `pack` functions convert coordinate and value arrays in COO format,
 * with nonzeros sorted lexicographically by their coordinates, to the
 * specified input format.
 *
 * The `unpack` function converts the specified output format to coordinate
 * and value arrays in COO format.
 *
 * For both, the `_COO_pos` arrays contain two elements, where the first is 0
 * and the second is the number of nonzeros in the tensor.
 */

int pack_C(taco_tensor_t *C, int* C_COO1_pos, int32_t* C_COO1_crd, int32_t* C_COO2_crd, float* C_COO_vals) {
  int C1_dimension = (int)(C->dimensions[0]);
  int* restrict C2_pos = (int*)(C->indices[1][0]);
  int* restrict C2_crd = (int*)(C->indices[1][1]);
  float* restrict C_vals = (float*)(C->vals);

  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    C2_pos[pC2] = 0;
  }
  int32_t C2_crd_size = 1048576;
  C2_crd = (int32_t*)malloc(sizeof(int32_t) * C2_crd_size);
  int32_t jC = 0;
  int32_t C_capacity = 1048576;
  C_vals = (float*)malloc(sizeof(float) * C_capacity);

  int32_t iC_COO = C_COO1_pos[0];
  int32_t pC_COO1_end = C_COO1_pos[1];

  while (iC_COO < pC_COO1_end) {
    int32_t i = C_COO1_crd[iC_COO];
    int32_t C_COO1_segend = iC_COO + 1;
    while (C_COO1_segend < pC_COO1_end && C_COO1_crd[C_COO1_segend] == i) {
      C_COO1_segend++;
    }
    int32_t pC2_begin = jC;

    int32_t jC_COO = iC_COO;

    while (jC_COO < C_COO1_segend) {
      int32_t j = C_COO2_crd[jC_COO];
      float C_COO_val = C_COO_vals[jC_COO];
      jC_COO++;
      while (jC_COO < C_COO1_segend && C_COO2_crd[jC_COO] == j) {
        C_COO_val += C_COO_vals[jC_COO];
        jC_COO++;
      }
      if (C_capacity <= jC) {
        C_vals = (float*)realloc(C_vals, sizeof(float) * (C_capacity * 2));
        C_capacity *= 2;
      }
      C_vals[jC] = C_COO_val;
      if (C2_crd_size <= jC) {
        C2_crd = (int32_t*)realloc(C2_crd, sizeof(int32_t) * (C2_crd_size * 2));
        C2_crd_size *= 2;
      }
      C2_crd[jC] = j;
      jC++;
    }

    C2_pos[i + 1] = jC - pC2_begin;
    iC_COO = C_COO1_segend;
  }

  int32_t csC2 = 0;
  for (int32_t pC20 = 1; pC20 < (C1_dimension + 1); pC20++) {
    csC2 += C2_pos[pC20];
    C2_pos[pC20] = csC2;
  }

  C->indices[1][0] = (int32_t*)(C2_pos);
  C->indices[1][1] = (int32_t*)(C2_crd);
  C->vals = (float*)C_vals;
  return 0;
}

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
    for (int32_t jA = A2_pos[i]; jA < A2_pos[i+1]; jA++) {
        int32_t j = A2_crd[jA];
        for (int32_t kB = B2_pos[j]; kB < B2_pos[j+1]; kB++) {
            int32_t k = B2_crd[kB];
            w_accumulator[w_accumulator_size].pos[0] = i;
            w_accumulator[w_accumulator_size].pos[1] = k;
            w_accumulator[w_accumulator_size].val = A_vals[jA] * B_vals[kB];
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
  print_array<int32_t>(C_COO1_pos,2);
  print_array<int32_t>(C_COO_crd[0],w_all_size);
  print_array<int32_t>(C_COO_crd[1],w_all_size);
  print_array<float>(C_COO_vals,w_all_size);
  pack_C(C, C_COO1_pos, C_COO_crd[0], C_COO_crd[1], C_COO_vals);
  return 0;
}

void CSR_CSR_4(const string& A_name, const string& B_name) {
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
    compute(&C,&A,&B);
    print_taco_tensor_DC(&A);
    print_taco_tensor_DC(&B);
    print_taco_tensor_DC(&C);
}