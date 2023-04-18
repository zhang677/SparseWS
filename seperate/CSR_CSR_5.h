#include "../utils/dataloader.h"
#include "../utils/lib.h"
#include "../utils/wspace.h"

// Hash function transforms a crd to a hash key
int hash_func(const int crd, const int tsize) {
    return (crd * HASH_SCAL) % tsize;
}

int32_t TryInsert(bool* insertFail, wspace* accumulator, int32_t accumulator_size, int32_t accumulator_capacity, int32_t* crds, float val) {
  int hashkey = hash_func(crds[0] + crds[1], accumulator_capacity);
  int step = 0;
  int bmap = 0; // The length depends on hash table size
  wspace tmp;
  tmp.crd[0] = crds[0];
  tmp.crd[1] = crds[1];
  tmp.val = val;
  do{
    if (accumulator[hashkey].crd[0] == -1){ // no conflict
      accumulator[hashkey] = tmp;
      accumulator_size ++;
      *insertFail = false;
      return accumulator_size;
    } else if (esc_cmp(&accumulator[hashkey],&tmp) == 0) {
      accumulator[hashkey].val += tmp.val;
      *insertFail = false;
      return accumulator_size;
    } else {
      bmap |= 1<<hashkey;
      hashkey = (hashkey + accumulator_capacity - 1) % accumulator_capacity; // linear-probe
    }
  }while(bmap & 1<<hashkey == 0); // Hash table still has space to insert
  *insertFail = true;
  return Sort(accumulator, accumulator_capacity, true);
}

int compute_coo_rev(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B, int32_t w_cap) {
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

  int32_t w_accumulator_capacity = w_cap; // Arbitrary
  int32_t w_accumulator_size = 0;
  wspace* restrict w_accumulator = 0;
  w_accumulator = (wspace*)malloc(sizeof(wspace) * w_accumulator_capacity);

  int32_t C_COO_capacity = w_accumulator_capacity; // Ensure enough space for merging
  int32_t* C_COO_crd[2];
  float* C_COO_vals;
  int32_t C_COO_size = 0;
  C_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * C_COO_capacity);
  C_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * C_COO_capacity);
  C_COO_vals = (float*)malloc(sizeof(float) * C_COO_capacity);
  bool insertFail;
  // Initialize hash key
  init_wspace(w_accumulator,0,w_accumulator_capacity);

  for (int32_t i = 0; i < A1_dimension; i++) {
    for (int32_t jA = A2_pos[i]; jA < A2_pos[i+1]; jA++) {
      int32_t j = A2_crd[jA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[j+1]; kB++) {
        int32_t k = B2_crd[kB];
        int32_t C_crds[2] = {k,i};
        // Try to insert to the hash array
        w_accumulator_size = TryInsert(&insertFail, w_accumulator, w_accumulator_size, w_accumulator_capacity, C_crds, A_vals[jA] * B_vals[kB]);
        //print_wspace(w_accumulator,w_accumulator_capacity);
        // No space to insert. Do contraction
        if (insertFail) {
          // Enlarge output and input arraies (optimization needed)
          if(w_accumulator_size + C_COO_size > C_COO_capacity) {
            C_COO_capacity = Enlarge(C_COO_crd,&C_COO_vals,C_COO_capacity);
          }
          // Merge
          C_COO_size = Merge(C_COO_crd, C_COO_vals, C_COO_size, w_accumulator, w_accumulator_size, true);
          //print_coo(C_COO_crd, C_COO_vals, C_COO_size);
          // Clear
          w_accumulator_size = 0;
          // Clear the hash array
          init_wspace(w_accumulator,0,w_accumulator_capacity);
          // Insert the wspace that conflicts
          w_accumulator_size = TryInsert(&insertFail, w_accumulator, w_accumulator_size, w_accumulator_capacity, C_crds, A_vals[jA] * B_vals[kB]);
        }
      }
    }
  }
  if (w_accumulator_size > 0) {
    // Sort
    w_accumulator_size = Sort(w_accumulator, w_accumulator_capacity, true);
    // Enlarge
    if(w_accumulator_size + C_COO_size > C_COO_capacity) {
      C_COO_capacity = Enlarge(C_COO_crd,&C_COO_vals,C_COO_capacity);
    }
    // Merge
    C_COO_size = Merge(C_COO_crd, C_COO_vals, C_COO_size, w_accumulator, w_accumulator_size, true);
    // Clear
    w_accumulator_size = 0;
  }
  int C_COO1_pos[2] = {0,C_COO_size};
  pack_C(C, C_COO1_pos, C_COO_crd[0], C_COO_crd[1], C_COO_vals);
  return 0;
}


void CSR_CSR_5(const string& A_name, const string& B_name, taco_tensor_t* C, int32_t w_cap,bool print = false) {
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
    init_taco_tensor_DC(C, nrow, ncol, {0,1});
    compute_coo_rev(C,&A,&B,w_cap);
    if (print) {
        print_taco_tensor_DC(&A);
        print_taco_tensor_DC(&B);
        print_taco_tensor_DC(C);
    }
}