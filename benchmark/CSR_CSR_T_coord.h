#include "../utils/dataloader.h"
#include "../utils/lib.h"
#include "../utils/wspace.h"

int32_t TryInsert_coord(bool* insertFail, wspace* accumulator, int32_t accumulator_size, int32_t accumulator_capacity, int32_t* crds, float val) {
  if (accumulator_size == accumulator_capacity) {
    *insertFail = true;
    return Sort(accumulator, accumulator_capacity, false);
  } else {
    accumulator[accumulator_size].crd[0] = crds[0];
    accumulator[accumulator_size].crd[1] = crds[1];
    accumulator[accumulator_size].val = val;
    *insertFail = false;
    accumulator_size ++;
    return accumulator_size;
  }
}

int Merge_coord(int32_t* COO1_crd, int32_t* COO2_crd, float* COO_vals, int32_t COO_size, wspace* accumulator, int32_t accumulator_size){
  if (COO_size == 0) {
      for (int i=0; i<accumulator_size; i++) {
        COO1_crd[i] = accumulator[i].crd[0];
        COO2_crd[i] = accumulator[i].crd[1];
        COO_vals[i] = accumulator[i].val;
      }
      return accumulator_size;
    }
    int32_t* tmp_COO_crd[2];
    float* tmp_COO_vals;
    for (int i=0; i<w_order; i++) {
      tmp_COO_crd[i] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
    }
    tmp_COO_vals = (float*)malloc(sizeof(float) * (accumulator_size + COO_size));
    int accumulator_pointer = 0;
    int content_pointer = 0;
    int target_pointer = 0;
    wspace tmp_con;
    while(accumulator_pointer < accumulator_size && content_pointer < COO_size) {
      tmp_con.crd[0] = COO1_crd[content_pointer];
      tmp_con.crd[1] = COO2_crd[content_pointer];
      if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) == 0) {
        for (int i=0; i<w_order; i++) {
          tmp_COO_crd[i][target_pointer] = accumulator[accumulator_pointer].crd[i];
        }
        tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val + COO_vals[content_pointer];
        accumulator_pointer ++;
        content_pointer ++;
        target_pointer ++;
      } else if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) < 0) {
        for (int i=0; i<w_order; i++) {
          tmp_COO_crd[i][target_pointer] = accumulator[accumulator_pointer].crd[i];
        }
        tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
        accumulator_pointer ++;
        target_pointer ++;
      } else {
        tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
        tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
        tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
        content_pointer ++;
        target_pointer ++;
      }
    }
    while(accumulator_pointer<accumulator_size) {
      for (int i=0; i<w_order; i++) {
        tmp_COO_crd[i][target_pointer] = accumulator[accumulator_pointer].crd[i];
      } 
      tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
      accumulator_pointer ++;
      target_pointer ++;
    }
    while(content_pointer<COO_size) {
      tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
      tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      content_pointer ++;
      target_pointer ++;
    }
    for (int i = 0; i < target_pointer; i++) {
      COO1_crd[i] = tmp_COO_crd[0][i];
      COO2_crd[i] = tmp_COO_crd[1][i];
      COO_vals[i] = tmp_COO_vals[i];
    }
    free(tmp_COO_crd[0]);
    free(tmp_COO_crd[1]);
    free(tmp_COO_vals);
    return target_pointer;
}

int compute_coo(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
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

  int32_t w_all_capacity = w_accumulator_capacity; // Ensure enough space for merging
  int32_t w_all_size = 0;
  int32_t* w1_crd = 0;
  int32_t* w2_crd = 0;
  float* w_vals = 0;
  w1_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  w2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  w_vals = (float*)malloc(sizeof(float) * w_all_capacity);
  bool* w_insertFail = 0;
  w_insertFail = (bool*)malloc(sizeof(bool) * 1);
  w_insertFail[0] = false;
  int32_t* w_point = 0;
  w_point = (int32_t*)malloc(sizeof(int32_t) * 2);
  for (int32_t i = 0; i < A1_dimension; i++) {
    for (int32_t jA = A2_pos[i]; jA < A2_pos[i+1]; jA++) {
      int32_t j = A2_crd[jA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[j+1]; kB++) {
        int32_t k = B2_crd[kB];
        w_point[0] = k;
        w_point[1] = i; // Transpose
        // Try to insert to the accumulator array
        w_accumulator_size = TryInsert_coord(w_insertFail,w_accumulator,w_accumulator_size,w_accumulator_capacity,w_point,A_vals[jA] * B_vals[kB]);
        if (w_insertFail[0]) {
          // Enlarge
          if(w_accumulator_size + w_all_size > w_all_capacity) {
            w_all_capacity = w_all_capacity * 2;
            w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
            w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
            w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
          }
          // Merge
          w_all_size = Merge_coord(w1_crd, w2_crd, w_vals, w_all_size, w_accumulator, w_accumulator_size);
          // Clear
          w_accumulator_size = 0;
          // Insert the wspace that conflicts
          w_accumulator_size = TryInsert_coord(w_insertFail, w_accumulator, w_accumulator_size, w_accumulator_capacity, w_point, A_vals[jA] * B_vals[kB]);
        }
      }
    }
  }
  if (w_accumulator_size > 0) {
    // Sort
    w_accumulator_size = Sort(w_accumulator, w_accumulator_size, false);
    // Enlarge
    if(w_accumulator_size + w_all_size > w_all_capacity) {
        w_all_capacity = w_all_capacity * 2;
        w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
        w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
        w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
    }
    // Merge
    w_all_size = Merge_coord(w1_crd, w2_crd, w_vals, w_all_size, w_accumulator, w_accumulator_size);
    // Clear
    w_accumulator_size = 0;
  }
  int w1_pos[2] = {0,w_all_size};
  pack_C(C, w1_pos, w1_crd, w2_crd, w_vals);
  return 0;
}

void CSR_CSR_T_coord(const string& A_name, const string& B_name, taco_tensor_t* C, bool print = false) {
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
    compute_coo(C,&A,&B);
    if (print) {
      print_taco_tensor_DC(&A);
      print_taco_tensor_DC(&B);
      print_taco_tensor_DC(C);
    }
}