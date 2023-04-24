#include "../utils/dataloader.h"

typedef struct {
  int32_t crd[2];
  float val;
}wspace ;
int w_cmp(const void *b, const void *a) {
  for (int i = 0; i < 2; i++) {
    if (((wspace*)b)->crd[i] == ((wspace*)a)->crd[i]) continue;
    return (((wspace*)b)->crd[i] - ((wspace*)a)->crd[i]);
  }
  return (((wspace*)b)->crd[1] - ((wspace*)a)->crd[1]);
}
int w_cmp_rev(const void *a, const void *b) {
  for (int i = 0; i < 2; i++) {
    if (((wspace*)b)->crd[i] == ((wspace*)a)->crd[i]) continue;
    return (((wspace*)b)->crd[i] - ((wspace*)a)->crd[i]);
  }
  return (((wspace*)b)->crd[1] - ((wspace*)a)->crd[1]);
}
int Merge_hash(int32_t* COO1_crd, int32_t* COO2_crd, float* COO_vals, int32_t COO_size, wspace* accumulator, int32_t accumulator_size) {
  if (COO_size == 0) {
    for (int i=0; i<accumulator_size; i++) {
      COO1_crd[i] = accumulator[accumulator_size-1-i].crd[0];
      COO2_crd[i] = accumulator[accumulator_size-1-i].crd[1];
      COO_vals[i] = accumulator[accumulator_size-1-i].val;
    }
    return accumulator_size;
  }
  int32_t* tmp_COO_crd[2];
  float* tmp_COO_vals;
  tmp_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_vals = (float*)malloc(sizeof(float) * (accumulator_size + COO_size));
  int accumulator_pointer = accumulator_size - 1;
  int content_pointer = 0;
  int target_pointer = 0;
  wspace tmp_con;
  while(accumulator_pointer >= 0 && content_pointer < COO_size) {
    tmp_con.crd[0] = COO1_crd[content_pointer];
    tmp_con.crd[1] = COO2_crd[content_pointer];
    if (w_cmp(&accumulator[accumulator_pointer], &tmp_con) == 0) {
      tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
      tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
      tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val + COO_vals[content_pointer];
      accumulator_pointer --;
      content_pointer ++;
      target_pointer ++;
    } else if (w_cmp(&accumulator[accumulator_pointer], &tmp_con) < 0) {
      tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
      tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
      tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
      accumulator_pointer --;
      target_pointer ++;
    } else {
      tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
      tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      content_pointer ++;
      target_pointer ++;
    }
  }
  while(accumulator_pointer >= 0) {
    tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
    tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
    tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
    accumulator_pointer --;
    target_pointer ++;
  }
  while(content_pointer < COO_size) {
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
int Sort(void* array, size_t size, bool rev) {
  if (rev) {
    qsort(array, size, sizeof(wspace), w_cmp_rev);
    for (int i = size - 1; i >= 0; i--) {
      if(((wspace*)array)[i].crd[0] != -1) {
        return i + 1;
      }
    }
  } 
  qsort(array, size, sizeof(wspace), w_cmp);
  return size;
}
const int HASH_SCAL = 1;
// Hash function transforms a crd to a hash key
int hash_func(const int crd, const int tsize) {
    return (crd * HASH_SCAL) % tsize;
}

int32_t TryInsert_hash(bool* insertFail, wspace* accumulator, int32_t accumulator_size, int32_t accumulator_capacity, int32_t* crds, float val, int32_t* bitmap) {
  int hashkey = hash_func(crds[0] + crds[1], accumulator_capacity);
  int step = 0;
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
    } else if (w_cmp(&accumulator[hashkey],&tmp) == 0) {
      accumulator[hashkey].val += tmp.val;
      *insertFail = false;
      return accumulator_size;
    } else {
      bitmap[hashkey / 32] |= 1<<(hashkey % 32);
      hashkey = (hashkey + accumulator_capacity - 1) % accumulator_capacity; // linear-probe
    }
  }while(bitmap[hashkey / 32] & 1<<(hashkey % 32) == 0); // Hash table still has space to insert
  *insertFail = true;
  return Sort(accumulator, accumulator_capacity, true);
}

bool init_wspace(wspace* w, int begin, int end) {
  for(int i = begin; i < end; i++){
    w[i].crd[0] = -1;
  }
  // memset(w, -1, sizeof(wspace) * (end - begin));
  return true;
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

  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    C2_pos[pC2] = 0;
  }
  int32_t C2_crd_size = 1048576;
  C2_crd = (int32_t*)malloc(sizeof(int32_t) * C2_crd_size);
  int32_t iC = 0;
  int32_t C_capacity = 1048576;
  C_vals = (float*)malloc(sizeof(float) * C_capacity);

  int32_t w_accumulator_capacity = 3;
  int32_t w_accumulator_size = 0;
  wspace* restrict w_accumulator = 0;
  w_accumulator = (wspace*)malloc(sizeof(wspace) * w_accumulator_capacity);
  int32_t w_all_capacity = w_accumulator_capacity;
  int32_t w_all_size = 0;
  int32_t* restrict w1_pos = 0;
  w1_pos = (int32_t*)malloc(sizeof(int32_t) * 2);
  int32_t* restrict w1_crd = 0;
  int32_t* restrict w2_crd = 0;
  w1_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  w2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  float* restrict w_vals = 0;
  w_vals = (float*)malloc(sizeof(float) * w_all_capacity);
  bool* restrict w_insertFail = 0;
  w_insertFail = (bool*)malloc(sizeof(bool) * 1);
  w_insertFail[0] = 0;
  int32_t* restrict w_point = 0;
  w_point = (int32_t*)malloc(sizeof(int32_t) * 2);
  init_wspace(w_accumulator,0,w_accumulator_capacity);
  int32_t* bitmap = (int32_t*)malloc(sizeof(int32_t) * (w_accumulator_capacity / 32 + 1));


  for (int32_t i = 0; i < A1_dimension; i++) {
    w_point[1] = i;
    int32_t iA = i;
    for (int32_t jA = A2_pos[i]; jA < A2_pos[(i + 1)]; jA++) {
      int32_t j = A2_crd[jA];
      int32_t jB = j;
      for (int32_t kB = B2_pos[j]; kB < B2_pos[(j + 1)]; kB++) {
        int32_t k = B2_crd[kB];
        w_point[0] = k;
        w_accumulator_size = TryInsert_hash(w_insertFail, w_accumulator, w_accumulator_size, w_accumulator_capacity, w_point, (A_vals[jA] * B_vals[kB]), bitmap);
        if (w_insertFail[0]) {
          if (w_accumulator_size + w_all_size > w_all_capacity) {
            w_all_capacity *= 2;
            w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
            w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
            w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
          }
          w_all_size = Merge_hash(w1_crd, w2_crd, w_vals, w_all_size, w_accumulator, w_accumulator_size);
          w_accumulator_size = 0;
          init_wspace(w_accumulator,0,w_accumulator_capacity);
          w_accumulator_size = TryInsert_hash(w_insertFail, w_accumulator, w_accumulator_size, w_accumulator_capacity, w_point, (A_vals[jA] * B_vals[kB]), bitmap);
        }
      }
    }
  }
  if (w_accumulator_size > 0) {
    w_accumulator_size = Sort(w_accumulator, w_accumulator_size, 1);
    if (w_accumulator_size + w_all_size > w_all_capacity) {
      w_all_capacity *= 2;
      w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
      w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
      w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
    }
    w_all_size = Merge_hash(w1_crd, w2_crd, w_vals, w_all_size, w_accumulator, w_accumulator_size);
    w_accumulator_size = 0;
  }
  w1_pos[0] = 0;
  w1_pos[1] = w_all_size;
  int32_t kw = w1_pos[0];
  int32_t pw1_end = w1_pos[1];

  while (kw < pw1_end) {
    int32_t k = w1_crd[kw];
    int32_t w1_segend = kw + 1;
    while (w1_segend < pw1_end && w1_crd[w1_segend] == k) {
      w1_segend++;
    }
    int32_t kC = k;
    int32_t pC2_begin = iC;

    for (int32_t iw = kw; iw < w1_segend; iw++) {
      int32_t i = w2_crd[iw];
      if (C_capacity <= iC) {
        C_vals = (float*)realloc(C_vals, sizeof(float) * (C_capacity * 2));
        C_capacity *= 2;
      }
      C_vals[iC] = w_vals[iw];
      if (C2_crd_size <= iC) {
        C2_crd = (int32_t*)realloc(C2_crd, sizeof(int32_t) * (C2_crd_size * 2));
        C2_crd_size *= 2;
      }
      C2_crd[iC] = i;
      iC++;
    }

    C2_pos[k + 1] = iC - pC2_begin;
    kw = w1_segend;
  }
  free(w_accumulator);
  free(w_vals);
  free(w_insertFail);
  free(w_point);
  free(w1_crd);
  free(w2_crd);
  free(w1_pos);
  free(bitmap);

  int32_t csC2 = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    csC2 += C2_pos[pC2];
    C2_pos[pC2] = csC2;
  }

  C->indices[1][0] = (int32_t*)(C2_pos);
  C->indices[1][1] = (int32_t*)(C2_crd);
  C->vals = (float*)C_vals;
  return 0;
}

void CSR_CSR_T_hash_ext(const string& A_name, const string& B_name, taco_tensor_t* C, bool print = false) {
  // C(k,i) = A(i,j) * B(j,k); C: CSR, A: CSR, B: CSR
  // Hash needs to be initialized (Different from Coord)
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
  compute(C,&A,&B);
  if (print) {
    print_taco_tensor_DC(&A);
    print_taco_tensor_DC(&B);
    print_taco_tensor_DC(C);
  }
}