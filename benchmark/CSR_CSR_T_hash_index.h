#include "../utils/dataloader.h"

#define w_accumulator_capacity 1048576

typedef struct {
  int32_t crd[2];
  float val;
}wspace;

typedef struct {
  int32_t* indices;
  wspace** values;
  int32_t numel;
  int32_t table_size;
  int32_t buffer_capacity;
  wspace* buffer;
}HashTable;

void print_hashTable(HashTable* table) {
  for (int i=0; i<table->table_size; i++) {
    std::cout<<"indices["<<i<<"] = "<<table->indices[i]<<", values["<<i<<"] = ";
    for (int j=0; j<table->indices[i]; j++) {
      std::cout<<"("<<table->values[i][j].crd[0]<<","<<table->values[i][j].crd[1]<<","<<table->values[i][j].val<<"),";
    }
    std::cout<<std::endl;
  }
  std::cout<<"Buffer: ";
  for (int i = 0; i < table->numel; i++) {
    std::cout<<"("<<table->buffer[i].crd[0]<<","<<table->buffer[i].crd[1]<<","<<table->buffer[i].val<<"),";
  }
  std::cout<<std::endl;
}

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
int Merge_hash(int32_t* COO1_crd, int32_t* COO2_crd, float* COO_vals, int32_t COO_size, HashTable* w) {
  wspace* accumulator = w->buffer;
  int accumulator_size = w->numel;
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
  tmp_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_vals = (float*)malloc(sizeof(float) * (accumulator_size + COO_size));
  int accumulator_pointer = 0;
  int content_pointer = 0;
  int target_pointer = 0;
  wspace tmp_con;
  while(accumulator_pointer < accumulator_size && content_pointer < COO_size) {
    tmp_con.crd[0] = COO1_crd[content_pointer];
    tmp_con.crd[1] = COO2_crd[content_pointer];
    if (w_cmp(&accumulator[accumulator_pointer], &tmp_con) == 0) {
      tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
      tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
      tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val + COO_vals[content_pointer];
      accumulator_pointer ++;
      content_pointer ++;
      target_pointer ++;
    } else if (w_cmp(&accumulator[accumulator_pointer], &tmp_con) < 0) {
      tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
      tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
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
    tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
    tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
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
int Sort(void* array, size_t size, bool rev) {
  if (rev) {
    qsort(array, size, sizeof(wspace), w_cmp_rev);
  } else {
    qsort(array, size, sizeof(wspace), w_cmp);
  }
  return size;
}
const int HASH_SCAL = 1;
// Hash function transforms a crd to a hash key
int hash_func(const int crd, const int tsize) {
    return (crd * HASH_SCAL) % tsize;
}

void copy_buffer(HashTable* accumulator) {
  int current_size = 0;
  int32_t* indices = accumulator->indices;
  wspace** values = accumulator->values;
  for (int i = 0; i < accumulator->table_size; i++) {
    if (indices[i] != 0) {
      memcpy(accumulator->buffer + current_size, values[i], sizeof(wspace) * indices[i]);
      current_size += indices[i];
    }
  }
}

int32_t TryInsert_hash(bool* insertFail, HashTable* accumulator, int32_t* crds, float val) {
  int hashkey = hash_func(crds[0] + crds[1], accumulator->table_size);
  wspace tmp;
  tmp.crd[0] = crds[0];
  tmp.crd[1] = crds[1];
  tmp.val = val;
  *insertFail = false;

  if (accumulator->indices[hashkey] == 0) {
    accumulator->indices[hashkey] = 1;
    accumulator->values[hashkey][0] = tmp;
    accumulator->numel ++;
    return accumulator->numel;
  } else {
    int32_t* indices = accumulator->indices;
    wspace** values = accumulator->values;
    int32_t size = indices[hashkey];
    for (int i = 0; i < size; i++) {
      if (w_cmp(&values[hashkey][i], &tmp) == 0) {
        values[hashkey][i].val += val;
        return accumulator->numel;
      }
    }
    if (accumulator->numel == accumulator->buffer_capacity) {
      copy_buffer(accumulator);
      *insertFail = true;
      return Sort(accumulator->buffer, accumulator->buffer_capacity, false);
    }
    values[hashkey] = (wspace*)realloc(values[hashkey], sizeof(wspace) * size * 2);
    values[hashkey][size] = tmp;
    indices[hashkey] ++;
    accumulator->numel ++;
    return accumulator->numel;
  }
}

void refresh_wspace(HashTable* w) {
  for(int i = 0; i < w->table_size; i++){
    free(w->values[i]);
    w->values[i] = (wspace*)malloc(sizeof(wspace*));
    w->indices[i] = 0;
  }
  w->numel = 0;
}

int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B, int32_t w_accumulator_capacity) {
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

  HashTable w_accumulator;
  w_accumulator.numel = 0;
  w_accumulator.table_size = w_accumulator_capacity;
  w_accumulator.buffer_capacity = w_accumulator_capacity * 2;
  w_accumulator.values = (wspace**)malloc(sizeof(wspace*) * w_accumulator.table_size);
  w_accumulator.indices = (int32_t*)calloc(sizeof(int32_t) * w_accumulator.table_size, sizeof(int32_t));
  w_accumulator.buffer = (wspace*)malloc(sizeof(wspace) * w_accumulator.buffer_capacity);

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
  refresh_wspace(&w_accumulator);

  for (int32_t i = 0; i < A1_dimension; i++) {
    w_point[1] = i;
    int32_t iA = i;
    for (int32_t jA = A2_pos[i]; jA < A2_pos[(i + 1)]; jA++) {
      int32_t j = A2_crd[jA];
      int32_t jB = j;
      for (int32_t kB = B2_pos[j]; kB < B2_pos[(j + 1)]; kB++) {
        int32_t k = B2_crd[kB];
        w_point[0] = k;
        TryInsert_hash(w_insertFail, &w_accumulator, w_point, (A_vals[jA] * B_vals[kB]));
        if (w_insertFail[0]) {
          if (w_accumulator.numel + w_all_size > w_all_capacity) {
            w_all_capacity *= 2;
            w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
            w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
            w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
          }
          w_all_size = Merge_hash(w1_crd, w2_crd, w_vals, w_all_size, &w_accumulator);
          refresh_wspace(&w_accumulator);
          TryInsert_hash(w_insertFail, &w_accumulator, w_point, (A_vals[jA] * B_vals[kB]));
        }
      }
    }
  }
  if (w_accumulator.numel > 0) {
    copy_buffer(&w_accumulator);
    Sort(w_accumulator.buffer, w_accumulator.numel, false);
    if (w_accumulator.numel + w_all_size > w_all_capacity) {
      w_all_capacity *= 2;
      w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
      w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
      w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
    }
    w_all_size = Merge_hash(w1_crd, w2_crd, w_vals, w_all_size, &w_accumulator);
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
  for (int i = 0; i < w_accumulator.table_size; i++) {
    free(w_accumulator.values[i]);
  }
  free(w_accumulator.indices);
  free(w_accumulator.buffer);
  free(w_vals);
  free(w_insertFail);
  free(w_point);
  free(w1_crd);
  free(w2_crd);
  free(w1_pos);

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

void CSR_CSR_T_hash(const string& A_name, const string& B_name, taco_tensor_t* C, int32_t w_cap, bool print = false) {
  // C(k,i) = A(i,j) * B(j,k); C: CSR, A: CSR, B: CSR
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
  compute(C,&A,&B, w_cap);
  if (print) {
    print_taco_tensor_DC(&A);
    print_taco_tensor_DC(&B);
    print_taco_tensor_DC(C);
  }
}