#include "../../utils/dataloader.h"
#include <math.h> 

const int32_t INTMAX = 2147483647;

typedef struct {
  int32_t crd[3];
  float val;
}wspace;

typedef struct {
  int32_t* values_size;
  int32_t* values_capacity;
  wspace** values;
  int32_t numel;
  int32_t table_size;
  int32_t buffer_capacity;
  wspace* buffer;
}HashTable;

int refresh_times = 0;

void print_hashTable(HashTable* table) {
  for (int i=0; i<table->table_size; i++) {
    if (table->values_size[i] > 0) {
      std::cout<<"indices["<<i<<"] = "<<table->values_size[i]<<", values["<<i<<"] = ";
      for (int j=0; j<table->values_size[i]; j++) {
        std::cout<<"("<<table->values[i][j].crd[0]<<","<<table->values[i][j].crd[1]<<","<<table->values[i][j].crd[1]<<","<<table->values[i][j].val<<"),";
      }
      std::cout<<std::endl;
    }
  }
  std::cout<<std::endl;
  // std::cout<<"Buffer: ";
  // for (int i = 0; i < table->numel; i++) {
  //   std::cout<<"("<<table->buffer[i].crd[0]<<","<<table->buffer[i].crd[1]<<","<<table->buffer[i].val<<"),";
  // }
  // std::cout<<std::endl;
}

int w_cmp(const void *b, const void *a) {
  for (int i = 0; i < 3; i++) {
    if (((wspace*)b)->crd[i] == ((wspace*)a)->crd[i]) continue;
    return (((wspace*)b)->crd[i] - ((wspace*)a)->crd[i]);
  }
  return (((wspace*)b)->crd[2] - ((wspace*)a)->crd[2]);
}
int w_cmp_rev(const void *a, const void *b) {
  for (int i = 0; i < 3; i++) {
    if (((wspace*)b)->crd[i] == ((wspace*)a)->crd[i]) continue;
    return (((wspace*)b)->crd[i] - ((wspace*)a)->crd[i]);
  }
  return (((wspace*)b)->crd[2] - ((wspace*)a)->crd[2]);
}
int Merge_hash(int32_t* COO1_crd, int32_t* COO2_crd, int32_t* COO3_crd, float* COO_vals, int32_t COO_size, HashTable* w) {
  wspace* accumulator = w->buffer;
  int accumulator_size = w->numel;
  // std::cout << "Acc size: " << accumulator_size << ", COO size: " << COO_size << std::endl;
  if (COO_size == 0) {
    for (int i=0; i<accumulator_size; i++) {
      COO1_crd[i] = accumulator[i].crd[0];
      COO2_crd[i] = accumulator[i].crd[1];
      COO3_crd[i] = accumulator[i].crd[2];
      COO_vals[i] = accumulator[i].val;
    }
    return accumulator_size;
  }
  int32_t* tmp_COO_crd[3];
  float* tmp_COO_vals;
  tmp_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_crd[2] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
  tmp_COO_vals = (float*)malloc(sizeof(float) * (accumulator_size + COO_size));
  int accumulator_pointer = 0;
  int content_pointer = 0;
  int target_pointer = 0;
  wspace tmp_con;
  while(accumulator_pointer < accumulator_size && content_pointer < COO_size) {
    tmp_con.crd[0] = COO1_crd[content_pointer];
    tmp_con.crd[1] = COO2_crd[content_pointer];
    tmp_con.crd[2] = COO3_crd[content_pointer];
    if (w_cmp(&accumulator[accumulator_pointer], &tmp_con) == 0) {
      tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
      tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
      tmp_COO_crd[2][target_pointer] = accumulator[accumulator_pointer].crd[2];
      tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val + COO_vals[content_pointer];
      accumulator_pointer ++;
      content_pointer ++;
      target_pointer ++;
    } else if (w_cmp(&accumulator[accumulator_pointer], &tmp_con) < 0) {
      tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
      tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
      tmp_COO_crd[2][target_pointer] = accumulator[accumulator_pointer].crd[2];
      tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
      accumulator_pointer ++;
      target_pointer ++;
    } else {
      tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
      tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
      tmp_COO_crd[2][target_pointer] = COO3_crd[content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      content_pointer ++;
      target_pointer ++;
    }
  }
  while(accumulator_pointer < accumulator_size) {
    tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
    tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
    tmp_COO_crd[2][target_pointer] = accumulator[accumulator_pointer].crd[2];
    tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
    accumulator_pointer ++;
    target_pointer ++;
  }
  while(content_pointer < COO_size) {
    tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
    tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
    tmp_COO_crd[2][target_pointer] = COO3_crd[content_pointer];
    tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
    content_pointer ++;
    target_pointer ++;
  }
  for (int i = 0; i < target_pointer; i++) {
    COO1_crd[i] = tmp_COO_crd[0][i];
    COO2_crd[i] = tmp_COO_crd[1][i];
    COO3_crd[i] = tmp_COO_crd[2][i];
    COO_vals[i] = tmp_COO_vals[i];
  }
  free(tmp_COO_crd[0]);
  free(tmp_COO_crd[1]);
  free(tmp_COO_crd[2]);
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


// Hash function transforms a crd to a hash key
int32_t hash_func(const int32_t crd, const int32_t tsize) {
    return crd < 0 ? ((INTMAX+crd) % tsize) : (crd % tsize);
}

void copy_buffer(HashTable* accumulator) {
  int current_size = 0;
  int32_t* values_size = accumulator->values_size;
  wspace** values = accumulator->values;
  for (int i = 0; i < accumulator->table_size; i++) {
    if (values_size[i] != 0) {
      memcpy(accumulator->buffer + current_size, values[i], sizeof(wspace) * values_size[i]);
      current_size += values_size[i];
    }
  }
}

int32_t TryInsert_hash(bool* insertFail, HashTable* accumulator, int32_t* crds, float val, int32_t* dims) {
  if (accumulator->numel == accumulator->buffer_capacity) {
    copy_buffer(accumulator);
    *insertFail = true;
    return Sort(accumulator->buffer, accumulator->buffer_capacity, false);
  }
  int32_t hashkey = hash_func(crds[0] * dims[0] * dims[1] + crds[1] * dims[1] + crds[2], accumulator->table_size);
  //std::cout << "Hashkey: " << hashkey << std::endl;
  wspace tmp;
  tmp.crd[0] = crds[0];
  tmp.crd[1] = crds[1];
  tmp.crd[2] = crds[2];
  tmp.val = val;
  *insertFail = false;

  int32_t* values_size = accumulator->values_size;
  int32_t* values_capacity = accumulator->values_capacity;
  wspace** values = accumulator->values;
  int32_t size = values_size[hashkey];
  int32_t capacity = values_capacity[hashkey];
  for (int i = 0; i < size; i++) {
    if (w_cmp(&values[hashkey][i], &tmp) == 0) {
      values[hashkey][i].val += val;
      return accumulator->numel;
    }
  }
  if (capacity == 0) {
    values[hashkey] = (wspace*)malloc(sizeof(wspace) * 2);
    values_capacity[hashkey] = 2;
  } else {
    if (size == capacity) {
      values[hashkey] = (wspace*)realloc(values[hashkey], sizeof(wspace) * capacity * 2);
      values_capacity[hashkey] = capacity * 2;
    }
  }
  values[hashkey][size] = tmp;
  values_size[hashkey] ++;
  accumulator->numel ++;
  //std::cout << "inserted " << accumulator->numel << std::endl;
  return accumulator->numel;
}

void init_hashTable(HashTable* w, int32_t w_accumulator_capacity) {
  w->numel = 0;
  w->table_size = w_accumulator_capacity;
  w->buffer_capacity = max(INT32_MAX, w_accumulator_capacity * 2);
  w->values = (wspace**)malloc(sizeof(wspace*) * w->table_size);
  w->values_size = (int32_t*)calloc(w->table_size, sizeof(int32_t));
  w->values_capacity = (int32_t*)calloc(w->table_size, sizeof(int32_t));
  w->buffer = (wspace*)malloc(sizeof(wspace) * w->buffer_capacity);
}

void refresh_wspace(HashTable* w) {
  // Save the size and capacity of each bucket
  // std::string filename1 = "/home/nfs_data/zhanggh/SparseWS/sizes" + std::to_string(refresh_times) + ".txt";
  // std::string filename2 = "/home/nfs_data/zhanggh/SparseWS/capacity" + std::to_string(refresh_times) + ".txt";
  // std::ofstream outfile1(filename1, std::ios_base::out);
  // std::ofstream outfile2(filename2, std::ios_base::out);

  for(int i = 0; i < w->table_size; i++){
    // outfile1 << w->values_size[i] << " ";
    // outfile2 << w->values_capacity[i] << " ";
    if (w->values_capacity[i] != 0) {
      
      free(w->values[i]);
      w->values_capacity[i] = 0;
      w->values_size[i] = 0;
    }
  }

  w->numel = 0;
  free(w->buffer);
  // Simple heuristic to avoid frequent reallocation
  if (w->buffer_capacity < 4194304) {
    w->buffer_capacity *= 2;
  } else if (w->buffer_capacity < 33554432) {
    w->buffer_capacity *= 1.5;
  }
  // std::cout << "Buff size: " << w->buffer_capacity << std::endl;

  w->buffer = (wspace*)malloc(sizeof(wspace) * w->buffer_capacity);
  refresh_times ++;
  // outfile1.close();
  // outfile2.close();
}

int compute(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C, int32_t w_accumulator_capacity) {
  int A2_dimension = (int)(A->dimensions[1]);
  int A3_dimension = (int)(A->dimensions[2]);
  int* restrict B1_pos = (int*)(B->indices[0][0]);
  int* restrict B1_crd = (int*)(B->indices[0][1]);
  int* restrict B2_pos = (int*)(B->indices[1][0]);
  int* restrict B2_crd = (int*)(B->indices[1][1]);
  int* restrict B3_pos = (int*)(B->indices[2][0]);
  int* restrict B3_crd = (int*)(B->indices[2][1]);
  float* restrict B_vals = (float*)(B->vals);
  int* restrict C1_pos = (int*)(C->indices[0][0]);
  int* restrict C1_crd = (int*)(C->indices[0][1]);
  int* restrict C2_pos = (int*)(C->indices[1][0]);
  int* restrict C2_crd = (int*)(C->indices[1][1]);
  float* restrict C_vals = (float*)(C->vals);

  HashTable w_accumulator;
  init_hashTable(&w_accumulator, w_accumulator_capacity);
  int32_t w_all_capacity = w_accumulator_capacity;
  int32_t w_all_size = 0;
  int32_t* restrict w1_crd = 0;
  int32_t* restrict w2_crd = 0;
  int32_t* restrict w3_crd = 0;
  w1_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  w2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  w3_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  float* restrict w_vals = 0;
  w_vals = (float*)malloc(sizeof(float) * w_all_capacity);
  bool* restrict w_insertFail = 0;
  w_insertFail = (bool*)malloc(sizeof(bool) * 1);
  w_insertFail[0] = 0;
  int32_t* restrict w_point = 0;
  w_point = (int32_t*)malloc(sizeof(int32_t) * 3);
  int32_t restrict w_dims[2] = {A2_dimension, 4};
  // int32_t restrict w_dims[2] = {A2_dimension, A3_dimension};

  // int iterations = 0;
  // std::vector<int> w_all_caps;
  // w_all_caps.push_back(w_all_capacity);
  for (int32_t kB = B1_pos[0]; kB < B1_pos[1]; kB++) {
    int32_t k = B1_crd[kB];
    int32_t iB = B2_pos[kB];
    int32_t pB2_end = B2_pos[(kB + 1)];
    int32_t iC = C1_pos[0];
    int32_t pC1_end = C1_pos[1];
    w_point[0] = k;
    // printf("k: %d\n", k);
    while (iB < pB2_end && iC < pC1_end) {
        int32_t iB0 = B2_crd[iB];
        int32_t iC0 = C1_crd[iC];
        int32_t i = TACO_MIN(iB0,iC0);
        if (iB0 == i && iC0 == i) {
            for (int32_t jB = B3_pos[iB]; jB < B3_pos[(iB + 1)]; jB++) {
                int32_t j = B3_crd[jB];
                int32_t jA = k * A2_dimension + j;
                w_point[1] = j;
                for (int32_t lC = C2_pos[iC]; lC < C2_pos[(iC + 1)]; lC++) {
                    int32_t l = C2_crd[lC];
                    int32_t lA = jA * A3_dimension + l;
                    w_point[2] = l;
                    TryInsert_hash(w_insertFail, &w_accumulator, w_point, (B_vals[jB] * C_vals[lC]), w_dims);
                    // iterations ++;
                    if (w_insertFail[0]) {
                        if (w_accumulator.numel + w_all_size > w_all_capacity) {
                            w_all_capacity = w_accumulator.numel + w_all_size;
                            w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
                            w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
                            w3_crd = (int32_t*)realloc(w3_crd, sizeof(int32_t) * w_all_capacity);
                            w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
                            // w_all_caps.push_back(w_all_capacity);
                        }
                        
                        w_all_size = Merge_hash(w1_crd, w2_crd, w3_crd, w_vals, w_all_size, &w_accumulator);
                        // print_array(w1_crd, w_all_size);
                        // print_array(w2_crd, w_all_size);
                        // print_array(w_vals, w_all_size);
                        refresh_wspace(&w_accumulator);
                        TryInsert_hash(w_insertFail, &w_accumulator, w_point, (B_vals[jB] * C_vals[lC]), w_dims);
                    }
                }
            }
        }
        iB += (int32_t)(iB0 == i);
        iC += (int32_t)(iC0 == i);
    }
  }
  
  if (w_accumulator.numel > 0) {
    copy_buffer(&w_accumulator);
    Sort(w_accumulator.buffer, w_accumulator.numel, false);
    if (w_accumulator.numel + w_all_size > w_all_capacity) {
      w_all_capacity = w_accumulator.numel + w_all_size;
      w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
      w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
      w3_crd = (int32_t*)realloc(w3_crd, sizeof(int32_t) * w_all_capacity);
      w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
      // w_all_caps.push_back(w_all_capacity);
    }
    w_all_size = Merge_hash(w1_crd, w2_crd, w3_crd, w_vals, w_all_size, &w_accumulator);
  }
  // float mean = (float)w_accumulator.numel / (float)w_accumulator.table_size;
  // float variance = 0.0f;
  for (int i = 0; i < w_accumulator.table_size; i++) {
    if (w_accumulator.values_capacity[i] >0) {
      // variance += (w_accumulator.values_size[i] - mean) * (w_accumulator.values_size[i] - mean);
      free(w_accumulator.values[i]);
    }
  }
  // variance /= w_accumulator.table_size;
  // std::cout << "Mean: " << mean << ", Variance: " << variance << std::endl;
  free(w_accumulator.values_size);
  free(w_accumulator.values_capacity);
  free(w_accumulator.buffer);
  //free(w_accumulator.values);
  free(w_insertFail);
  free(w_point);

  A->indices[1][0] = (int32_t*)(w1_crd);
  A->indices[2][0] = (int32_t*)(w2_crd);
  A->indices[3][0] = (int32_t*)(w3_crd);
  A->vals = (float*)w_vals;
  A->vals_size = w_all_size;

  // std::cout << "Iterations: " << iterations << std::endl;
  // std::cout << "w_all_caps: ";
  // for (int i = 0; i < w_all_caps.size(); i++) {
  //   std::cout << w_all_caps[i] << ",";
  // }
  // std::cout << std::endl;
  return 0;
}

double COO_CSF_DCSR_hash(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t w_cap, int32_t warmup, int32_t repeat, bool bench = false, bool print = false) {
  for (int i = 0; i < warmup; i++) {
    compute(A,B,C,w_cap);
    if (bench) {
      free(A->vals);
      free(A->indices[1][0]);
      free(A->indices[2][0]);
      free(A->indices[3][0]);
    }
  }
  double start = clock();
  for (int i = 0; i < repeat; i++) {
    compute(A,B,C,w_cap);
    if (bench && i != repeat - 1) {
      free(A->vals);
      free(A->indices[1][0]);
      free(A->indices[2][0]);
      free(A->indices[3][0]);
    }
  }
  double end = clock();
  double duration = (double)(end - start) / (CLOCKS_PER_SEC * repeat);
  if (print) {
    print_taco_tensor_DC(A);
    print_taco_tensor_DC(B);
    print_taco_tensor_DC(C);
  }
  return duration;
}