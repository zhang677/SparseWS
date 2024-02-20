#include "../../utils/dataloader.h"
#include <math.h> 

#ifndef CAP
 #define CAP 4
#endif

const int32_t INTMAX = 2147483647;

typedef struct {
  int32_t crd[2];
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

// int refresh_times = 0;

void print_hashTable(HashTable* table) {
  for (int i=0; i<table->table_size; i++) {
      if (table->values_size[i] > 0) {
      std::cout<<"indices["<<i<<"] = "<<table->values_size[i]<<", values["<<i<<"] = ";
      for (int j=0; j<table->values_size[i]; j++) {
        std::cout<<"("<<table->values[i][j].crd[0]<<","<<table->values[i][j].crd[1]<<","<<table->values[i][j].val<<"),";
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
  //std::cout << "Acc size: " << accumulator_size << ", COO size: " << COO_size << std::endl;
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
  while(accumulator_pointer < accumulator_size) {
    tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
    tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
    tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
    accumulator_pointer ++;
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
  } else {
    qsort(array, size, sizeof(wspace), w_cmp);
  }
  return size;
}


// Hash function transforms a crd to a hash key
int32_t hash_func(const int32_t crd, const int32_t tsize) {
    return crd < 0 ? (abs(INTMAX+crd) % tsize) : (crd % tsize);
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

int32_t TryInsert_hash(int32_t* insertFail, HashTable* accumulator, int32_t i, int32_t k, float val, int32_t dimension) {
  if (accumulator->numel == accumulator->buffer_capacity) {
    copy_buffer(accumulator);
    return Sort(accumulator->buffer, accumulator->buffer_capacity, false);
  }
  int32_t hashkey = hash_func(i * dimension + k, accumulator->table_size);
  //std::cout << "Hashkey: " << hashkey << std::endl;
  wspace tmp;
  tmp.crd[0] = i;
  tmp.crd[1] = k;
  tmp.val = val;
  *insertFail = 0;

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
  w->buffer_capacity = w_accumulator_capacity * 2;
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
  } else {
    w->buffer_capacity *= 1.5;
  }
  w->buffer = (wspace*)malloc(sizeof(wspace) * w->buffer_capacity);
  // refresh_times ++;
  // outfile1.close();
  // outfile2.close();
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

  HashTable w_accumulator_h[CAP]; 
  // HashTable w_accumulator;
  for (int i = 0; i < CAP; i++) {
    init_hashTable(w_accumulator_h + i, w_accumulator_capacity);
  }

//   int32_t* restrict w_point = 0;
//   w_point = (int32_t*)malloc(sizeof(int32_t) * 2);
  int32_t w_insertFail_h[CAP] = {0};
  int32_t flag_sets = CAP;
  omp_set_max_active_levels(2);
  // printf("TACO threads = {%d}\n", omp_get_max_threads());
  // printf("A1_dimension = {%d}\n", A1_dimension);

  #pragma omp parallel sections num_threads( 2 )
  {
  #pragma omp section
  #pragma omp parallel num_threads( CAP )
  {
    int thread_id = omp_get_thread_num();
    // #pragma omp for schedule(runtime)
    for (int32_t i = thread_id; i < A1_dimension; i+=CAP) {
      int32_t iA = i;
      for (int32_t jA = A2_pos[i]; jA < A2_pos[(i + 1)]; jA++) {
        int32_t j = A2_crd[jA];
        int32_t jB = j;
        for (int32_t kB = B2_pos[j]; kB < B2_pos[(j + 1)]; kB++) {
          int32_t k = B2_crd[kB];
          // printf("Thread %d is inserting\n", thread_id);
          TryInsert_hash(w_insertFail_h + thread_id, w_accumulator_h + thread_id, i, k, (A_vals[jA] * B_vals[kB]), B2_dimension);
          if (w_accumulator_h[thread_id].numel == w_accumulator_h[thread_id].buffer_capacity){
            // printf("Thread %d insertion fails\n", omp_get_thread_num());
            #pragma omp atomic write
            w_insertFail_h[thread_id] = 1;
            while (w_insertFail_h[thread_id] == 1) {;}
            // printf("Thread %d is refreshing\n", thread_id);
            refresh_wspace(w_accumulator_h + thread_id);
            TryInsert_hash(w_insertFail_h + thread_id, w_accumulator_h + thread_id, i, k, (A_vals[jA] * B_vals[kB]), B2_dimension);
          }
        }
      }
    }
    // printf("Thread %d insertion done\n", thread_id);
    if (w_accumulator_h[thread_id].numel > 0) {
      copy_buffer(w_accumulator_h + thread_id);
      Sort(w_accumulator_h[thread_id].buffer, w_accumulator_h[thread_id].numel, false);
      // printf("Thread %d finishes sorting\n", thread_id);
      #pragma omp atomic write
      w_insertFail_h[thread_id] = 2;
    } else {
      #pragma omp atomic write
      w_insertFail_h[thread_id] = 3;
    }
    while(!(w_insertFail_h[thread_id] == 4)) {;}
    // printf("Thread %d finishes final merge\n", thread_id);
    for (int l = 0; l < w_accumulator_h[thread_id].table_size; l++) {
      if (w_accumulator_h[thread_id].values_capacity[l] >0) {
        free(w_accumulator_h[thread_id].values[l]);
      }
    }
    free(w_accumulator_h[thread_id].values_size);
    free(w_accumulator_h[thread_id].values_capacity);
    free(w_accumulator_h[thread_id].buffer);
    free(w_accumulator_h[thread_id].values);
    // printf("Thread %d finishes freeing\n", thread_id);
  }

  #pragma omp section
  #pragma omp parallel num_threads( 1 )
  { 
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
    // printf("w_all_capacity: %d\n", w_all_capacity);
    
    while(flag_sets > 0) {
      for (int i = 0; i < CAP; i++) {
        if (w_insertFail_h[i] == 1) {
          // printf("Thread %d merges Thread %d\n", omp_get_thread_num(), i);
          if (w_accumulator_h[i].numel + w_all_size > w_all_capacity) {
              w_all_capacity = w_accumulator_h[i].numel + w_all_size;
              // printf("w_all_capacity: %d\n", w_all_capacity);
              w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
              w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
              w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
          }
          // printf("Pre w_all_size: %d\n", w_all_size);
          w_all_size = Merge_hash(w1_crd, w2_crd, w_vals, w_all_size, w_accumulator_h + i);
          // printf("Post w_all_size: %d\n", w_all_size);
          #pragma omp atomic write
          w_insertFail_h[i] = 0;
          // printf("Thread %d merge done\n", i);
        } else if (w_insertFail_h[i] == 2) {
          // printf("Thread %d merges Thread %d\n", omp_get_thread_num(), i);
          if (w_accumulator_h[i].numel + w_all_size > w_all_capacity) {
              w_all_capacity = w_accumulator_h[i].numel + w_all_size;
              // printf("w_all_capacity: %d\n", w_all_capacity);
              w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
              w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
              w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
          }
          // printf("Pre w_all_size: %d\n", w_all_size);
          w_all_size = Merge_hash(w1_crd, w2_crd, w_vals, w_all_size, w_accumulator_h + i);
          flag_sets --;
          // printf("Post w_all_size: %d\n", w_all_size);
          #pragma omp atomic write
          w_insertFail_h[i] = 4;
        } else if (w_insertFail_h[i] == 3) {
          flag_sets --;
          #pragma omp atomic write
          w_insertFail_h[i] = 4;
        }
      }
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
      C2_pos[k + 1] = w1_segend - kw;
      kw = w1_segend;
    }

    free(w1_crd);
    free(w1_pos);

    int32_t csC2 = 0;
    for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
      csC2 += C2_pos[pC2];
      C2_pos[pC2] = csC2;
    }

    C->indices[1][0] = (int32_t*)(C2_pos);
    C->indices[1][1] = (int32_t*)(w2_crd);
    C->vals = (float*)w_vals;
  }
  }
  return 0;
}

double CSR_CSR_hash(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t w_cap, int32_t warmup, int32_t repeat, bool bench = false, bool print = false) {
  // std::cout << "Capacity: " << w_cap << std::endl;
  for (int i = 0; i < warmup; i++) {
    compute(C,A,B,w_cap);
    if (bench) {
      free(C->vals);
      free(C->indices[1][0]);
      free(C->indices[1][1]);
    }
  }
  double start = clock();
  for (int i = 0; i < repeat; i++) {
    compute(C,A,B,w_cap);
    if (bench && i != repeat - 1) {
      free(C->vals);
      free(C->indices[1][0]);
      free(C->indices[1][1]);
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