#include "../utils/dataloader.h"
#include "../utils/lib.h"
using namespace std;
const int w_order = 2;
const int HASH_SCAL = 1;
struct wspace {
  int32_t crd[2];
  float val;
};
// Initialize wspace
void init_wspace(wspace* w) {
    w->crd[0] = -1;
    //for (int i=0; i<w_order; i++) {
    //    w->crd[i] = -1;
    //}
}
// cmp function is generated when knowing the dimension of pos in wspace
int esc_cmp(const void *b, const void *a) {
  const int psize = w_order;
  for (int i=0; i<psize; i++) {
    if (((wspace*)b)->crd[i] == ((wspace*)a)->crd[i]) continue;
    return (((wspace*)b)->crd[i] - ((wspace*)a)->crd[i]);
  }
  return (((wspace*)b)->crd[psize-1] - ((wspace*)a)->crd[psize-1]);
}
// Hash function transforms a crd to a hash key
int hash_func(const int crd, const int tsize) {
    return (crd * HASH_SCAL) % tsize;
}
int esc_cmp_rev(const void *a, const void *b) {
  const int psize = w_order;
  for (int i=0; i<psize; i++) {
    if (((wspace*)b)->crd[i] == ((wspace*)a)->crd[i]) continue;
    return (((wspace*)b)->crd[i] - ((wspace*)a)->crd[i]);
  }
  return (((wspace*)b)->crd[psize-1] - ((wspace*)a)->crd[psize-1]);
}
void print_wspace(wspace* w, const int len) {
  for(int i=0; i<len; i++) {
    cout<<"("<<w[i].crd[0]<<","<<w[i].crd[1]<<","<<w[i].val<<"),";
  }
  cout<<endl;
}
void print_coo(int32_t** C_COO_crd, float* C_COO_vals,const int len) {
  for(int i=0; i<len; i++) {
    cout<<"("<<C_COO_crd[0][i]<<","<<C_COO_crd[1][i]<<","<<C_COO_vals[i]<<"),";
  }
  cout<<endl;
}
int Enlarge(int32_t** C_COO_crd, float** C_COO_vals, int C_COO_capacity) {
  C_COO_capacity = C_COO_capacity * 2;
  C_COO_crd[0] = (int32_t*)realloc(C_COO_crd[0], sizeof(int32_t) * C_COO_capacity);
  C_COO_crd[1] = (int32_t*)realloc(C_COO_crd[1], sizeof(int32_t) * C_COO_capacity);
  *C_COO_vals = (float*)realloc(*C_COO_vals, sizeof(float) * C_COO_capacity);
  return C_COO_capacity;
}
int merge_coo(int32_t** COO_crd, float* COO_vals, int32_t COO_size, wspace* accumulator, int32_t accumulator_size) {
  if (COO_size == 0) {
    for (int i=0; i<accumulator_size; i++) {
      COO_crd[0][i] = accumulator[i].crd[0];
      COO_crd[1][i] = accumulator[i].crd[1];
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
    tmp_con.crd[0] = COO_crd[0][content_pointer];
    tmp_con.crd[1] = COO_crd[1][content_pointer];
    if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) == 0) {
      tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
      tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
      tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val + COO_vals[content_pointer];
      accumulator_pointer ++;
      content_pointer ++;
      target_pointer ++;
    } else if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) < 0) {
      tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
      tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
      tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
      accumulator_pointer ++;
      target_pointer ++;
    } else {
      tmp_COO_crd[0][target_pointer] = COO_crd[0][content_pointer];
      tmp_COO_crd[1][target_pointer] = COO_crd[1][content_pointer];
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
    tmp_COO_crd[0][target_pointer] = COO_crd[0][content_pointer];
    tmp_COO_crd[1][target_pointer] = COO_crd[1][content_pointer];
    tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
    content_pointer ++;
    target_pointer ++;
  }
  for (int i = 0; i < target_pointer; i++) {
    COO_crd[0][i] = tmp_COO_crd[0][i];
    COO_crd[1][i] = tmp_COO_crd[1][i];
    COO_vals[i] = tmp_COO_vals[i];
  }
  free(tmp_COO_crd[0]);
  free(tmp_COO_crd[1]);
  free(tmp_COO_vals);
  return target_pointer;
}

int Merge(int32_t** COO_crd, float* COO_vals, int32_t COO_size, wspace* accumulator, int32_t accumulator_size, bool rev=false) {
  if (rev) {
    if (COO_size == 0) {
      for (int i=0; i<accumulator_size; i++) {
        COO_crd[0][i] = accumulator[accumulator_size-1-i].crd[0];
        COO_crd[1][i] = accumulator[accumulator_size-1-i].crd[1];
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
      tmp_con.crd[0] = COO_crd[0][content_pointer];
      tmp_con.crd[1] = COO_crd[1][content_pointer];
      if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) == 0) {
        tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
        tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
        tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val + COO_vals[content_pointer];
        accumulator_pointer --;
        content_pointer ++;
        target_pointer ++;
      } else if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) < 0) {
        tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
        tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
        tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
        accumulator_pointer --;
        target_pointer ++;
      } else {
        tmp_COO_crd[0][target_pointer] = COO_crd[0][content_pointer];
        tmp_COO_crd[1][target_pointer] = COO_crd[1][content_pointer];
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
      tmp_COO_crd[0][target_pointer] = COO_crd[0][content_pointer];
      tmp_COO_crd[1][target_pointer] = COO_crd[1][content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      content_pointer ++;
      target_pointer ++;
    }
    for (int i = 0; i < target_pointer; i++) {
      COO_crd[0][i] = tmp_COO_crd[0][i];
      COO_crd[1][i] = tmp_COO_crd[1][i];
      COO_vals[i] = tmp_COO_vals[i];
    }
    free(tmp_COO_crd[0]);
    free(tmp_COO_crd[1]);
    free(tmp_COO_vals);
    return target_pointer;
  }
  else {
    if (COO_size == 0) {
      for (int i=0; i<accumulator_size; i++) {
        COO_crd[0][i] = accumulator[i].crd[0];
        COO_crd[1][i] = accumulator[i].crd[1];
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
      tmp_con.crd[0] = COO_crd[0][content_pointer];
      tmp_con.crd[1] = COO_crd[1][content_pointer];
      if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) == 0) {
        tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
        tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
        tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val + COO_vals[content_pointer];
        accumulator_pointer ++;
        content_pointer ++;
        target_pointer ++;
      } else if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) < 0) {
        tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
        tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
        tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
        accumulator_pointer ++;
        target_pointer ++;
      } else {
        tmp_COO_crd[0][target_pointer] = COO_crd[0][content_pointer];
        tmp_COO_crd[1][target_pointer] = COO_crd[1][content_pointer];
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
      tmp_COO_crd[0][target_pointer] = COO_crd[0][content_pointer];
      tmp_COO_crd[1][target_pointer] = COO_crd[1][content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      content_pointer ++;
      target_pointer ++;
    }
    for (int i = 0; i < target_pointer; i++) {
      COO_crd[0][i] = tmp_COO_crd[0][i];
      COO_crd[1][i] = tmp_COO_crd[1][i];
      COO_vals[i] = tmp_COO_vals[i];
    }
    free(tmp_COO_crd[0]);
    free(tmp_COO_crd[1]);
    free(tmp_COO_vals);
    return target_pointer;
  }
}

int Sort(void * array, size_t size, bool rev) {
  if (rev) {
    qsort(array, size, sizeof(wspace), esc_cmp_rev);
    for (int i = size - 1; i >= 0; i--) {
      if(((wspace*)array)[i].crd[0] != -1) {
        return i + 1;
      }
    }
  } 
  qsort(array, size, sizeof(wspace), esc_cmp);
  return size;
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

  const int w_accumulator_capacity = 3; // Arbitrary
  int w_accumulator_size = 0;
  wspace* restrict w_accumulator = 0; // hash accumulator
  w_accumulator = (wspace*)malloc(sizeof(wspace) * w_accumulator_capacity);
  for (int32_t i = 0; i<w_accumulator_capacity; i++) {
    init_wspace(&w_accumulator[i]);
  }

  int32_t w_all_capacity = w_accumulator_capacity; // Ensure enough space for merging
  int32_t w_all_size[2] = {0,0};
  wspace* restrict w_all[2];
  w_all[0] = (wspace*)malloc(sizeof(wspace) * w_all_capacity);
  w_all[1] = (wspace*)malloc(sizeof(wspace) * w_all_capacity);
  int w_all_current = 0;// Current content wspace array, use ^ to get target wspace array
  print_wspace(w_accumulator,w_accumulator_capacity);
  for (int32_t i = 0; i < A1_dimension; i++) {
    for (int32_t jA = A2_pos[i]; jA < A2_pos[i+1]; jA++) {
      int32_t j = A2_crd[jA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[j+1]; kB++) {
        int32_t k = B2_crd[kB];
        int w_hashkey = hash_func(i * B2_dimension + k, w_accumulator_capacity);
        int w_step = 0;
        wspace w_tmp;
        w_tmp.crd[0] = i;
        w_tmp.crd[1] = k;
        w_tmp.val = A_vals[jA] * B_vals[kB];
        cout<<i<<","<<k<<","<<w_hashkey<<endl;
        // Try insertion to the hash array
        do{
          if (w_accumulator[w_hashkey].crd[0] == -1){ // no conflict
            w_accumulator[w_hashkey] = w_tmp;
            w_accumulator_size ++;
            break;
          } else if (esc_cmp(&w_accumulator[w_hashkey],&w_tmp) == 0) {
            w_accumulator[w_hashkey].val += w_tmp.val;
            break;
          } else {
            w_hashkey = (w_hashkey + 1) % w_accumulator_capacity; // linear-probe
            w_step ++;
          }
        }while(w_step < w_accumulator_capacity);
        print_wspace(w_accumulator,w_accumulator_capacity);

        // No space to insert. Do contraction
        if (w_step == w_accumulator_capacity) {
          // Sort
          qsort(w_accumulator, w_accumulator_size, sizeof(wspace), esc_cmp);
          
          // Enlarge output and input arraies (optimization needed)
          if(w_accumulator_size + w_all_size[w_all_current] > w_all_capacity) {
              w_all_capacity = w_all_capacity * 2;
              w_all[0] = (wspace*)realloc(w_all[0], sizeof(wspace) * w_all_capacity);
              w_all[1] = (wspace*)realloc(w_all[1], sizeof(wspace) * w_all_capacity);
          }

          // Clear the target array
          w_all_size[w_all_current ^ 1] = 0;

          // Warm up
          if (w_all_size[w_all_current] == 0) {
            for (int wi=0; wi<w_accumulator_size; wi++) {
              w_all[w_all_current ^ 1][wi] = w_accumulator[wi];
              w_all_size[w_all_current ^ 1] ++;
            }
          }

          // Merge sort
          else { 
            int w_accumulator_pointer = 0;
            int w_content_pointer = 0;
            int w_target_pointer = 0;
            wspace* w_content = w_all[w_all_current];
            wspace* w_target = w_all[w_all_current ^ 1];
            int w_current_size = w_all_size[w_all_current];

            while(w_accumulator_pointer<w_accumulator_size && w_content_pointer<w_current_size) {
              if(esc_cmp(&w_accumulator[w_accumulator_pointer],&w_content[w_content_pointer]) == 0) {
                w_target[w_target_pointer] = w_accumulator[w_accumulator_pointer];
                //init_wspace(&w_accumulator[w_accumulator_pointer]); // Clear the accumulator after usage
                w_target[w_target_pointer].val += w_content[w_content_pointer].val;
                w_accumulator_pointer ++;
                w_content_pointer ++;
                w_target_pointer ++;
              } else if(esc_cmp(&w_accumulator[w_accumulator_pointer],&w_content[w_content_pointer]) < 0) {
                w_target[w_target_pointer] = w_accumulator[w_accumulator_pointer];
                //init_wspace(&w_accumulator[w_accumulator_pointer]);
                w_accumulator_pointer ++;
                w_target_pointer ++;
              } else {
                w_target[w_target_pointer] = w_content[w_content_pointer];
                w_content_pointer ++;
                w_target_pointer ++;
              }
            }
            while(w_accumulator_pointer<w_accumulator_size) {
              w_target[w_target_pointer] = w_accumulator[w_accumulator_pointer];
              //init_wspace(&w_accumulator[w_accumulator_pointer]);
              w_accumulator_pointer ++;
              w_target_pointer ++;
            }
            while(w_content_pointer<w_current_size) {
              w_target[w_target_pointer] = w_content[w_content_pointer];
              w_content_pointer ++;
              w_target_pointer ++;
            }
            w_all_size[w_all_current ^ 1] = w_target_pointer;
          }
          // Change the content index
          w_all_current ^= 1;
          w_accumulator_size = 0;
          // Clear the hash array
          for (int32_t wi = 0; wi<w_accumulator_capacity; wi++) {
            init_wspace(&w_accumulator[wi]);
          }
          // Insert the wspace that conflicts
          w_accumulator[hash_func(i * B2_dimension + k, w_accumulator_capacity)] = w_tmp;
          w_accumulator_size ++;

        }
      }
    }
  }
  if (w_accumulator_size > 0) {
    qsort(w_accumulator, w_accumulator_size, sizeof(wspace), esc_cmp);
    // Enlarge output and input arraies (optimization needed)
    if(w_accumulator_capacity + w_all_size[w_all_current] > w_all_capacity) {
        w_all_capacity = w_all_capacity * 2;
        w_all[0] = (wspace*)realloc(w_all[0], sizeof(wspace) * w_all_capacity);
        w_all[1] = (wspace*)realloc(w_all[1], sizeof(wspace) * w_all_capacity);
    }

    // Clear the target array
    w_all_size[w_all_current ^ 1] = 0;

    // Warm up
    if (w_all_size[w_all_current] == 0) {
      for (int wi=0; wi<w_accumulator_capacity; wi++) {
        w_all[w_all_current ^ 1][wi] = w_accumulator[wi];
        w_all_size[w_all_current ^ 1] ++;
      }
    }

    // Merge sort
    else { 
      int w_accumulator_pointer = 0;
      int w_content_pointer = 0;
      int w_target_pointer = 0;
      //merge(current w_all, w_accumulator) -> new w_all
      wspace* w_content = w_all[w_all_current];
      wspace* w_target = w_all[w_all_current ^ 1];
      int w_current_size = w_all_size[w_all_current];

      while(w_accumulator_pointer<w_accumulator_size && w_content_pointer<w_current_size) {
        if(esc_cmp(&w_accumulator[w_accumulator_pointer],&w_content[w_content_pointer]) == 0) {
          w_target[w_target_pointer] = w_accumulator[w_accumulator_pointer];
          w_target[w_target_pointer].val += w_content[w_content_pointer].val;
          w_accumulator_pointer ++;
          w_content_pointer ++;
          w_target_pointer ++;
        } else if(esc_cmp(&w_accumulator[w_accumulator_pointer],&w_content[w_content_pointer]) < 0) {
          w_target[w_target_pointer] = w_accumulator[w_accumulator_pointer];
          w_accumulator_pointer ++;
          w_target_pointer ++;
        } else {
          w_target[w_target_pointer] = w_content[w_content_pointer];
          w_content_pointer ++;
          w_target_pointer ++;
        }
      }
      while(w_accumulator_pointer<w_accumulator_size) {
        w_target[w_target_pointer] = w_accumulator[w_accumulator_pointer];
        w_accumulator_pointer ++;
        w_target_pointer ++;
      }
      while(w_content_pointer<w_current_size) {
        w_target[w_target_pointer] = w_content[w_content_pointer];
        w_content_pointer ++;
        w_target_pointer ++;
      }
      w_all_size[w_all_current ^ 1] = w_target_pointer;
    }
    // Change the content index
    w_all_current ^= 1;
  }

  int w_current_size = w_all_size[w_all_current];
  wspace* w_content = w_all[w_all_current];
  int C_COO1_pos[2] = {0,w_current_size};
  int32_t* C_COO_crd[2];
  for (int wi = 0; wi < w_order; wi++) {
    C_COO_crd[wi] = (int32_t*)malloc(sizeof(int32_t) * w_current_size);
    for (int wk = 0; wk < w_current_size; wk++) {
      C_COO_crd[wi][wk] = w_content[wk].crd[wi];
    }
  }
  float* C_COO_vals = (float*)malloc(sizeof(float) * w_current_size);
  for (int wk = 0; wk < w_current_size; wk++) {
    C_COO_vals[wk] = w_content[wk].val;
  }

  pack_C(C, C_COO1_pos, C_COO_crd[0], C_COO_crd[1], C_COO_vals);
  return 0;
}

int32_t TryInsert_oldhash(bool* insertFail, wspace* accumulator, int32_t accumulator_size, int32_t accumulator_capacity, int32_t* crds, float val, int hashkey) {
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
    } else if (esc_cmp(&accumulator[hashkey],&tmp) == 0) {
      accumulator[hashkey].val += tmp.val;
      *insertFail = false;
      return accumulator_size;
    } else {
      hashkey = (hashkey + 1) % accumulator_capacity; // linear-probe
      step ++;
    }
  }while(step < accumulator_capacity);
  *insertFail = true;
  return accumulator_size;
}

int32_t TryInsert(bool* insertFail, wspace* accumulator, int32_t accumulator_size, int32_t accumulator_capacity, int32_t* crds, float val) {
  int hashkey = hash_func(crds[0] + crds[1], accumulator_capacity);
  int step = 0;
  int bmap = 0;
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
  }while(bmap & 1<<hashkey != 0);
  *insertFail = true;
  return Sort(accumulator, accumulator_capacity, true);
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

  int32_t C_COO_capacity = w_accumulator_capacity; // Ensure enough space for merging
  int32_t* C_COO_crd[2];
  float* C_COO_vals;
  int32_t C_COO_size = 0;
  C_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * C_COO_capacity);
  C_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * C_COO_capacity);
  C_COO_vals = (float*)malloc(sizeof(float) * C_COO_capacity);
  bool insertFail;
  // Initialize hash key
  for (int32_t i = 0; i<w_accumulator_capacity; i++) {
    init_wspace(&w_accumulator[i]);
  }

  for (int32_t i = 0; i < A1_dimension; i++) {
    for (int32_t jA = A2_pos[i]; jA < A2_pos[i+1]; jA++) {
      int32_t j = A2_crd[jA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[j+1]; kB++) {
        int32_t k = B2_crd[kB];
        int32_t C_crds[2] = {i,k};
        // Try to insert to the hash array
        int w_hashkey = hash_func(i * B2_dimension + k, w_accumulator_capacity);
        w_accumulator_size = TryInsert_oldhash(&insertFail, w_accumulator, w_accumulator_size, w_accumulator_capacity, C_crds, A_vals[jA] * B_vals[kB], w_hashkey);
        // No space to insert. Do contraction
        if (insertFail) {
          // Sort
          // Enlarge output and input arraies (optimization needed)
          if(w_accumulator_size + C_COO_size > C_COO_capacity) {
            C_COO_capacity = Enlarge(C_COO_crd,&C_COO_vals,C_COO_capacity);
          }
          // Merge
          C_COO_size = merge_coo(C_COO_crd, C_COO_vals, C_COO_size, w_accumulator, w_accumulator_size);
          // Clear
          w_accumulator_size = 0;
          // Clear the hash array
          for (int32_t wi = 0; wi<w_accumulator_capacity; wi++) {
            init_wspace(&w_accumulator[wi]);
          }
          // Insert the wspace that conflicts
          w_accumulator_size = TryInsert_oldhash(&insertFail, w_accumulator, w_accumulator_size, w_accumulator_capacity, C_crds, A_vals[jA] * B_vals[kB], w_hashkey);
        }
      }
    }
  }
  if (w_accumulator_size > 0) {
    // Sort
    qsort(w_accumulator, w_accumulator_size, sizeof(wspace), esc_cmp);
    // Enlarge
    if(w_accumulator_size + C_COO_size > C_COO_capacity) {
      C_COO_capacity = Enlarge(C_COO_crd,&C_COO_vals,C_COO_capacity);
    }
    // Merge
    C_COO_size = merge_coo(C_COO_crd, C_COO_vals, C_COO_size, w_accumulator, w_accumulator_size);
    // Clear
    w_accumulator_size = 0;
  }
  int C_COO1_pos[2] = {0,C_COO_size};
  pack_C(C, C_COO1_pos, C_COO_crd[0], C_COO_crd[1], C_COO_vals);
  return 0;
} 

int compute_coo_rev(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
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

  int32_t C_COO_capacity = w_accumulator_capacity; // Ensure enough space for merging
  int32_t* C_COO_crd[2];
  float* C_COO_vals;
  int32_t C_COO_size = 0;
  C_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * C_COO_capacity);
  C_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * C_COO_capacity);
  C_COO_vals = (float*)malloc(sizeof(float) * C_COO_capacity);
  bool insertFail;
  // Initialize hash key
  for (int32_t i = 0; i<w_accumulator_capacity; i++) {
    init_wspace(&w_accumulator[i]);
  }

  for (int32_t i = 0; i < A1_dimension; i++) {
    for (int32_t jA = A2_pos[i]; jA < A2_pos[i+1]; jA++) {
      int32_t j = A2_crd[jA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[j+1]; kB++) {
        int32_t k = B2_crd[kB];
        int32_t C_crds[2] = {i,k};
        // Try to insert to the hash array
        w_accumulator_size = TryInsert(&insertFail, w_accumulator, w_accumulator_size, w_accumulator_capacity, C_crds, A_vals[jA] * B_vals[kB]);
        print_wspace(w_accumulator,w_accumulator_capacity);
        // No space to insert. Do contraction
        if (insertFail) {
          // Enlarge output and input arraies (optimization needed)
          if(w_accumulator_size + C_COO_size > C_COO_capacity) {
            C_COO_capacity = Enlarge(C_COO_crd,&C_COO_vals,C_COO_capacity);
          }
          // Merge
          C_COO_size = Merge(C_COO_crd, C_COO_vals, C_COO_size, w_accumulator, w_accumulator_size, true);
          print_coo(C_COO_crd, C_COO_vals, C_COO_size);
          // Clear
          w_accumulator_size = 0;
          // Clear the hash array
          for (int32_t wi = 0; wi<w_accumulator_capacity; wi++) {
            init_wspace(&w_accumulator[wi]);
          }
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


void CSR_CSR_5(const string& A_name, const string& B_name) {
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
    compute_coo_rev(&C,&A,&B);
    print_taco_tensor_DC(&A);
    print_taco_tensor_DC(&B);
    print_taco_tensor_DC(&C);
}