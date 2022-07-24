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

void print_wspace(wspace* w, const int len) {
  for(int i=0; i<len; i++) {
    cout<<"("<<w[i].crd[0]<<","<<w[i].crd[1]<<","<<w[i].val<<"),";
  }
  cout<<endl;
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
          w_accumulator_size = 0;
          // Sort
          qsort(w_accumulator, w_accumulator_capacity, sizeof(wspace), esc_cmp);
          
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
            wspace* w_content = w_all[w_all_current];
            wspace* w_target = w_all[w_all_current ^ 1];
            int w_current_size = w_all_size[w_all_current];

            while(w_accumulator_pointer<w_accumulator_capacity && w_content_pointer<w_current_size) {
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
            while(w_accumulator_pointer<w_accumulator_capacity) {
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
      wspace* w_content = w_all[w_all_current];
      wspace* w_target = w_all[w_all_current ^ 1];
      int w_current_size = w_all_size[w_all_current];

      while(w_accumulator_pointer<w_accumulator_capacity && w_content_pointer<w_current_size) {
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
      while(w_accumulator_pointer<w_accumulator_capacity) {
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
    compute(&C,&A,&B);
    print_taco_tensor_DC(&A);
    print_taco_tensor_DC(&B);
    print_taco_tensor_DC(&C);
}