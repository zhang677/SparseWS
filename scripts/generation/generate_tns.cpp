#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <time.h>
#include <chrono>
#include <vector>
#include <numeric>
#include <functional>
#include <cstring>

using namespace std;

class Timer {
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count(); }
private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

std::vector<std::string> split(const std::string &str, const std::string &delim, bool keepDelim) {
  std::vector<std::string> results;
  size_t prev = 0;
  size_t next = 0;

  while ((next = str.find(delim, prev)) != std::string::npos) {
    if (next - prev != 0) {
      std::string substr = ((keepDelim) ? delim : "")
                         + str.substr(prev, next-prev);
      results.push_back(substr);
    }
    prev = next + delim.size();
  }

  if (prev < str.size()) {
    string substr = ((keepDelim) ? delim : "") + str.substr(prev);
    results.push_back(substr);
  }

  return results;
}

int main(int argc, char** argv) {
    std::string tensorDims = "1000_1000_1000", tensorNnzs = "100_100_100", folder = "/home/zgh23/code/SparseWS/data/origin/tns/";
    for (int i = 1; i < argc; i++) {
#define INT_ARG(argname, varname) do {      \
          if (!strcmp(argv[i], (argname))) {  \
            varname = atoi(argv[++i]);      \
            continue;                       \
          } } while(0);
#define STRING_ARG(argname, varname) do {      \
          if (!strcmp(argv[i], (argname))) {  \
            varname = std::string(argv[++i]);      \
            continue;                       \
          } } while(0);
    STRING_ARG("-nnzs", tensorNnzs);
    STRING_ARG("-dims", tensorDims);
    STRING_ARG("-folder", folder);
#undef INT_ARG
#undef STRING_ARG
  }
    auto dimsStr = split(tensorDims, "_", false /* keepDelim */);
    std::vector<int> dims;
    for (auto it : dimsStr) {
        dims.push_back(atoi(it.c_str()));
    }
    auto nnzsStr = split(tensorNnzs, "_", false /* keepDelim */);
    std::vector<int> nnzs;
    for (auto it : nnzsStr) {
        nnzs.push_back(atoi(it.c_str()));
    }
    int nnz = 1;
    int modes = dims.size();
    for (int i = 0; i < modes; i++) {
        nnz *= nnzs[i];
    }
    int coords = 1;
    for (int i = 0; i < modes; i++) {
        coords *= dims[i];
    }
    double density = 1.0 * nnz / coords;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << density;
    std::ofstream file;
    std::cout << "Generating..." << std::endl;
    file.open(folder + "cube_" + tensorDims + "_" + tensorNnzs + "_" + oss.str() + ".tns", std::ios_base::app);
    srand(42);

    std::vector<std::vector<int>> idxs(modes);
    for (int i = 0; i < modes; i++) {
        idxs[i].resize(dims[i]);
    }
    // A recursive function to generate the coordinates. Mode i Randomly pick nnzs[i] coordinates from range [0, modes[i] - 1] without replacement
    std::vector<int> coord(modes);
    std::function <void(int, std::vector<int>&)> gen_coords = [&](int i, std::vector<int>& coord) {
        std::iota(idxs[i].begin(), idxs[i].end(), 0);
        std::random_shuffle(idxs[i].begin(), idxs[i].end());
        idxs[i].resize(nnzs[i]);
        std::sort(idxs[i].begin(), idxs[i].end());
        if (i == modes - 1) {
            for (int j = 0; j < nnzs[i]; j++) {
                coord[i] = idxs[i][j];
                for (int k = 0; k < modes; k++) {
                    file << coord[k] + 1 << " ";
                }
                file << 1.0 * rand() / RAND_MAX << std::endl;
            }
            coord.clear();
        } else {
            for (int j = 0; j < nnzs[i]; j++) {
                coord[i] = idxs[i][j];
                gen_coords(i + 1, coord);
            }
        }
    };
    gen_coords(0, coord);

}