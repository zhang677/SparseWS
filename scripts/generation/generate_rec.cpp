#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <time.h>
#include <chrono>
#include <vector>
#include <numeric>


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

int main(int argc, char** argv) {
    // Take in arguments from command line
    const int row = (argc > 1) ? stoi(argv[1]) : 25;
    const int col = (argc > 2) ? stoi(argv[2]) : 200;
    const int nz_row = (argc > 3) ? stoi(argv[3]) : 5;
    const int nz_col = (argc > 4) ? stoi(argv[4]) : 10;
    const string folder = (argc > 5) ? argv[5] : "/home/zgh23/code/SparseWS/data/origin/rec/";
    const int seed = (argc > 6) ? stoi(argv[6]) : 42;

    int nnz = nz_row * nz_col;
    double density = 1.0 * nnz / (row * col);
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << density;
    std::ofstream file;
    std::cout << "Generating..." << std::endl;
    file.open(folder + "rec_" + std::to_string(nz_row) + "_" + std::to_string(row) + "_" + std::to_string(nz_col) + "_" + std::to_string(col) + "_" + std::to_string(seed) + ".mtx", std::ios::trunc);
    file << "%%MatrixMarket matrix coordinate real general" << std::endl;
    file << row << " " << col << " " << nnz << std::endl;
    std::cout << row << " " << col << " " << nnz << std::endl;
    Timer timer;
    timer.reset();
    // Set random seed to 42
    srand(seed);
    // Randomly pick nz_row from range [0, row - 1] without replacement
    std::vector<int> col_idx(col);
    std::iota(col_idx.begin(), col_idx.end(), 0);
    std::random_shuffle(col_idx.begin(), col_idx.end());
    col_idx.resize(nz_col);
    std::sort(col_idx.begin(), col_idx.end());
    std::vector<int> row_idx(row);
    for (int j = 0; j < nz_col; j++) {
        // Randomly pick nz_col from range [0, col - 1] without replacement
        std::iota(row_idx.begin(), row_idx.end(), 0);
        std::random_shuffle(row_idx.begin(), row_idx.end());
        row_idx.resize(nz_row);
        std::sort(row_idx.begin(), row_idx.end());
        for (int i = 0; i < nz_row; i++) {
            file << row_idx[i] + 1 << " " << col_idx[j] + 1 << " " << 1.0 * rand() / RAND_MAX << std::endl;
        }
    }
    file.close();
    std::cout << "Time elapsed: " << timer.elapsed() << "s" << std::endl;
    return 0;
}