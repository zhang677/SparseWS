#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <time.h>
#include <chrono>

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
    // arg0: row
    // arg1: band width
    // arg2: output file name
    const int row = (argc > 1) ? stoi(argv[1]) : 10;
    const int band_width = (argc > 2) ? stoi(argv[2]) : 1;
    const string folder = (argc > 3) ? argv[3] : "/home/zgh23/code/SparseWS/data/origin/band/";
    // Generate a band matrix
    int nnz = row * (2 * band_width - 1) - band_width * (band_width - 1); 
    double density = 1.0 * nnz / row / row;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << density;
    std::ofstream file;
    std::cout << "Generating..." << std::endl;
    file.open(folder + "band_" + std::to_string(row) + "_" + std::to_string(band_width) + "_" + oss.str() + ".mtx", std::ios_base::app);
    file << "%%MatrixMarket matrix coordinate real general" << std::endl;
    file << row << " " << row << " " << nnz << std::endl;
    std::cout << row << " " << row << " " << nnz << std::endl;
    Timer timer;
    timer.reset();
    // Generate a random band matrix
    for (int i = 0; i < row; i++) {
        for (int j = std::max(0, i - band_width + 1); j < std::min(row, i + band_width); j++) {
            file << (i + 1) << " " << (j + 1) << " " << (double) rand() / RAND_MAX << std::endl; // mtx format is 1-based
        }
    }
    file.close();
    // Write to file
    std::cout << "Time elapsed: " << timer.elapsed() << "s" << std::endl;
    return 0;
}