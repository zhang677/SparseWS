//#include "seperate/CSR_CSC_0.h"
#include "seperate/CSR_CSR_4.h"
//#include "include/CSR_CSR_1.h"
//#include "include/DCSR_CSR_1.h"
//#include "include/DCSR_DCSC_0.h"
//#include "include/DCSR_DCSC_1.h"
//#include "include/DCSR_DCSR_1.h"
//#include "CSR_CSR_ESC.h"
#include <iostream>
using namespace std;
int main() {
    string filename = "./data/test.mtx";
    CSR_CSR_4(filename, filename);
    return 0;    
}