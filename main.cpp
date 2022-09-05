//#include "seperate/CSR_CSC_0.h" // failed
#include "seperate/CSR_CSR_4.h"
//#include "evaluate/CSR_CSR_1.h"
//#include "seperate/DCSR_DCSR_4.h"
//#include "seperate/DCSR_DCSC_4.h"
//#include "evaluate/DCSR_DCSC_1.h"
//#include "seperate/CSR_CSC_4.h"
//#include "seperate/CSR_CSR_5.h"
#include <iostream>
using namespace std;
int main() {
    string filename = "./data/test5.mtx";
    CSR_CSR_4(filename, filename);
    return 0;    
}