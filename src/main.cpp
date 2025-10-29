#include <iostream>
//#include "pid.h"
#include "matrix.h"

int main() {
    
    Matrix mops;
    matrix A = mops.zeros(2, 2);

    mops.print(A);

    A(0,0) = 1.0;
    A(1,1) = 1.0;


    double num = A(0,0);

    std::cout << "first element of matrix A: " << num << "\n";

    mops.print(A);

    

    matrix B = mops.zeros(2,2);

    B(0,0) = 2.0;
    B(1,1) = 2.0;

    matrix C = A * B;

    std::cout << " A * B = C: " << "\n";
    mops.print(C);

    matrix D = A + B;

    std::cout << " A + B = D: " << "\n";
    mops.print(D);
    std::cout << "test " << 2*3 << "\n";

    delete[] A.data;
    delete[] C.data;
    delete[] B.data;


    return 0;

}