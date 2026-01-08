// main.cpp


#include <iostream>
#include <vector> 
#include "matrix.h"


int main() {
    Matrix mops;
    matrix A(2, 2);
    A(0,0) = 1;
    A(1,0) = 2;
    A(0,1) = 3;
    A(1,1) = 4;

    mops.print(A);
    matrix At(2, 2); // How to egt rid of these ??
    At = A.transpose();
    std::cout << "Transposed:" << std::endl;
    mops.print(At);
    matrix A2(2, 2); // How to egt rid of these ??
    A2 = A * A;
    std::cout << "Multiplied:" << std::endl;
    mops.print(A2);
    return 0;
}