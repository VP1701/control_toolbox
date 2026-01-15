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
    matrix At(2, 2);
    At = A.T();
    std::cout << "Transposed:" << std::endl;
    mops.print(At);
    matrix A2(2, 2);
    A2 = A * A;
    std::cout << "Multiplied:" << std::endl;
    mops.print(A2);
    std::cout << "Row 2 multiplied by 2:" << std::endl;
    mops.print(A2.multiply_row(1, 2));
    std::cout << "swap rows 1 and 2:" << std::endl;
    mops.print(A2.swap_rows(1, 0));
    std::cout << "add 2 times row1 to row 2" << std::endl;
    mops.print(A2.add_multiple_of_row(1, 0, 2));

    return 0;
}