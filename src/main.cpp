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
    return 0;
}