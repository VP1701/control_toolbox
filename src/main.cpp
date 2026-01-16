// main.cpp


#include <iostream>
#include <vector> 
#include "matrix.h"
#include "simplex.h"


void run_test(const char* name, const matrix& A, const matrix& b, const matrix& c, 
              const std::vector<ConstraintType>& types) {
    std::cout << "\n";
    std::cout << "==================================================\n";
    std::cout << "TEST: " << name << "\n";
    std::cout << "==================================================\n";

    // Call the updated constructor
    
    
    Simplex solver(A, b, c, types); 
    
    // solver.solve() is called in constructor

    
}

int main() {
    Matrix mops;
    

    {
        matrix A = mops.zeros(2,2);
        A(0,0)=2; A(0,1)=3;
        A(1,0)=-1; A(1,1)=1;

        matrix b = mops.zeros(2,1);
        b(0,0)=6; b(1,0)=1;

        matrix c = mops.zeros(1,2);
        c(0,0)=-1; c(0,1)=-3;

        std::vector<ConstraintType> types = {LEQ, LEQ};

        run_test("lecture notes", A, b, c, types);
        std::cout << "answer should be: " << "x_transpose = [0.6, 1.6], opt_cost = -5.4" << "\n";
    }

    return 0;
}