#include <iostream>
#include <vector> 
#include <iomanip>
#include "matrix.h"
#include "simplex.h" 
Matrix mops;


void run_test(const char* name, const matrix& A, const matrix& b, const matrix& c, 
              const std::vector<ConstraintType>& types, bool maxim) {
    std::cout << "\n";
    std::cout << "==================================================\n";
    std::cout << "TEST: " << name << "\n";
    std::cout << "==================================================\n";

    // Call the updated constructor
    Simplex solver(A, b, c, types); 
    solver.solve();
    //
    if (maxim) {
         matrix res = solver.return_solution();
        mops.print(res);
        double solution = solver.return_opt_cost() * (-1.0);
        std::cout << "optimal value: " << solution << "\n";
    } else {
        solver.print_solution();
    } 
   


    
}


int main() {
    std::cout << std::setprecision(15); 

    matrix A = mops.zeros(2,2);
    A(0,0)=4.5; A(0,1)=3.2;
    A(1,0)=1.3; A(1,1)=2.4;

    matrix b = mops.zeros(2,1);
    b(0,0)=70.0; b(1,0)=30.0;

    matrix c = mops.zeros(1,2);
    c(0,0)=-65.0; c(0,1)=-48.0;

    std::vector<ConstraintType> types = {LEQ, LEQ};

    run_test("lecture notes", A, b, c, types, true);
    return 0;
}