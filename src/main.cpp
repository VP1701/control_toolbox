// main.cpp


#include <iostream>
#include <vector> 
#include "matrix.h"
#include "simplex.h"
#include "mpc.h"

Matrix mops;

void run_test(const char* name, const matrix& A, const matrix& b, const matrix& c, 
              const std::vector<ConstraintType>& types) {
    std::cout << "\n";
    std::cout << "==================================================\n";
    std::cout << "TEST: " << name << "\n";
    std::cout << "==================================================\n";

    // Call the updated constructor
    
    
    Simplex solver(A, b, c, types); 
    int iters = 1;
    long long time_total;
    
    for (int i = 0; i < iters; ++i) {
        typedef std::chrono::high_resolution_clock Clock;
        auto t1 = Clock::now();
        solver.solve();
        auto t2 = Clock::now();
        auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1);
        time_total += ns.count();
        std::cout << "iteration" << i << "\n";
    }
    double avg_ns = time_total / iters;
    double avg_ms = avg_ns / 1000000.0;
    std::cout << "time used for solving "<< iters << " LP problems: " << time_total / 1000000.0 << " ms\n";
    std::cout << "average time for on LP problem: " << avg_ms << " ms\n";
    matrix sol(3,1);
    sol = solver.return_solution();
    std::cout << "printing soltuion outside class" << "\n";
    mops.print(sol);
    
}

int main() {
    
    

    {

        // model of 
        /*
        matrix A(3,3);
        A(0,0) = 1;  A(0,1) = 0.0002;  A(0,2) = 0.00003;
        A(1,0) = 0;  A(1,1) = 1;  A(1,2) = 0.1726;
        A(2,0) = 0;  A(2,1) = 0;  A(2,2) = 1;

        matrix B = mops.zeros(3,1);
        B(0,0)=0.00003; 
        B(1,0)=0.1726;
        B(2,0)=0.0;

        matrix C = mops.zeros(1,3);
        C(0,0)=1; C(0,1)=0; C(0,2)=0;

        matrix x0(3,1);
        x0(0,0)=2; x0(1,0)=2; x0(2,0)=1;
        */
        std::vector<ConstraintType> types = {LEQ, LEQ, LEQ};
        
        matrix A(2,2);
        A(0,0) = 0.9;  A(0,1) = -0.5;
        A(1,0) = 0;  A(1,1) = 0.8;
        

        matrix B = mops.zeros(2,2);
        B(0,0)= 1.0; B(0,1) = 1;
        B(1,0)= -2.0; B(1,1) = 0.0;
        

        matrix C = mops.zeros(2,2);
        C(0,0)=2; C(0,1)=0.5; 
        C(1,0) = -1.0; C(1,1) = 1.0;
        matrix W_x = mops.eye(4);
        matrix W_u = mops.eye(4) * 2;
        matrix W_du = mops.eye(4) * 3;

        matrix du_max(2,1);
        du_max(0,0) = 1.0;
        du_max(1,0) = 2.0;

        matrix du_min(2,1);
        du_min(0,0) = -1.0;
        du_min(1,0) = -3.0;

        matrix u_max(2,1);
        u_max(0,0) = 1.0;
        u_max(1,0) = 4.0;

        matrix u_min(2,1);
        u_min(0,0) = -1.0;
        u_min(1,0) = -5.0;

        matrix x_max(2,1);
        x_max(0,0) = 1.0;
        x_max(1,0) = 6.0;

        matrix x_min(2,1);
        x_min(0,0) = -1.0;
        x_min(1,0) = -7.0;
        matrix x0(2,1);
        x0(0,0) = 0.0;
        x0(1,0) = 0.0;

        
        MPC test_solver(A, B, C, W_x, W_u, W_du,
                        du_max, du_min, u_max, 
                        u_min, x_max, x_min, x0, 4);

        
        /*run_test("lecture notes", A, B, C, types);
        std::cout << "answer should be: " << "x_transpose = [0.6, 1.6], opt_cost = -5.4" << "\n";
        std::cout << "A = " << "\n";
        mops.print(A);
        std::cout << "A + A = " << "\n";
        mops.print(A+A);
        matrix D = A * C;
        std::cout << "A * x = " << "\n";
        mops.print(A*x);
        std::cout << "C * x = " << "\n";
        mops.print(C*x);
        std::cout << "C= " << "\n";
        mops.print(C);

        std::cout << "C * B = " << "\n";
        mops.print(C*B);
        std::cout << "C * A * B = " << "\n";
        mops.print(C*A*B);
        std::cout << "C * A * A * B = " << "\n";
        mops.print(C*A*A*B);
        MPC test(A,B,C,3);
        */

        /*int time_steps = 100;
        matrix u(1,1);
        u(0,0)=0.1;
        for (int i = 0; i < time_steps; i++) {
            x0 = A*x0+B*u;
            std::cout << "state value: " << "\n";
            mops.print(x0);
        }
        
    
    }

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
    }*/
    }

    return 0;
}