// simplex.h

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <iostream>
#include <chrono>
#include <matrix.h>
#include <vector>

enum ConstraintType {LEQ, EQ, GEQ};

class Simplex {
    private:
        int m; // amount of constraints, rows
        int n; // amount of original svariables, columns
        int n_big; // amount of variables in big-Ms
        int nb;
        int mnrc_idx;
        int leaving_row;

        int n_leq = 0;
        int n_geq = 0;
        int n_eq = 0;

        bool infeasible = false;
        bool DEBUG = false;

        std::vector<int> basis_idx;
        std::vector<int> non_basis_idx;


        matrix A;
        matrix b;
        matrix c_trans;
        matrix A_big;
        matrix c_trans_big;
        matrix x;
       
        matrix B;
        matrix B_inv;
        matrix c_B;
        matrix c_N;
        matrix N;
        matrix a_j;
        matrix a_j_hat;
        std::vector<int> artificial_indices;

        double M = 1e9; // Big-M value
        
        std::vector<ConstraintType> constraint_types;

        mutable Matrix mops;
        const int MAX_ITERATIONS = 1000;
        // functions

        void get_basis();
        void calculate_current_solution();
        void print_solution() const;
        void sherman_morrison(int leaving_row, int entering_column, const matrix& a_j_hat);
        void construct_b(const matrix& b_in, const std::vector<ConstraintType>& constraint_types_in);

    public:
        // constructor
        Simplex(const matrix& A_in,
                const matrix& b_in,
                const matrix& c_trans_in,
                const std::vector<ConstraintType>& constraint_types_in);
        
        // Solves the optimization problem
        void solve();
        matrix return_solution() const;
        double return_opt_cost() const;



};


#endif // SIMPLEX_H