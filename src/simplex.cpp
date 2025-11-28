// simplex.cpp

#include <algorithm>
#include <iostream>
#include <simplex.h>
#include <vector>

void Simplex::get_basis() {

    
    for (int i = 0; i < m; ++i) {
        int A_col = basis_idx[i];
        for (int j = 0; j < m; ++j) {
            B(j,i) = A(j, A_col);
        }
    }

    for (int i = 0; i < nb; ++i) {
        int A_col = non_basis_idx[i];
        for (int j = 0; j < m; ++j) {
            N(j,i) = A(j, A_col);
        }
    }


}

void Simplex::calculate_current_solution() {
    // calculate xB = B^(-1) * b

    get_basis();
    B_inv = mops.inverse(B);
    matrix xB = B_inv * b;
    x = mops.zeros(n_big,1);
    for (int i = 0; i < m; ++i) {
        int idx = basis_idx[i];
        x(idx, 0) = xB(i, 0);
    }

}

void Simplex::print_solution() const {

    std::cout << "Printing solution to the simplex" << "\n";


    mops.print(x);
    
    matrix opt = c_trans * x;
    double opt_val = opt(0,0);
    std::cout << "Optimal value: " << opt_val <<"\n";
    //mops.print(opt);
}

void Simplex::solve() {
    // todo
    std::cout << "Solver started!" << "\n";
    for (int i = 0; i < MAX_ITERATONS; ++i) {
        calculate_current_solution();

        for (int i = 0; i < m; ++i) {
            c_B(0,i) = c_trans(0, basis_idx[i]);
            
        }

        for (int i = 0; i < nb; ++i) {
            c_N(0,i) = c_trans(0, non_basis_idx[i]);
        }

        

        //reduced_cost = c_N - c_B*B_inv*N
        matrix reduced_cost = c_N - c_B*B_inv*N;
        //std::cout << "Reduced cost" << "\n";
        mops.print(reduced_cost);

        double mnrc = -1e-12; // most negative reduces cost. start with barely negative value
        mnrc_idx = -1; // most negative reduces cost index

        // check if all values in reduced_cost are negative. If they are a solution is found
        for (int i = 0; i < nb; ++i) {
            double temp = reduced_cost(0, i);
            if (temp < mnrc) {
                mnrc = temp;
                mnrc_idx = i;
            }
        }


        // check optimality
        if (mnrc_idx == -1) {
        
            // cheack feasibility
            infeasible = false;
            for (int i = 0; i < m; ++i) {
                int col  = basis_idx[i];

                // Check if this column is an Artificial Variable
                bool is_artificial = false;
                for(int art_idx : artificial_indices) {
                    if(col == art_idx) {
                        is_artificial = true;
                        break;
                    }
                }

                if (is_artificial && x(col,0) > 1e-8) {
                    infeasible = true;
                    break;
                }
            }

            if (infeasible) {
                std::cout << "Problem is infeasible \n";
            } else {
                std::cout << "optimal solution found!\n";
                print_solution();
            }
            return;
        }

        int entering_col = non_basis_idx[mnrc_idx];

        // Get correct column from A
        a_j = mops.get_column(A, entering_col);

        a_j_hat = B_inv * a_j;

        piv_row = -1;
        double min_ratio = 1e20;

        // check witch columns to swap
        bool degenerate = false;
        for (int i = 0; i < m; ++i) {
            double div = a_j_hat(i, 0);
            double x_i = x(basis_idx[i], 0);
            
            // degeneracy check
            /*
            if some value in x_i is zero or smaller the solution on degenerate
            in this case we
            
            */
            if (x_i < 1e-10) {  
                degenerate = true;
                std::cout << "DEGENERACY DETECTED: basic var x" << basis_idx[i] 
                        << " = " << x_i << "\n";
            }

            if (div > 1e-10) {
                double ratio = x_i / div;
                if (ratio < min_ratio) {
                    min_ratio = ratio;
                    piv_row = i;
                }
            }
        }

        if (piv_row == -1) {
            std::cout << "Problem is unbounded \n";
            return; 
        }

        int leaving_col = basis_idx[piv_row];

        // update basis indices
        basis_idx[piv_row] = entering_col;
        auto it = std::find(non_basis_idx.begin(), non_basis_idx.end(), leaving_col);
        if (it != non_basis_idx.end()) {
            *it = entering_col;  // replace the leaving one with the entering one
        }
    }
}

Simplex::Simplex(const matrix& A_in, const matrix& b_in, const matrix& c_trans_in, const std::vector<ConstraintType>& constraint_types_in) {
    // constructor for the simplex

    this->m = A_in.rows;
    this->n = A_in.columns;
    this->constraint_types = constraint_types_in;

    // initialize matrices for optimization
    A = mops.zeros(m,n);
    b = mops.zeros(m,1);
    c_trans = mops.zeros(1,n);

    // check that inputs have valid dimensions
    if (b_in.rows != m) {
        std::cout << "Invalid size of A or b. different amount of rows" << "\n";
    }

    if (b_in.columns != 1) {
        std::cout << "Invalid amount of columns on b. must have one column!" << "\n";
    }

    if (c_trans_in.rows != 1) {
        std::cout << "Invalid amount of rows in c_trans! Must have one row!" << "\n";
    }

    if (c_trans_in.columns != n) {
        std::cout << "Invalid amount of columns in A or c_trans! Must ahve same amount of columns!" << "\n";
    }

    // copy data from input matrices
    // REDO after implementing matrix struct in a smarter way

    // copy A_in data to A
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i,j) = A_in(i,j);
        }
    }
    
    // copy b_in data to b
    for (int i = 0; i < m; ++i) {
        b(i,0) = b_in(i,0);
    }
    
    // copy c_trans_in data to c_trans
    for (int i = 0; i < n; ++i) {
        c_trans(0,i) = c_trans_in(0,i);
    }
    



    if (DEBUG) {
        std::cout << "After copying inputs \n"; 
        std::cout << "A: \n";
        mops.print(A);
        std::cout << "b: \n";
        mops.print(b);
        std::cout << "c_trans: \n";
        mops.print(c_trans);
    }
    // make sure that elements of b are non negative
    for (int i = 0; i < m; ++i) {
        if (b(i,0) < 0.0) {
            b(i,0) *= -1.0;
            for (int j = 0; j < n; ++j) {
                A(i, j) *= -1.0;
            }
            // flip constrain type
            if (constraint_types[i] == LEQ) {
                constraint_types[i] = GEQ;
            } else if (constraint_types[i] == GEQ) {
                constraint_types[i] = LEQ;
            }
        }
    }

    for (int i = 0; i < m; ++i) {
        if (constraint_types[i] == LEQ) {
                n_leq += 1;
            } else if (constraint_types[i] == GEQ) {
                n_geq += 1;
            } else {
                n_eq += 1;
            }
    }

    // take into account if cosntraints are equality, greater than or lesser than
    // if LEQ add only a slack variable
    // if GEQ add slack and artificial big-M
    // if EQ add artificial big-M

    // construct big-M matrices

    n_big = n + n_leq + n_eq + 2 * n_geq; 
    A_big = mops.zeros(m, n_big);
    c_trans_big = mops.zeros(1, n_big);

    // copy A into left side of A_big
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A_big(i,j) = A(i,j);
        }
    }

    for (int j = 0; j < n; ++j) {
        c_trans_big(0, j) = c_trans_in(0, j);
    }

    // add original non absis indices
    for (int j = 0; j < n; ++j) {
        non_basis_idx.push_back(j);
    }

    int current_col = n;
    for (int i = 0; i < m; ++i) {
        if (constraint_types[i] == LEQ) {
            A_big(i, current_col) = 1.0;
            c_trans_big(0, current_col) = 0.0;

            basis_idx.push_back(current_col);

            current_col++;
        } else if (constraint_types[i] == GEQ) {
            A_big(i, current_col) = -1.0;
            
            c_trans_big(0, current_col) = 0.0;

            non_basis_idx.push_back(current_col);

            current_col++;   

            A_big(i, current_col) = 1.0;
            
            c_trans_big(0, current_col) = M;

            basis_idx.push_back(current_col);
            artificial_indices.push_back(current_col);

            current_col++;

        } else {
            A_big(i, current_col) = 1.0;
            c_trans_big(0, current_col) = M;

            basis_idx.push_back(current_col);
            artificial_indices.push_back(current_col);

            current_col++;   
        }
    }

    A = A_big;
    c_trans = c_trans_big;

    if (DEBUG) {
        std::cout << "After slacks and artificial variables \n";
    std::cout << "A: \n";
    mops.print(A);
    std::cout << "c_trans: \n";
    mops.print(c_trans);
    }
    

    nb = n_big - m;

    B = mops.zeros(m, m);
    N = mops.zeros(m, n_big - m);
    x = mops.zeros(n_big, 1);
    c_B = mops.zeros(1, m);
    c_N = mops.zeros(1, n_big - m);
    a_j = mops.zeros(m, 1);
    std::cout << "Simplex initialized" << "\n";

    solve();

    //calculate_current_solution();
    //print_solution();
    
}
