// simplex.cpp

#include "simplex.h"
#include <limits>
#include <cmath>
#include <iostream>

static Matrix ops;   // we only need one global instance for row operations

static int find_pivot_column(const matrix& tableau, int total_vars) {
    int pivot_col = -1;
    double most_negative = -1e-12;   // tolerance

    for (int j = 0; j < total_vars; ++j) {
        double val = tableau(tableau.rows - 1, j);
        if (val < most_negative) {
            most_negative = val;
            pivot_col = j;
        }
    }
    return pivot_col;
}

static int find_pivot_row(const matrix& tableau, int pivot_col, int m) {
    int pivot_row = -1;
    double min_ratio = std::numeric_limits<double>::infinity();

    for (int i = 0; i < m; ++i) {
        double a = tableau(i, pivot_col);
        if (a <= 1e-12) continue;                // not strictly positive → skip

        double ratio = tableau(i, tableau.columns - 1) / a;
        if (ratio < min_ratio - 1e-12 ||                          // strictly better
            (std::abs(ratio - min_ratio) < 1e-12 && pivot_row != -1 && 
             tableau(i, tableau.columns - 1) < tableau(pivot_row, tableau.columns - 1))) {
            min_ratio = ratio;
            pivot_row = i;
        }
    }
    return pivot_row;
}

static void pivot(matrix& tableau, int pivot_row, int pivot_col) {
    double pivot_elem = tableau(pivot_row, pivot_col);

    // Make pivot = 1
    ops.multiply_row(tableau, pivot_row, 1.0 / pivot_elem);

    // Eliminate column everywhere else (including objective row)
    for (int i = 0; i < tableau.rows; ++i) {
        if (i == pivot_row) continue;
        double factor = tableau(i, pivot_col);
        ops.add_multiple_of_row(tableau, i, pivot_row, -factor);
    }
}


Simplex::Solution Simplex::solve(const std::vector<double>&              c,
                                 const std::vector<std::vector<double>>& A,
                                 const std::vector<double>&              b)
{
    int m = (int)A.size();       // constraints
    int n = (int)c.size();       // decision variables

    if (m == 0 || n == 0 || A[0].size() != n || b.size() != m) {
        return {};
    }

    // Tableau layout:
    // rows:    m  constraint rows + 1 objective row
    // columns: n decision + m slack + 1 RHS
    matrix tableau = ops.zeros(m + 1, n + m + 1);

    // ----- Fill A -----
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            tableau(i, j) = A[i][j];

    // ----- Slack variables (identity) -----
    for (int i = 0; i < m; ++i)
        tableau(i, n + i) = 1.0;

    // ----- RHS (b) -----
    for (int i = 0; i < m; ++i) {
        if (b[i] < 0) {                     // make sure RHS ≥ 0 for initial basis
            for (int j = 0; j < tableau.columns; ++j)
                tableau(i, j) = -tableau(i, j);
        } else {
            tableau(i, tableau.columns - 1) = b[i];
        }
    }

    // ----- Objective row: -c (maximization → move to left side) -----
    for (int j = 0; j < n; ++j)
        tableau(m, j) = -c[j];
    // constant term stays 0

    const int MAX_ITER = 10000;
    int iter = 0;

    while (++iter <= MAX_ITER) {
        int col = find_pivot_column(tableau, n + m);
        if (col == -1) break;                 // optimality reached

        int row = find_pivot_row(tableau, col, m);
        if (row == -1) {
            return { true, false, 0.0, {} };      // unbounded
        }

        pivot(tableau, row, col);
    }

    if (iter > MAX_ITER) {
        std::cerr << "Simplex: too many iterations!\n";
        return {};
    }

    // ----- Extract solution -----
    std::vector<double> x(n, 0.0);
    double objective_value = -tableau(m, tableau.columns - 1);  // because we used -c

    // Identify basic columns
    for (int j = 0; j < n; ++j) {               // only original variables
        int basic_row = -1;
        bool is_basic = true;

        for (int i = 0; i < m; ++i) {
            double val = tableau(i, j);
            if (std::abs(val - 1.0) < 1e-10) {
                if (basic_row != -1) { is_basic = false; break; }
                basic_row = i;
            } else if (std::abs(val) > 1e-10) {
                is_basic = false;
                break;
            }
        }
        if (is_basic && basic_row != -1) {
            x[j] = tableau(basic_row, tableau.columns - 1);
        }
    }

    // Final feasibility check (should always be true here, but be safe)
    bool feasible = true;
    for (int i = 0; i < m; ++i) {
        if (tableau(i, tableau.columns - 1) < -1e-10) {
            feasible = false;
            break;
        }
    }

    return { feasible, true, objective_value, x };
}