// dual_simplex.hpp
#ifndef DUAL_SIMPLEX_HPP
#define DUAL_SIMPLEX_HPP

#include "matrix.h"
#include "status.hpp"
#include <vector>

constexpr double Infinity = 1.0e30;

enum class VariableStatus : std::uint8_t {
    AtLowerBOund,
    AtUpperBound,
    Fixed,
    Basic
};

class DualSimplex {
    private:
        matrix A;
        matrix b;
        matrix c;
        matrix l; // lower bound
        matrix u; // upper bounds
        matrix A_b;
        matrix x_b;
        matrix c_hat;

        std::vector<int> B; // basic indices
        std::vector<int> N; // non basic indices
        std::vector<VariableStatus> status_;

        Status find_dual_feasible_basis();
        Status calculate_current_solution();

    public:
        
        Status solve();
        Status reset();
        Status return_opt_cost(double& opt_cost);
        Status print_solution();
        Status return_solution(matrix& solution);

        DualSimplex(const matrix& A, const matrix& b, const matrix& c, const matrix& l, const matrix& u);
        ~DualSimplex() = default;
};

#endif // dual_simplex.hpp