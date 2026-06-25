// simplex.h
#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "matrix.h"
#include <vector>

class Simplex {
public:
    struct Solution {
        bool   feasible = false;   
        bool   bounded  = false;   
        double optimum  = 0.0;     
        std::vector<double> x;     
    };

    // Maximize cᵀx  subject to  Ax ≤ b,  x ≥ 0
    Solution solve(const std::vector<double>&              c,
                   const std::vector<std::vector<double>>& A,
                   const std::vector<double>&              b);
};

#endif // SIMPLEX_H