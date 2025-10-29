// matrix.cpp
#include "matrix.h"



matrix Matrix::zeros(int r, int c) {
    matrix m;
    m.rows = r;
    m.columns = c;
    m.data = new double[m.rows * m.columns];

    for (int i = 0; i < m.rows * m.columns; i++) {
        m.data[i] = 0.0;
    }

    return m;
};

void Matrix::print(matrix A)  {
    for (int i = 0; i < A.rows; i++)  {
        for (int j = 0; j < A.columns; j++) {
            std::cout << A(i,j) << " ";
        }
        std::cout << "\n";
    }
}

matrix Matrix::multiply(matrix A, matrix B) {
    matrix C = zeros(A.rows, B.columns);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < B.columns; j++)  {
            for (int k = 0; k < A.columns; k++) {
                C(i,j) += A(i,k) * B(k, j);
            }
        }
    }
    return C;
}


matrix Matrix::addition(matrix A, matrix B) {
    matrix C = zeros(A.rows, A.columns);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.columns; j++) {
            C(i,j) = A(i,j) + B(i,j);
        }
    }
    return C;
}
