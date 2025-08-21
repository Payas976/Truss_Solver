#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
private:
    std::vector<std::vector<double>> data;
    int rows, cols;
public:
    Matrix(int r, int c);
    double& operator()(int i, int j);
    Matrix operator*(const Matrix& other);
    std::vector<double> solve(const std::vector<double>& b); // Gaussian elimination
};

#endif