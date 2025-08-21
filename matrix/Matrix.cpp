#include "Matrix.h"
#include <stdexcept>
#include <cmath>

// Constructor: Initialize matrix with r rows and c columns, filled with zeros
Matrix::Matrix(int r, int c) : rows(r), cols(c) {
    if(r <= 0 || c <= 0) {
        throw std::invalid_argument("Matrix dimensions must be positive");
    }
    data.resize(r, std::vector<double>(c, 0.0));
}

// Element access operator: Access element at (i, j)
double& Matrix::operator()(int i, int j) {
    if(i < 0 || i >= rows || j < 0 || j >= cols) {
        throw std::out_of_range("Matrix index out of bounds");
    }
    return data[i][j];
}

// Matrix multiplication operator
Matrix Matrix::operator*(const Matrix& other) {
    if(cols != other.rows) {
        throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
    }
    Matrix result(rows, other.cols);
    for(int i = 0; i < rows; ++i) {
        for(int j = 0; j < other.cols; ++j) {
            for(int k = 0; k < cols; ++k) {
                result(i, j) += data[i][k] * other.data[k][j];
            }
        }
    }
    return result;
}

// Gaussian elimination to solve Ax = b
std::vector<double> Matrix::solve(const std::vector<double>& b) {
    if(rows != cols) {
        throw std::invalid_argument("Matrix must be square for Gaussian elimination");
    }
    if(b.size() != rows) {
        throw std::invalid_argument("Vector b size must match matrix rows");
    }

    // Create augmented matrix [A|b]
    std::vector<std::vector<double>> aug(rows, std::vector<double>(cols + 1));
    for(int i = 0; i < rows; ++i) {
        for(int j = 0; j < cols; ++j) {
            aug[i][j] = data[i][j];
        }
        aug[i][cols] = b[i];
    }

    // Forward elimination with partial pivoting
    for(int i = 0; i < rows; ++i) {
        // Find pivot
        int maxRow = i;
        double maxVal = std::abs(aug[i][i]);
        for(int k = i + 1; k < rows; ++k) {
            if(std::abs(aug[k][i]) > maxVal) {
                maxVal = std::abs(aug[k][i]);
                maxRow = k;
            }
        }
        if(std::abs(maxVal) < 1e-10) {
            throw std::runtime_error("Matrix is singular or nearly singular");
        }

        // Swap rows
        if(maxRow != i) {
            std::swap(aug[i], aug[maxRow]);
        }

        // Eliminate column
        for(int k = i + 1; k < rows; ++k) {
            double factor = aug[k][i] / aug[i][i];
            for(int j = i; j <= cols; ++j) {
                aug[k][j] -= factor * aug[i][j];
            }
        }
    }

    // Back substitution
    std::vector<double> x(rows);
    for(int i = rows - 1; i >= 0; --i) {
        double sum = 0.0;
        for(int j = i + 1; j < cols; ++j) {
            sum += aug[i][j] * x[j];
        }
        if(std::abs(aug[i][i]) < 1e-10) {
            throw std::runtime_error("Matrix is singular or nearly singular");
        }
        x[i] = (aug[i][cols] - sum) / aug[i][i];
    }

    return x;
}