#pragma once

#include <vector>
#include <memory>

typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<double> Vector;

std::unique_ptr<Matrix> createTridiagonalMatrix(
	const std::unique_ptr<Vector>& a,
	const std::unique_ptr<Vector>& b,
	const std::unique_ptr<Vector>& c);
std::unique_ptr<Matrix> createMatrix(int cap1, int cap2);
std::unique_ptr<Vector> createVector(int capacity);
void printMatrix(const Matrix& mat);
void printMatrix(const std::unique_ptr<Matrix>& mat);
void printVector(const Vector& vec);
void printVector(const std::unique_ptr<Vector>& vec);

void testFDMutils();