#pragma once

#include <vector>
#include <memory>

typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<double> Vector;

std::shared_ptr<Matrix> createTridiagonalMatrix(
	const std::shared_ptr<Vector>& a,
	const std::shared_ptr<Vector>& b,
	const std::shared_ptr<Vector>& c);
std::shared_ptr<Matrix> createMatrix(int length1, int length2);
std::shared_ptr<Vector> createVector(int length);
void printMatrix(const Matrix& mat);
void printMatrix(const std::shared_ptr<Matrix>& mat);
void printVector(const Vector& vec);
void printVector(const std::shared_ptr<Vector>& vec);

void testFDMutils();