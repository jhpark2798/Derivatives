#include "FDMutils.h"

#include <iostream>
#include <memory>

using std::cout;
using std::endl;

void testFDMutils()
{
	std::shared_ptr<Matrix> mat = createMatrix(2, 2);
	printMatrix(mat);
	std::shared_ptr<Vector> vec = createVector(2);
	(*vec)[0] = 1; (*vec)[1] = 2;
	printVector(vec);
	std::shared_ptr<Vector> b = createVector(1);
	std::shared_ptr<Vector> c = createVector(1);
	(*b)[0] = 1; (*c)[0] = 1;
	std::shared_ptr<Matrix> tri = createTridiagonalMatrix(vec, b, c);
	printMatrix(tri);
}

// a : diagonal vector
// b,c : subdiagonal vector
std::shared_ptr<Matrix> createTridiagonalMatrix(
	const std::shared_ptr<Vector>& a, const std::shared_ptr<Vector>& b, const std::shared_ptr<Vector>& c)
{
	int size = a->size();
	std::shared_ptr<Matrix> res = createMatrix(size, size);
	for (int i = 0; i < size; ++i)	(*res)[i][i] = (*a)[i];
	for (int i = 0; i < size - 1; ++i) (*res)[i + 1][i] = (*b)[i];
	for (int i = 0; i < size - 1; ++i) (*res)[i][i + 1] = (*c)[i];
	return res;
}
	
// create vector<vector<double>> with length1, length2
std::shared_ptr<Matrix> createMatrix(int length1, int length2)
{
	auto mat = std::make_shared<Matrix>();
	mat->reserve(length1);
	for (int i = 0; i < length1; ++i) {
		std::shared_ptr<Vector> tempVec = createVector(length2);
		mat->push_back(*tempVec);	
	}
	return mat;
}

// create vector<double> with length
std::shared_ptr<Vector> createVector(int length) {
	auto vec = std::make_shared<Vector>(length, 0.0);
	return vec;
}

void printMatrix(const Matrix& mat) {
	cout << endl;
	for (int i = 0; i < mat.size(); ++i) {
		for (int j = 0; j < mat[i].size(); ++j) {
			cout << mat[i][j] << ", ";
		}
		cout << endl;
	}
}

void printMatrix(const std::shared_ptr<Matrix>& mat) {
	printMatrix(*mat);
}

void printVector(const Vector& vec) {
	cout << endl;
	for (int i = 0; i < vec.size(); ++i)
		cout << vec[i] << ", ";
	cout << endl;
}

void printVector(const std::shared_ptr<Vector>& vec) {
	printVector(*vec);
}