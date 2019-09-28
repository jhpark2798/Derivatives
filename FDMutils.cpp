#include "FDMutils.h"

#include <iostream>
#include <memory>

using std::cout;
using std::endl;

void testFDMutils()
{
	std::unique_ptr<Matrix> mat = createMatrix(2, 2);
	(*mat)[0].push_back(1);
	(*mat)[0].push_back(2);
	(*mat)[1].push_back(3);
	(*mat)[1].push_back(4);
	printMatrix(mat);

	std::unique_ptr<Vector> vec = createVector(2);
	vec->push_back(1);
	vec->push_back(2);
	printVector(vec);

	std::unique_ptr<Vector> b = createVector(1);
	b->push_back(1);
	std::unique_ptr<Vector> c = createVector(1);
	c->push_back(1);

	std::unique_ptr<Matrix> tri = createTridiagonalMatrix(vec, b, c);
	printMatrix(tri);
}

// a : diagonal vector
// b,c : subdiagonal vector
std::unique_ptr<Matrix> createTridiagonalMatrix(
	const std::unique_ptr<Vector>& a, const std::unique_ptr<Vector>& b, const std::unique_ptr<Vector>& c)
{
	int size = a->size();
	std::unique_ptr<Matrix> res = createMatrix(size, size);
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			(*res)[i].push_back(0);
	for (int i = 0; i < size; ++i)	(*res)[i][i] = (*a)[i];
	for (int i = 0; i < size - 1; ++i) (*res)[i + 1][i] = (*b)[i];
	for (int i = 0; i < size - 1; ++i) (*res)[i][i + 1] = (*c)[i];
	return res;
}
	
// create vector<vector<double>> with capacity1, capacity2
std::unique_ptr<Matrix> createMatrix(int capacity1, int capacity2)
{
	auto mat = std::make_unique<Matrix>();
	mat->reserve(capacity1);
	for (int i = 0; i < capacity1; ++i) {
		Vector temp;
		temp.reserve(capacity2);
		mat->push_back(temp);
	}
	return mat;
}

// create vector<double> with capacity
std::unique_ptr<Vector> createVector(int capacity)
{
	// Vector* vec = new Vector();
	auto vec = std::make_unique<Vector>();
	vec->reserve(capacity);
	return vec;
}

void printMatrix(const Matrix& mat)
{
	cout << endl;
	for (int i = 0; i < mat.size(); ++i) {
		for (int j = 0; j < mat[i].size(); ++j) {
			cout << mat[i][j] << ", ";
		}
		cout << endl;
	}
}

void printMatrix(const std::unique_ptr<Matrix>& mat)
{
	printMatrix(*mat);
}

void printVector(const Vector& vec)
{
	cout << endl;
	for (int i = 0; i < vec.size(); ++i)
		cout << vec[i] << ", ";
	cout << endl;
}

void printVector(const std::unique_ptr<Vector>& vec)
{
	printVector(*vec);
}