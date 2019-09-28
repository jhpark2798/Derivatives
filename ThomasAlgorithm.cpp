#include "ThomasAlgorithm.h"
#include "FDMutils.h"

#include <vector>
#include <iostream>

using std::cout;
using std::endl;

void testThomasAlgorithm()
{
	Matrix mat;
	mat.push_back({ 3, -1, 0, 0, 0 });
	mat.push_back({ 1, -4, 2, 0, 0 });
	mat.push_back({ 0, -1, 4, 2, 0 });
	mat.push_back({ 0, 0, 3, 5, -1 });
	mat.push_back({ 0, 0, 0, 1, 2 });
	Vector vec =
	{
		1, 2, 3, 4, 5,
	};
	std::unique_ptr<Matrix> A = std::make_unique<Matrix>(mat);
	std::unique_ptr<Vector> d = std::make_unique<Vector>(vec);
	cout << "Matrix A is " << endl;
	printMatrix(A);
	cout << "Vector d is " << endl;
	printVector(d);

	ThomasAlgorithm ta(A, d);
	std::unique_ptr<Vector> x = ta.solve();
	cout << "solution vector x is " << endl;
	printVector(x);
}

ThomasAlgorithm::ThomasAlgorithm()
{
}

ThomasAlgorithm::ThomasAlgorithm(const std::unique_ptr<Matrix>& A, const std::unique_ptr<Vector>& d)
	: A(std::make_unique<Matrix>(*(std::move(A)))), d(std::make_unique<Vector>(*(std::move(d))))
{
}


ThomasAlgorithm::~ThomasAlgorithm()
{
}

int ThomasAlgorithm::eqtSize()
{
	return A->size();
}

// Core Function
std::unique_ptr<Vector> ThomasAlgorithm::solve()
{
	int size = eqtSize();
	std::unique_ptr<Vector> a_ = createVector(size);	// a'
	std::unique_ptr<Vector> d_ = createVector(size);	// d'
	std::unique_ptr<Vector> x = createVector(size);		// x (solution)

	// calculate a' and d'
	(*a_).push_back((*A)[0][0]);
	(*d_).push_back((*d)[0]);
	for (int i = 1; i < size; ++i) {
		(*a_).push_back((*A)[i][i] - (*A)[i][i - 1] * (*A)[i - 1][i] / (*a_)[i - 1]);
		(*d_).push_back((*d)[i] - (*A)[i][i - 1] * (*d_)[i - 1] / (*a_)[i - 1]);
	}

	// calculate x
	for (int i = 0; i < size; ++i) x->push_back(0);
	(*x)[size - 1] = (*d_)[size - 1] / (*a_)[size - 1];
	for (int i = size - 2; i >= 0; --i) {
		(*x)[i] = (*d_)[i] / (*a_)[i] - (*A)[i][i + 1] / (*a_)[i] * (*x)[i + 1];
	}

	return x;
}

bool ThomasAlgorithm::isSquare()
{
	// 미완성
	return true;
}

bool ThomasAlgorithm::isTridiagonal()
{
	// 미완성
	return true;
}

bool ThomasAlgorithm::isStrictlyDiagonallyDominant()
{
	// 미완성
	return true;
}

