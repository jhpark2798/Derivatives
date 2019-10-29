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
	std::shared_ptr<Matrix> A = std::make_unique<Matrix>(mat);
	std::shared_ptr<Vector> d = std::make_unique<Vector>(vec);
	cout << "Matrix A is " << endl;
	printMatrix(A);
	cout << "Vector d is " << endl;
	printVector(d);

	ThomasAlgorithm ta(A, d);
	std::shared_ptr<Vector> x = ta.solve();
	cout << "solution vector x is " << endl;
	printVector(x);
}

ThomasAlgorithm::ThomasAlgorithm() {}

ThomasAlgorithm::ThomasAlgorithm(const std::shared_ptr<Matrix>& A, const std::shared_ptr<Vector>& d)
	: A(A), d(d) {}

ThomasAlgorithm::~ThomasAlgorithm() {}

int ThomasAlgorithm::eqtSize() {
	return A->size();
}

// Core Function
std::shared_ptr<Vector> ThomasAlgorithm::solve() {
	int size = eqtSize();
	std::shared_ptr<Vector> a_ = createVector(size);	// a'
	std::shared_ptr<Vector> d_ = createVector(size);	// d'
	std::shared_ptr<Vector> x = createVector(size);		// x (solution)

	// calculate a' and d'
	(*a_)[0] = (*A)[0][0];
	(*d_)[0] = (*d)[0];
	for (int i = 1; i < size; ++i) {
		(*a_)[i] = (*A)[i][i] - (*A)[i][i - 1] * (*A)[i - 1][i] / (*a_)[i - 1];
		(*d_)[i] = (*d)[i] - (*A)[i][i - 1] * (*d_)[i - 1] / (*a_)[i - 1];
	}
	// calculate x
	(*x)[size - 1] = (*d_)[size - 1] / (*a_)[size - 1];
	for (int i = size - 2; i >= 0; --i)
		(*x)[i] = (*d_)[i] / (*a_)[i] - (*A)[i][i + 1] / (*a_)[i] * (*x)[i + 1];
	return x;
}

bool ThomasAlgorithm::isSquare()
{
	// �̿ϼ�
	return true;
}

bool ThomasAlgorithm::isTridiagonal()
{
	// �̿ϼ�
	return true;
}

bool ThomasAlgorithm::isStrictlyDiagonallyDominant()
{
	// �̿ϼ�
	return true;
}

