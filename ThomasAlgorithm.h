#pragma once

#include <vector>
#include "FDMutils.h"

/*
���� ������ġ�� �ȸ��������. std::vector �޴� �����ڿ���
A�� d�� ����� ��ġ�ϴ��� Ȯ���� ���ϰ� �ְ�
����� Tridiagonal ����, Square����, Strictly diagonally dominant ���� Ȯ�� �ؾ���.
*/

class ThomasAlgorithm {
public:
	ThomasAlgorithm();
	ThomasAlgorithm(const std::shared_ptr<Matrix>& A, const std::shared_ptr<Vector>& d);
	~ThomasAlgorithm();

	int eqtSize();
	std::shared_ptr<Vector> solve();

private:
	std::shared_ptr<Matrix> A;
	std::shared_ptr<Vector> d;

	// ���� �̱��� �� ������ġ �Լ���
	bool isSquare();
	bool isTridiagonal();
	bool isStrictlyDiagonallyDominant();
};

void testThomasAlgorithm();
