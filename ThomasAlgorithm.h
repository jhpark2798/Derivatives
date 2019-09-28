#pragma once

#include <vector>
#include "FDMutils.h"

/*
아직 안전장치가 안만들어졌음. std::vector 받는 생성자에서
A와 d의 사이즈가 일치하는지 확인을 안하고 있고
행렬이 Tridiagonal 인지, Square인지, Strictly diagonally dominant 인지 확인 해야함.
*/

class ThomasAlgorithm {
public:
	ThomasAlgorithm();
	ThomasAlgorithm(const std::unique_ptr<Matrix>& A, const std::unique_ptr<Vector>& d);
	~ThomasAlgorithm();

	int eqtSize();
	std::unique_ptr<Vector> solve();

private:
	std::unique_ptr<Matrix> A;
	std::unique_ptr<Vector> d;

	// 아직 미구현 된 안전장치 함수들
	bool isSquare();
	bool isTridiagonal();
	bool isStrictlyDiagonallyDominant();
};

void testThomasAlgorithm();
