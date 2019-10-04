#pragma once

#include "FDMutils.h"
#include <deque>

// 이거 완성하고, 나중에 클래스 구조를 다 바꾸자.
// pricing method로 class로 만들지 말고 상품 종류로 만들어서
// 그 안에 fdm으로 푸는 함수, mc로 푸는 함수 등등 이렇게 만드는게 나을 것 같음.
// 어차피 K, TTM, r, sigma 이런거는 상품에 종속적인 변수들이니까.
// 일단 과제 해야되니까 ELS 프라이싱 FDM 구현부터 해보자.

// price 함수 만들어야 함.
// solve 함수 안에를 좀 더 다른 함수로 만들어서 깔끔하게 만들 필요가 있음.
class ELStwo {
public:
	ELStwo();
	ELStwo(Vector S0, double sigma1, double sigma2, double rho, double r, double div1, double div2,
		double TTM, Vector cr, Vector K, double KI, int pp, int Nx, int Ny,
		double S1max, double S2max, double S1min, double S2min);

	double price(Vector S);

private:
	// parameters
	Vector S0;
	double sigma1;
	double sigma2;
	double rho;
	double r;
	double div1;
	double div2;

	// payoff of ELS
	double TTM;
	Vector cr;
	Vector K;
	double KI;	// Knock-In barrier level

	// Spacing
	int pp; // number of time points in each 6month
	int Nx;
	int Ny;
	int M;
	double dx;
	double dy;
	double dt;
	double S1max;
	double S2max;
	double S1min;
	double S2min;

	std::unique_ptr<Matrix> solve(bool knockIn);
	std::unique_ptr<Matrix> createCoefMatrix(int N, double sigma);
	std::unique_ptr<Vector> createFirstConstVector(int y, const std::unique_ptr<Matrix>& vLow);
	std::unique_ptr<Vector> createSecondConstVector(int x, const std::unique_ptr<Matrix>& vMid);
	void writeBoundaryCondition(const std::unique_ptr<Matrix>& v);
	std::deque<std::unique_ptr<Matrix>> initPayoff(bool knockIn);
	void updatePayoff(std::unique_ptr<Matrix>& vLow, std::deque<std::unique_ptr<Matrix>>& bc);
	void printDegreeOfProcess(int k, int pp);
};

void testELStwo();