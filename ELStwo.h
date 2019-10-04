#pragma once

#include "FDMutils.h"
#include <deque>

// �̰� �ϼ��ϰ�, ���߿� Ŭ���� ������ �� �ٲ���.
// pricing method�� class�� ������ ���� ��ǰ ������ ����
// �� �ȿ� fdm���� Ǫ�� �Լ�, mc�� Ǫ�� �Լ� ��� �̷��� ����°� ���� �� ����.
// ������ K, TTM, r, sigma �̷��Ŵ� ��ǰ�� �������� �������̴ϱ�.
// �ϴ� ���� �ؾߵǴϱ� ELS �����̽� FDM �������� �غ���.

// price �Լ� ������ ��.
// solve �Լ� �ȿ��� �� �� �ٸ� �Լ��� ���� ����ϰ� ���� �ʿ䰡 ����.
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