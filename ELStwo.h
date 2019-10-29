#pragma once

#include "FDMutils.h"
#include "date.h"
#include <deque>
#include <utility>

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
	ELStwo(Vector S0, Vector S, double sigma1, double sigma2, double rho, 
		double r, double div1, double div2,
		Date evalDate, std::vector<Date> redemDate,
		double TTM, Vector cr, Vector K, double KI, int pp, int Nx, int Ny,
		double S1max, double S2max, double S1min, double S2min);
	ELStwo(const ELStwo& els);

	void setS1(double S1);
	void setS2(double S2);
	void setVol1(double sigma1);
	void setVol2(double sigma2);
	void setInterestRate(double r);
	void setCorrelation(double rho);
	void setDiv1(double div1);
	void setDiv2(double div2);

	Matrix getPriceMatrix();

	double price(Vector S);
	double delta(Vector S);
	std::vector<std::pair<double, double>> deltaX(Vector S);
	double deltaX();
	std::vector<std::pair<double, double>> deltaY(Vector S);
	double deltaY();
	std::vector<std::pair<double, double>> gammaX(Vector S);
	double gammaX();
	std::vector<std::pair<double, double>> gammaY(Vector S);
	double gammaY();
	// std::vector<std::pair<double, double>> gammaXY(std::vector<std::pair<double, double>> S);
	std::vector<std::pair<double, double>> vegaX(Vector vol);
	std::vector<std::pair<double, double>> vegaY(Vector vol);

private:
	// parameters
	Vector S0;
	Vector S;
	double sigma1;
	double sigma2;
	double rho;
	double r;
	double div1;
	double div2;

	// payoff of ELS
	Date evalDate;
	std::vector<Date> redemDate;
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

	std::shared_ptr<Matrix> knockInMatrix;
	std::shared_ptr<Matrix> noKnockInMatrix;

	void solve();
	void updatePrice();
	std::shared_ptr<Matrix> createCoefMatrix(int N, double sigma, double div);
	std::shared_ptr<Vector> createFirstConstVector(int y, const std::shared_ptr<Matrix>& vLow);
	std::shared_ptr<Vector> createSecondConstVector(int x, const std::shared_ptr<Matrix>& vMid);
	void writeBoundaryCondition(const std::shared_ptr<Matrix>& v);
	std::deque<std::shared_ptr<Matrix>> initPayoff(bool knockIn);
	void updatePayoff(std::shared_ptr<Matrix>& vLow, std::deque<std::shared_ptr<Matrix>>& bc);
	void printDegreeOfProcess(Date currDate);
	void reset();
};

void testELStwo();
void HW();

extern "C" __declspec(dllexport)
double _stdcall kiwoomPrice(int date, double S1, double S2);


