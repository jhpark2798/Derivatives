#include "BinomialModel_CRR.h"

#include <vector>
#include <cmath>
#include <algorithm>

using std::vector;
using std::exp;
using std::sqrt;
using std::pow;
using std::max;

BinomialModel_CRR::BinomialModel_CRR(const BinomialModel_CRR& p)
	: T(p.T), S(p.S), r(p.r), sigma(p.sigma), q(p.q), n(p.n), call(p.call)
{
}

BinomialModel_CRR::~BinomialModel_CRR()
{
}

BinomialModel_CRR& BinomialModel_CRR::operator=(const BinomialModel_CRR& p)
{
	if (this != &p) {
		T = p.T;
		S = p.S;
		r = p.r;
		sigma = p.sigma;
		q = p.q;
		n = p.n;
		call = p.call;
	}
	return *this;
}

BinomialModel_CRR::BinomialModel_CRR(double T, double S, double r, double sigma, double q, int n, bool call)
	: T(T), S(S), r(r), sigma(sigma), q(q), n(n), call(call)
{
}

void BinomialModel_CRR::computePriceStep(int i, int j, double K, vec& prices, double p_u, double p_d, double u)
{
	prices[i][j] = prices[i + 1][j] * p_u + prices[i + 1][j + 1] * p_d;
}

double BinomialModel_CRR::optionPriceForStrike(double K)
{
	// use Cox, Ross, Rubinstein Binomial Tree
	
	double delta = T / n; // size of each step
	double u = exp(sigma * sqrt(delta));
	double d = exp(-sigma * sqrt(delta));

	double p_u = (exp((r - q)*delta) - d) / (u - d);
	double p_d = 1 - p_u;

	vec prices;
	for (int i = 0; i < n+1; ++i) {
		vector<double> v(n+1);
		prices.push_back(v);
	}

	for (int j = 0; j < n+1; ++j) {
		double stockPrice = S * pow(u, n - j) * pow(d, j);
		if (call)
			prices[n][j] = max(0.0, stockPrice - K);
		else
			prices[n][j] = max(0.0, K - stockPrice);
	}

	for (int i = n - 1; i >= 0; --i) {
		for (int j = 0; j <= i; ++j)	
			computePriceStep(i, j, K, prices, p_u, p_d, u);
	}

	return prices[0][0];
}

// console example

//#include <iostream>
//using std::cout;
//using std::endl;
//
//int main()
//{
//	BinomialModel_CRR bm(1, 50, 0.05, 0.2, 0, 2, true);
//	cout << bm.optionPriceForStrike(50) << endl;
//
//	return 0;
//}