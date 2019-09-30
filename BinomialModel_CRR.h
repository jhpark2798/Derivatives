// Option pricing using CRR Binomial Model

#pragma once

#include <vector>

class BinomialModel_CRR
{
	typedef std::vector<std::vector<double>> vec;

public:
	BinomialModel_CRR(const BinomialModel_CRR& p);
	~BinomialModel_CRR();
	BinomialModel_CRR& operator=(const BinomialModel_CRR& p);

	BinomialModel_CRR(
		double T,
		double S,
		double r,
		double sigma,
		double q,
		int n,
		bool call
	);

	double optionPriceForStrike(double K);
	
protected:
	//double getStockPrice();

private:
	double T;	// expiration time
	double S;	// stock price
	double r;	// interest rate
	double sigma;	// volatility
	double q;	// dividend yield
	int n;	// number of steps
	double call;	// is call?

	void computePriceStep(int i, int j, double K, vec& prices, double p_u, double p_d, double u);

};