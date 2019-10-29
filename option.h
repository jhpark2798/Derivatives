#pragma once

#include "date.h"
#include "GBMProcess.h"
#include "payoff.h"

enum OptionType { Call = 1, Put = -1 };

class Option {
public:
	Option(Date expiration, double strike, OptionType type) :
		expiration_(expiration), strike_(strike), type_(type) {}
	virtual ~Option() { 
		delete payoff_; 
	}

	Date getExpiration() { return expiration_; }
	Date getEvalDate() { return evalDate_; }
	double getStrike() { return strike_; }
	OptionType getType() { return type_; }

	double getSpot() { return s_; }
	double getVolatility() { return sigma_; }
	double getInterestRate() { return r_; }

	void setProcess(GBMProcess p);
	void setEvalDate(Date d);
	void setSpot(double s);
	void setVolatility(double sigma) { sigma_ = sigma; }
	void setInterestRate(double yield) { r_ = yield; }

	virtual double price() = 0;
	virtual double mcprice(int numOfSimulation);
	virtual double bntprice(unsigned int nsteps);
	virtual double delta();
	virtual double gamma();
	virtual double vega();
	virtual double rho();
	virtual double theta();
	virtual double impliedVol(double mktPrice, double init=0.5, double tol=1e-6);
	virtual double impliedVolBisection(double mktPrice, double init=0.5, double tol=1e-6);

protected:
	double getd1();
	double getd2();


	Payoff* payoff_;
	Date evalDate_;
	Date expiration_;
	double strike_;
	OptionType type_;
	GBMProcess p_;
	double s_, r_, q_, sigma_, t_;
};
