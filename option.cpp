#include "option.h"
#include <cmath>
#include <random>
#include <stdexcept>

void Option::setEvalDate(Date d) {
	evalDate_ = d;
	t_ = daysBetween(evalDate_, expiration_) / DAYS_OF_YEAR;
}

void Option::setProcess(GBMProcess p) {
	p_ = p;
	s_ = p_.getSpot();
	r_ = p_.getRf(expiration_);
	q_ = p_.getDiv(expiration_);
	sigma_ = p_.getVol(expiration_);
}

double Option::getd1() {
	return (log(s_ / strike_) + (r_ - q_ + 0.5 * sigma_ * sigma_) * t_) / (sigma_ * sqrt(t_));
}

double Option::getd2() {
	return getd1() - sigma_*sqrt(t_);
}

double Option::mcprice(int numOfSimulation) {
	double sumOfPayoff = 0;
	double df = exp(-r_ * t_);
	std::mt19937_64 gen;
	std::normal_distribution<double> engine(0.0, 1.0);
	gen.seed(std::random_device{}());
	double es = s_ * exp((r_ - q_ - 0.5 * sigma_ * sigma_) * t_);
	double diffution = sigma_ * sqrt(t_);
	for (unsigned int i = 0; i < numOfSimulation; ++i) {
		double e = engine(gen);
		for (int j = 0; j < 2; ++j) {
			double st = es * exp(diffution * (1 - j * 2) * e);
			double p = (*payoff_)(st);
			sumOfPayoff += df * p;
		}
	}
	return sumOfPayoff / numOfSimulation / 2.0;
}

double Option::bntprice(unsigned int nsteps) {
	double dt = t_ / nsteps;
	double u = exp(sigma_ * sqrt(dt));
	double d = 1 / u;
	double p = (exp((r_ - q_) * dt) - d) / (u - d);
	double df = exp(-r_ * dt);
	std::vector<double> v(nsteps + 1, 0.0);
	for (int j = 0; j <= nsteps; ++j) {
		double st = s_ * pow(u, nsteps - j) * pow(d, j);
		v[j] = (*payoff_)(st);
	}
	for (int i = nsteps - 1; i >= 0; --i) {
		for (int j = 0; j <= i; ++j)
			v[j] = df * (v[j] * p + v[j + 1] * (1 - p));
	}
	return v[0];
}

double Option::impliedVol(double mktPrice, double init, double tol) {
	double origin = sigma_;
	double x = init;
	double e = 1;
	while (e > tol) {
		sigma_ = x;
		double diff = price() - mktPrice;
		e = abs(diff);
		x = x - diff / vega();
	}
	sigma_ = origin;
	return x;
}

double Option::impliedVolBisection(double mktPrice, double init, double tol) {
	double origin = sigma_;
	double x0 = 0.01;
	double x1 = init;
	sigma_ = x0;
	double price0 = price();
	sigma_ = x1;
	double price1 = price();
	if ((price0 - mktPrice) * (price1 - mktPrice) > 0)
		throw std::runtime_error(std::string(__FUNCTION__) + std::string("initial value is wrong"));
	double e = 1;
	while (e > tol) {
		sigma_ = (x0 + x1) / 2;
		price0 = price();
		if (price0 - mktPrice > 0) x1 = sigma_;
		else x0 = sigma_;
		sigma_ = x0;
		e = abs(price() - mktPrice);
	}
	sigma_ = origin;
	return x0;
}

void Option::setSpot(double s) {
	p_.setSpot(s);
	s_ = s;
}

double Option::delta() {
	double origin = s_;
	double ds = 0.01;
	s_ = origin * (1+ds);
	double upperPrice = price();
	s_ = origin * (1-ds);
	double lowerPrice = price();
	double delta = (upperPrice - lowerPrice) / (2 * ds * origin);
	s_ = origin;
	return delta;
}

double Option::gamma() {
	double origin = s_;
	double ds = 0.01;
	double currPrice = price();
	s_ = origin * (1 + ds);
	double upperPrice = price();
	s_ = origin * (1 - ds);
	double lowerPrice = price();
	double gamma = (upperPrice - 2 * currPrice + lowerPrice) / (ds * ds* origin * origin);
	s_ = origin;
	return gamma;
}

double Option::vega() {
	double origin = sigma_;
	double dv = 0.01;
	sigma_ = origin + dv;
	double upperPrice = price();
	sigma_ = origin - dv;
	double lowerPrice = price();
	double vega = (upperPrice - lowerPrice) / (2 * dv);
	sigma_ = origin;
	return vega;
}

double Option::rho() {
	double origin = r_;
	double dr = 0.0001;
	r_ = origin + dr;
	double upperPrice = price();
	r_ = origin - dr;
	double lowerPrice = price();
	double rho = (upperPrice - lowerPrice) / (2 * dr);
	r_ = origin;
	return rho;
}

double Option::theta() {
	double origin = t_;
	double dt = 0.0001;
	double currPrice = price();
	t_ = origin - dt;
	double upperPrice = price();
	double theta = (upperPrice - currPrice) / dt;
	t_ = origin;
	return theta;
}

