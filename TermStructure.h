#pragma once

#include <vector>

#include "date.h"

#define DAYS_OF_YEAR 365.0

class TermStructure {
public:
	TermStructure() {}
	~TermStructure() {}
	TermStructure(std::vector<Date> dates, std::vector<double> rates);
	double value(Date d);
protected:
	std::vector<Date> dates_;
	std::vector<double> rates_;
};

class YieldTermStructure : public TermStructure {
public:
	YieldTermStructure() {}
	~YieldTermStructure() {}
	YieldTermStructure(std::vector<Date> dates, std::vector<double> rates)
		: TermStructure(dates, rates) {}
	double discount(Date d);
	double forwardRate(Date d1, Date d2);
};

class VolatilityTermStructure : public TermStructure {
public:
	VolatilityTermStructure() {}
	~VolatilityTermStructure() {}
	VolatilityTermStructure(std::vector<Date> dates, std::vector<double> rates)
		: TermStructure(dates, rates) {}
	double variance(Date d);
};