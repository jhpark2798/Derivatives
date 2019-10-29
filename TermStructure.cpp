#include "TermStructure.h"
#include "date.h"

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>

TermStructure::TermStructure(std::vector<Date> dates, std::vector<double> rates)
	: dates_(dates), rates_(rates) {
	// dates가 시간 순으로 정렬되어 있는지 확인
	for (auto iter = dates_.begin() + 1; iter != dates_.end(); ++iter)
		if (*(iter - 1) >= *iter) throw std::runtime_error("termstructure's dates is invalid");
}

double TermStructure::value(Date d) {
	if (d < dates_.front())
		throw std::out_of_range(std::string(__FUNCTION__) + std::string(": invalid Date"));
	auto iter = std::find_if(dates_.begin(), dates_.end(),	[d](Date date) { return date > d; });
	if (iter == dates_.end()) --iter;	// dates_의 모든 원소가 d보다 작은 경우 처리
	int t = std::distance(dates_.begin(), iter);
	double slope = (rates_[t] - rates_[t-1]) / daysBetween(*(iter - 1), *iter);
	return rates_[t-1] + slope * daysBetween(*(iter - 1), d);
}

double YieldTermStructure::discount(Date d) {
	//return exp(-value(d) * daysBetween(dates_.front(), d) / DAYS_OF_YEAR);
	// daily compounding
	return 1 / std::pow(1 + value(d) / DAYS_OF_YEAR, daysBetween(dates_.front(), d));
}

double YieldTermStructure::forwardRate(Date d1, Date d2) {
	return DAYS_OF_YEAR / daysBetween(d1, d2) * (discount(d1) / discount(d2) - 1);
}

double VolatilityTermStructure::variance(Date d) {
	return daysBetween(dates_.front(), d) / DAYS_OF_YEAR * std::pow(value(d), 2);
}


