#pragma once
#include <string>
#include <vector>

class Date {
public:
	Date() {}
	Date(int y, int m, int d) 
		: y_(y), m_(m), d_(d), monthsThirtyOneDays({ 1,3,5,7,8,10,12 }) {	}
	Date(std::string ymd);

	int year() { return y_; }
	int month() { return m_; }
	int day() { return d_; }
	int daysFrom(Date d);
	Date addDays(int days);

	void print();
	std::string to_str();
	bool isLeapYear();

	Date& operator++();
	Date operator++(int);
	Date operator+(int days);
	bool operator>(Date& rhs) const;
	bool operator>(const Date& rhs) const;
	bool operator<(Date& rhs) const;
	bool operator<(const Date& rhs) const;
	bool operator>=(Date& rhs) const;
	bool operator>=(const Date& rhs) const;
	bool operator<=(Date& rhs) const;
	bool operator<=(const Date& rhs) const;
	bool operator==(const Date& rhs) const;
	bool operator==(Date& rhs) const;
	bool operator!=(const Date& rhs) const;
	bool operator!=(Date& rhs) const;

private:
	int y_, m_, d_;
	std::vector<int> monthsThirtyOneDays;
};

int daysBetween(Date d1, Date d2);