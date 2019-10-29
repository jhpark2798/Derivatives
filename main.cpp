#include <iostream>
#include <random>
#include <vector>
#include <iomanip>
#include <numeric>
#include <string>

#include "binary_option.h"
#include "option.h"
#include "plainvanilla_option.h"
#include "date.h"
#include "TermStructure.h"

using std::cout;
using std::endl;

void print_option(Option* x) {
	std::cout << "price = " << x->price() << std::endl;
	std::cout << "mc = " << x->mcprice(1000) << std::endl;
	std::cout << "bnt = " << x->bntprice(100) << std::endl;
	std::cout << std::string(30, '-') << std::endl;
}

double vectorSum(std::vector<double> vec) {
	return std::accumulate(vec.begin(), vec.end(), 0.0);
}

enum Position { Long = 1, Short = -1 };
enum OptionProduct{Vanilla, Binary};

std::ostream& operator<<(std::ostream& os, Position pos) {
	switch (pos) {
	case Long:
		os << "Long";
		break;
	case Short:
		os << "Short";
		break;
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, OptionProduct op) {
	switch (op) {
	case Vanilla:
		os << "Vanilla";
		break;
	case Binary:
		os << "Binary";
		break;
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, OptionType type) {
	switch (type) {
	case Call:
		os << "Call";
		break;
	case Put:
		os << "Put";
		break;
	}
	return os;
}

void printGreeks(std::vector<Option*> inst, std::vector<OptionProduct> product,	
	std::vector<Position> position, std::vector<int> quantity) {
	// 현재 주가에서 가격, 그릭 출력
	std::vector<double> priceVec;
	std::vector<double> deltaVec;
	std::vector<double> gammaVec;
	std::vector<double> vegaVec;
	std::vector<double> rhoVec;
	std::vector<double> thetaVec;
	cout << "All number are calculated considering quantity and position" << endl;
	cout << "Spot : " << inst[0]->getSpot()
		<< "\t Evaluation Date : " << inst[0]->getEvalDate().to_str() << endl << endl;
	for (int i = 0; i < inst.size(); ++i) {
		cout << std::setprecision(5) << std::left << std::setw(8) << product[i]
			<< std::left << std::setw(6) << position[i]
			<< std::right << std::setw(4) << quantity[i] << ' '
			<< std::left << std::setw(6) << inst[i]->getType()
			<< std::right << std::setw(8) << "strike : "
			<< std::left << std::setw(5) << inst[i]->getStrike()
			<< std::right << std::setw(10) << "expiration : "
			<< std::left << std::setw(10) << inst[i]->getExpiration().to_str() << endl;
		double price = inst[i]->price();
		double delta = inst[i]->delta();
		double gamma = inst[i]->gamma();
		double vega = inst[i]->vega();
		double rho = inst[i]->rho();
		double theta = inst[i]->theta();
		priceVec.push_back(price * quantity[i] * position[i]);
		deltaVec.push_back(delta * quantity[i] * position[i]);
		gammaVec.push_back(gamma * quantity[i] * position[i]);
		vegaVec.push_back(vega * quantity[i] * position[i]);
		rhoVec.push_back(rho * quantity[i] * position[i]);
		thetaVec.push_back(theta * quantity[i] * position[i]);
		cout << "price : " << std::left << std::setw(10) << priceVec[i]
			<< std::right << std::setw(10) << "delta : "
			<< std::left << std::setw(10) << deltaVec[i]
			<< std::right << std::setw(10) << "gamma : "
			<< std::left << std::setw(10) << gammaVec[i]
			<< std::right << std::setw(10) << "vega : "
			<< std::left << std::setw(10) << vegaVec[i]
			<< std::right << std::setw(10) << "rho : "
			<< std::left << std::setw(10) << rhoVec[i]
			<< std::right << std::setw(10) << "theta : "
			<< std::left << std::setw(10) << thetaVec[i] << endl << endl;
	}
	cout << "==================== Book Greeks ====================" << endl;
	cout << "price : " << std::left << std::setw(10) << vectorSum(priceVec)
		<< std::right << std::setw(10) << "delta : "
		<< std::left << std::setw(10) << vectorSum(deltaVec)
		<< std::right << std::setw(10) << "gamma : "
		<< std::left << std::setw(10) << vectorSum(gammaVec)
		<< std::right << std::setw(10) << "vega : "
		<< std::left << std::setw(10) << vectorSum(vegaVec)
		<< std::right << std::setw(10) << "rho : "
		<< std::left << std::setw(10) << vectorSum(rhoVec)
		<< std::right << std::setw(10) << "theta : "
		<< std::left << std::setw(10) << vectorSum(thetaVec) << endl << endl;
}

double portfolioPrice(std::vector<Option*> inst, std::vector<OptionProduct> product,	
std::vector<Position> position, std::vector<int> quantity) {
	std::vector<double> priceVec;
	for (int i = 0; i < product.size(); ++i) {
		double price = inst[i]->price();
		priceVec.push_back(price * quantity[i] * position[i]);
	}
	return vectorSum(priceVec);
}

double portfolioDelta(std::vector<Option*> inst, std::vector<OptionProduct> product,
	std::vector<Position> position, std::vector<int> quantity) {
	std::vector<double> deltaVec;
	for (int i = 0; i < product.size(); ++i) {
		double price = inst[i]->delta();
		deltaVec.push_back(price * quantity[i] * position[i]);
	}
	return vectorSum(deltaVec);
}

double portfolioGamma(std::vector<Option*> inst, std::vector<OptionProduct> product,
	std::vector<Position> position, std::vector<int> quantity) {
	std::vector<double> gammaVec;
	for (int i = 0; i < product.size(); ++i) {
		double price = inst[i]->gamma();
		gammaVec.push_back(price * quantity[i] * position[i]);
	}
	return vectorSum(gammaVec);
}

double portfolioVega(std::vector<Option*> inst, std::vector<OptionProduct> product,
	std::vector<Position> position, std::vector<int> quantity) {
	std::vector<double> vegaVec;
	for (int i = 0; i < product.size(); ++i) {
		double price = inst[i]->vega();
		vegaVec.push_back(price * quantity[i] * position[i]);
	}
	return vectorSum(vegaVec);
}

double portfolioRho(std::vector<Option*> inst, std::vector<OptionProduct> product,
	std::vector<Position> position, std::vector<int> quantity) {
	std::vector<double> rhoVec;
	for (int i = 0; i < product.size(); ++i) {
		double price = inst[i]->rho();
		rhoVec.push_back(price * quantity[i] * position[i]);
	}
	return vectorSum(rhoVec);
}

int main() {
	// parameter data
	Date evalDate(2019, 9, 30);
	double spot = 200;
	std::vector<Date> dates{ evalDate, Date(2019,10,30), Date(2019,11,30), Date(2019,12,30),
		Date(2020,1,30), Date(2020,2,28), Date(2020,3,30), Date(2020,4,30) };
	std::vector<double> yieldRate{ 0.015, 0.015, 0.017, 0.0185, 0.0195, 0.0205, 0.0213, 0.022 };
	std::vector<double> divRate{ 0.0, 0.0, 0.0, 0.03, 0.03, 0.03, 0.04, 0.04 };
	std::vector<double> volRate{ 0.1, 0.11, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145 };
	TermStructure divTs(dates, divRate);
	YieldTermStructure yieldTs(dates, yieldRate);
	VolatilityTermStructure volTs(dates, volRate);
	GBMProcess process(spot, divTs, yieldTs, volTs);

	// option data
	std::vector<OptionProduct> product{ Vanilla, Vanilla, Vanilla, Vanilla, Vanilla, Vanilla,
		Binary, Binary, Binary, Binary, Binary, Binary };
	std::vector<Position> position{ Long, Short, Long, Short, Short, Long,
		Short, Long, Short, Long, Long, Short };
	std::vector<int> quantity{ 2, 1, 3, 2, 1, 2, 
		10, 25, 10, 10, 20, 20 };
	std::vector<OptionType> type{ Call, Call, Call, Put, Put, Put, 
		Call, Call, Put, Put, Put, Call };
	std::vector<double> strike{ 200, 205, 195, 200, 210, 190, 
		200, 220, 200, 210, 190, 205 };
	// 2019년 3월 15일 만기 바닐라 풋옵션의 만기를 2020년 3월 15일로 수정
	std::vector<Date> expiration{ 
		Date(2020,1,10), Date(2019,12,12), Date(2020,3,15), Date(2019,12,12), Date(2020,3,15), Date(2020,1,10), 
		Date(2019,11,25), Date(2020,3,20),	Date(2020,2,18), Date(2019,12,19), Date(2020,1,15), Date(2020,2,15) };
	
	// create option vector 
	std::vector<Option*> inst;
	for (int i = 0; i < product.size(); ++i) {
		switch (product[i]) {
		case Vanilla:
			inst.push_back(new PlainVanillaOption(expiration[i], strike[i], type[i]));
			break;
		case Binary:
			inst.push_back(new BinaryOption(expiration[i], strike[i], type[i]));
			break;
		}
		inst[i]->setEvalDate(evalDate);
		inst[i]->setProcess(process);
	}

	// plain vanilla option data
	std::vector<OptionProduct> productPV{ Vanilla, Vanilla, Vanilla, Vanilla, Vanilla, Vanilla };
	std::vector<Position> positionPV{ Long, Short, Long, Short, Short, Long };
	std::vector<int> quantityPV{ 2, 1, 3, 2, 1, 2 };
	std::vector<OptionType> typePV{ Call, Call, Call, Put, Put, Put };
	std::vector<double> strikePV{ 200, 205, 195, 200, 210, 190 };
	// 2019년 3월 15일 만기 바닐라 풋옵션의 만기를 2020년 3월 15일로 수정
	std::vector<Date> expirationPV{ 
		Date(2020,1,10), Date(2019,12,12), Date(2020,3,15), Date(2019,12,12), Date(2020,3,15), Date(2020,1,10) };

	// create option vector 
	std::vector<Option*> instPV;
	for (int i = 0; i < productPV.size(); ++i) {
		switch (productPV[i]) {
		case Vanilla:
			instPV.push_back(new PlainVanillaOption(expiration[i], strike[i], type[i]));
			break;
		case Binary:
			instPV.push_back(new BinaryOption(expiration[i], strike[i], type[i]));
			break;
		}
		instPV[i]->setEvalDate(evalDate);
		instPV[i]->setProcess(process);
	}

	// 주가에 따른 Greeks 출력
	printGreeks(inst, product, position, quantity);
	for (auto iter = inst.begin(); iter != inst.end(); ++iter) (*iter)->setSpot(spot * 1.1);
	printGreeks(inst, product, position, quantity);
	for (auto iter = inst.begin(); iter != inst.end(); ++iter) (*iter)->setSpot(spot * 0.9); 
	printGreeks(inst, product, position, quantity);
	for (auto iter = inst.begin(); iter != inst.end(); ++iter) (*iter)->setSpot(spot);
	
	// 시나리오 별 Plain Vanilla Option 포트폴리오 가격 변동 출력
	cout << "============ Plain vanilla option portfolio price change with various scenarios ============" << endl;
	double ds = spot * 0.01;
	for (auto iter = instPV.begin(); iter != instPV.end(); ++iter) (*iter)->setSpot(spot + ds);
	double upperPrice = portfolioPrice(instPV, productPV, positionPV, quantityPV);
	for (auto iter = instPV.begin(); iter != instPV.end(); ++iter) (*iter)->setSpot(spot - ds);
	double lowerPrice = portfolioPrice(instPV, productPV, positionPV, quantityPV);
	for (auto iter = instPV.begin(); iter != instPV.end(); ++iter) (*iter)->setSpot(spot);
	double currPrice = portfolioPrice(instPV, productPV, positionPV, quantityPV);
	double delta = portfolioDelta(instPV, productPV, positionPV, quantityPV);
	double gamma = portfolioGamma(instPV, productPV, positionPV, quantityPV);
	cout << std::left << std::setw(25) << "Scenario"
		<< std::left << std::setw(15) << "Price Change"
		<< std::left << std::setw(25) << "Price Change by Greeks"
		<< std::left << std::setw(10) << "error" << endl;
	cout << std::left << std::setw(25) << "stock + 1%"
		<< std::left << std::setw(15) << upperPrice - currPrice
		<< std::left << std::setw(25) << delta * ds + 0.5 * gamma * ds * ds
		<< std::left << std::setw(10) << upperPrice - currPrice - (delta * ds + 0.5 * gamma * ds * ds) << endl;
	cout << std::left << std::setw(25) << "stock - 1%"
		<< std::left << std::setw(15) << lowerPrice - currPrice
		<< std::left << std::setw(25) << delta * -ds + 0.5 * gamma * -ds * -ds
		<< std::left << std::setw(10) << lowerPrice - currPrice - (delta * -ds + 0.5 * gamma * -ds * -ds) << endl;

	//변동성에 따른 가격 변동 출력
	double dv = 0.01;
	std::vector<double> vol;
	for (auto iter = instPV.begin(); iter != instPV.end(); ++iter) vol.push_back((*iter)->getVolatility());
	for (int i = 0; i < instPV.size(); ++i) instPV[i]->setVolatility(vol[i] + dv);
	upperPrice = portfolioPrice(instPV, productPV, positionPV, quantityPV);
	for (int i = 0; i < instPV.size(); ++i) instPV[i]->setVolatility(vol[i] - dv);
	lowerPrice = portfolioPrice(instPV, productPV, positionPV, quantityPV);
	for (int i = 0; i < instPV.size(); ++i) instPV[i]->setVolatility(vol[i]);
	currPrice = portfolioPrice(instPV, productPV, positionPV, quantityPV);
	double vega = portfolioVega(instPV, productPV, positionPV, quantityPV);
	cout << std::left << std::setw(25) << "vol + 1%"
		<< std::left << std::setw(15) << upperPrice - currPrice
		<< std::left << std::setw(25) << vega * dv
		<< std::left << std::setw(10) << upperPrice - currPrice - (vega * dv) << endl;
	cout << std::left << std::setw(25) << "vol - 1%"
		<< std::left << std::setw(15) << lowerPrice - currPrice
		<< std::left << std::setw(25) << vega * -dv
		<< std::left << std::setw(10) << lowerPrice - currPrice - (vega * -dv) << endl;

	// 금리에 따른 가격 변동 출력
	double dr = 0.001;
	std::vector<double> yield;
	for (auto iter = instPV.begin(); iter != instPV.end(); ++iter) yield.push_back((*iter)->getInterestRate());
	for (int i = 0; i < instPV.size(); ++i) instPV[i]->setInterestRate(yield[i] + dr);
	upperPrice = portfolioPrice(instPV, productPV, positionPV, quantityPV);
	for (int i = 0; i < instPV.size(); ++i) instPV[i]->setInterestRate(yield[i] - dr);
	lowerPrice = portfolioPrice(instPV, productPV, positionPV, quantityPV);
	for (int i = 0; i < instPV.size(); ++i) instPV[i]->setInterestRate(yield[i]);
	currPrice = portfolioPrice(instPV, productPV, positionPV, quantityPV);
	double rho = portfolioRho(instPV, productPV, positionPV, quantityPV);
	cout << std::left << std::setw(25) << "interest rate + 10bp"
		<< std::left << std::setw(15) << upperPrice - currPrice
		<< std::left << std::setw(25) << rho * dr
		<< std::left << std::setw(10) << upperPrice - currPrice - rho * dr << endl;
	cout << std::left << std::setw(25) << "interest rate - 10bp"
		<< std::left << std::setw(15) << lowerPrice - currPrice
		<< std::left << std::setw(25) << rho * -dr
		<< std::left << std::setw(10) << lowerPrice - currPrice - rho * -dr << endl;
	cout << endl;

	// Plain Vanilla 옵션 내재변동성 출력
	cout << "============ Plain Vanilla Option Implied Volatility ============" << endl;
	std::vector<double> mktPrice{ 5,3,10,5,15,2 };
	cout << "Spot : " << inst[0]->getSpot()
		<< "\t Evaluation Date : " << inst[0]->getEvalDate().to_str() << endl << endl;
	for (int i = 0; i < 6; ++i) {
		cout << std::setprecision(5) << std::left << std::setw(8) << product[i]
			<< std::left << std::setw(10) << type[i]
			<< std::right << std::setw(8) << "strike : "
			<< std::left << std::setw(5) << inst[i]->getStrike()
			<< std::right << std::setw(10) << "expiration : "
			<< std::left << std::setw(10) << inst[i]->getExpiration().to_str()
			<< std::right << std::setw(10) << "mkt Price : "
			<< std::left << std::setw(5) << mktPrice[i] << endl
			<< std::left << std::setw(25) << "Implied Vol by Newton-Raphson : "
			<< std::left << std::setw(10) << inst[i]->impliedVol(mktPrice[i])
			<< std::left << std::setw(25) << "Implied Vol by Bisection : "
			<< std::left << std::setw(10) << inst[i]->impliedVolBisection(mktPrice[i])
			<< endl << endl;
	}

	for (int i = 0; i < inst.size(); ++i)
		delete inst[i];
	for (int i = 0; i < instPV.size(); ++i)
		delete instPV[i];

	return 0;
}

