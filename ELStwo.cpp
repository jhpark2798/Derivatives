#include "ThomasAlgorithm.h"
#include "FDMutils.h"
#include "ELStwo.h"
#include "date.h"

#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))

#include <time.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <memory>
#include <deque>
#include <algorithm>
#include <stdexcept>

using std::cout;
using std::endl;

void HW() {
	// KOSPI200, EuroStoxx50
	Vector S0{ 244.39, 2997.55 };
	Vector S{ 244.39, 2997.55 };	
	double sigma1 = 0.133588, sigma2 = 0.249192;
	double rho = 0.354973;
	double r = 0.01703;	// È¸»çÃ¤ ¹ÎÆòÆò±Õ
	double div1 = 0.0144, div2 = 0.0316;
	double TTM = 3;
	Date evalDate(2016, 5, 30);
	std::vector<Date> redemDate{
		Date(2016,11,30), Date(2017,6,2), Date(2017,11,29), Date(2018,5,30), Date(2018,11,29), Date(2019,5,30) };
	Vector cr{ 0.037, 0.074, 0.111, 0.148, 0.185, 0.222 };
	Vector K{ 0.90, 0.90, 0.85, 0.85, 0.80, 0.80 };
	double KI = 0.55;
	int pp = 100;
	int Nx =300, Ny = 300;
	double S1max = S0[0] * 3, S1min = 0;
	double S2max = S0[1] * 3, S2min = 0;

	ELStwo els(S0, S, sigma1, sigma2, rho, r, div1, div2,
		evalDate, redemDate,
		TTM, cr, K, KI, pp, Nx, Ny,
		S1max, S2max, S1min, S2min);
	cout << els.price(S) << endl;
	//els.delta(S0);

	// cega
	//els.setCorrelation(rho);
	//double currentPrice = els.price(S0);
	//els.setCorrelation(rho + 0.0001);
	//double upperPrice = els.price(S0);
	//cout << "cega = " << (upperPrice - currentPrice) / 0.0001;
	
	// current price
	//Matrix priceMatrix = els.getPriceMatrix();

	 // price with various scenarios
	 //1. volX
	//std::vector<std::pair<double,double>> p;
	//double volX = 0; 
	//for (int i = 0; i <= 20; ++i) {
	//	volX = 0.05 + i * 0.01;
	//	els.setVol1(volX);
	//	p.push_back(std::pair<double,double>(volX, els.price(S0)));
	//}
	//cout << "volX, price" << endl;
	//for (auto iter = p.begin(); iter != p.end(); ++iter)
	//	cout << iter->first << ", " << iter->second << endl;
	//cout << endl;

	// 2. volY
	//p.clear();
	//double volY = 0;
	//for (int i = 0; i <= 30; ++i) {
	//	volY = 0.10 + i * 0.01;
	//	els.setVol2(volY);
	//	p.push_back(std::pair<double, double>(volY, els.price(S0)));
	//}
	//cout << "volY, price" << endl;
	//for (auto iter = p.begin(); iter != p.end(); ++iter)
	//	cout << iter->first << ", " << iter->second << endl;
	//cout << endl;

	// 3. interest rate
	//p.clear();
	//double rf = 0;
	//for (int i = 0; i <= 20; ++i) {
	//	rf = 0.005 + i * 0.001;
	//	els.setInterestRate(rf);
	//	p.push_back(std::pair<double, double>(rf, els.price(S0)));
	//}
	//cout << "interest rate, price" << endl;
	//for (auto iter = p.begin(); iter != p.end(); ++iter)
	//	cout << iter->first << ", " << iter->second << endl;
	//cout << endl;

	// 4. correlation rate
	//p.clear();
	//double corr = 0;
	//for (int i = 0; i <= 20; ++i) {
	//	corr = -0.305 + i * 0.05;
	//	els.setCorrelation(corr);
	//	p.push_back(std::pair<double, double>(corr, els.price(S0)));
	//}
	//cout << "correlation, price" << endl;
	//for (auto iter = p.begin(); iter != p.end(); ++iter)
	//	cout << iter->first << ", " << iter->second << endl;

	// 5. dividend rate
	//p.clear();
	//for (int i = 0; i < 10; ++i) {
	//	double dividend1 = 0 + i * 0.006;
	//	els.setDiv1(dividend1);
	//	p.push_back(std::pair<double, double>(dividend1, els.price(S0)));
	//}
	//cout << "div1, price" << endl;
	//for (auto iter = p.begin(); iter != p.end(); ++iter)
	//	cout << iter->first << ", " << iter->second << endl;

	//p.clear();
	//for (int i = 0; i < 10; ++i) {
	//	double dividend2 = 0 + i * 0.006;
	//	els.setDiv2(dividend2);
	//	p.push_back(std::pair<double, double>(dividend2, els.price(S0)));
	//}
	//cout << "div2, price" << endl;
	//for (auto iter = p.begin(); iter != p.end(); ++iter)
	//	cout << iter->first << ", " << iter->second << endl;

	// current price
	//cout << "price Matrix is " << endl;
	//printMatrix(priceMatrix);
	//cout << endl;
}

void testELStwo()
{
	Vector S0{ 3500, 3000 }; // eurostoxx and sp500
	Vector S{ 3500, 3000 };
	double sigma1 = 0.25, sigma2 = 0.2;
	double rho = 0.2;
	double r = 0.02;
	double div1 = 0, div2 = 0;
	double TTM = 3;
	Vector cr{ 0.0285, 0.057, 0.0855, 0.114, 0.1425, 0.171 };
	Vector K{ 0.85, 0.85, 0.85, 0.85, 0.85, 0.85 };
	Date evalDate(2016, 1, 1);
	std::vector<Date> redemDate{
		Date(2016,7,1), Date(2017,1,1), Date(2017,7,1), Date(2018,1,1), Date(2018,7,1), Date(2019,1,1) };
	double KI = 0.55;
	int pp = 100;
	int Nx = 300, Ny = 300;
	double S1max = S0[0] *3, S1min = 0;
	double S2max = S0[1] *3, S2min = 0;

	ELStwo els(S0, S, sigma1, sigma2, rho, r, div1, div2, 
		evalDate, redemDate,
		TTM, cr, K, KI, pp, Nx, Ny,
		S1max, S2max, S1min, S2min);

	clock_t start, end;
	start = clock();
	cout << "price = " << els.price(S) << endl;
	end = clock();
	cout << "time = " << double(end - start) << "ms" << endl;

	// print deltaX
	//Vector S;
	//for (int i = 0; i <= 450; ++i) S.push_back(1000 + i * 10);
	//std::vector<std::pair<double, double>> delta = els.deltaX(S);
	//cout << "deltaX" << endl << "Sx price, delta" << endl;
	//for (auto iter = delta.begin(); iter != delta.end(); ++iter) {
	//	cout << iter->first << ", " << iter->second << endl;
	//}

	// print deltaY
	//S.clear(); delta.clear();
	//for (int i = 0; i <= 400; ++i) S.push_back(1000 + i * 10);
	//delta = els.deltaY(S);
	//cout << "deltaY" << endl << "Sy price, delta" << endl;
	//for (auto iter = delta.begin(); iter != delta.end(); ++iter) {
	//	cout << iter->first << ", " << iter->second << endl;
	//}

	// print gammaX
	//S.clear(); delta.clear();
	//for (int i = 0; i <= 400; ++i) S.push_back(1000 + i * 10);
	//delta = els.gammaX(S);
	//cout << "gammaX" << endl << "Sx price, gammaX" << endl;
	//for (auto iter = delta.begin(); iter != delta.end(); ++iter) {
	//	cout << iter->first << ", " << iter->second << endl;
	//}

	// print gammaY
	//S.clear(); delta.clear();
	//for (int i = 0; i <= 400; ++i) S.push_back(1000 + i * 10);
	//delta = els.gammaY(S);
	//cout << "gammaY" << endl << "Sy price, gammaY" << endl;
	//for (auto iter = delta.begin(); iter != delta.end(); ++iter) {
	//	cout << iter->first << ", " << iter->second << endl;
	//}
}

double ELStwo::price(Vector S) { 
	if (!knockInMatrix) solve();
	int x = int(S[0] / dx);
	int y = int(S[1] / dy);
	return (*noKnockInMatrix)[x][y];
}

double ELStwo::delta(Vector S) {
	if (!knockInMatrix) solve();
	int x = int(S[0] / dx);
	int y = int(S[1] / dy);
	double currentPrice = (*noKnockInMatrix)[x][y];
	double upperPriceX = (*noKnockInMatrix)[x+1][y];
	double lowerPriceX = (*noKnockInMatrix)[x - 1][y];
	double upperPriceY = (*noKnockInMatrix)[x][y+1];
	double lowerPriceY = (*noKnockInMatrix)[x][y - 1];
	cout << "deltaX = " << (upperPriceX - currentPrice) / dx << endl;
	cout << "deltaY = " << (upperPriceY - currentPrice) / dy << endl;
	cout << "gammaX = " << (((*noKnockInMatrix)[x + 2][y] - (*noKnockInMatrix)[x + 1][y]) / dx - ((*noKnockInMatrix)[x + 1][y] - (*noKnockInMatrix)[x][y]) / dx) / dx;
	cout << "gammaY = " << (((*noKnockInMatrix)[x][y+2] - (*noKnockInMatrix)[x][y+1]) / dy - ((*noKnockInMatrix)[x][y+1] - (*noKnockInMatrix)[x][y]) / dy) / dy;
	cout << "gammaXY = " << ((*noKnockInMatrix)[x + 1][y + 1] + (*noKnockInMatrix)[x - 1][y - 1] -
		(*noKnockInMatrix)[x + 1][y - 1] - (*noKnockInMatrix)[x - 1][y + 1]) / (4 * dx * dy);
	return 0;
}

std::vector<std::pair<double,double>> ELStwo::deltaX(Vector S) {
	if (!noKnockInMatrix) solve();
	std::vector<std::pair<double, double>> res;
	for (auto iter = S.begin(); iter != S.end(); ++iter) {
		int x = int(*iter / dx);
		int y = int(S0[1] / dy);
		double currentPrice = (*noKnockInMatrix)[x][y];
		double upperPrice = (*noKnockInMatrix)[x + 1][y];
		double delta = (upperPrice - currentPrice) / dx;
		res.push_back(std::pair<double, double>(*iter, delta));
	}
	return res;
}

double ELStwo::deltaX() {
	if (!noKnockInMatrix) solve();
	int x = int(S[0] / dx);
	int y = int(S[1] / dy);
	double upperPrice = (*noKnockInMatrix)[x + 1][y];
	double lowerPrice = (*noKnockInMatrix)[x - 1][y];
	double delta = (upperPrice - lowerPrice) / (2*dx);
	return delta;
}

double ELStwo::deltaY() {
	if (!noKnockInMatrix) solve();
	int x = int(S[0] / dx);
	int y = int(S[1] / dy);
	double upperPrice = (*noKnockInMatrix)[x][y+1];
	double lowerPrice = (*noKnockInMatrix)[x][y-1];
	double delta = (upperPrice - lowerPrice) / (2 * dy);
	return delta;
}

double ELStwo::gammaX() {
	if (!noKnockInMatrix) solve();
	int x = int(S[0] / dx);
	int y = int(S[1] / dy);
	double currentPrice = (*noKnockInMatrix)[x][y];
	double upperPrice = (*noKnockInMatrix)[x + 1][y];
	double lowerPrice = (*noKnockInMatrix)[x - 1][y];
	return (upperPrice + lowerPrice - 2 * currentPrice) / (dx * dx);
}

double ELStwo::gammaY() {
	if (!noKnockInMatrix) solve();
	int x = int(S[0] / dx);
	int y = int(S[1] / dy);
	double currentPrice = (*noKnockInMatrix)[x][y];
	double upperPrice = (*noKnockInMatrix)[x][y+1];
	double lowerPrice = (*noKnockInMatrix)[x][y-1];
	return (upperPrice + lowerPrice - 2 * currentPrice) / (dy * dy);
}

std::vector<std::pair<double, double>> ELStwo::deltaY(Vector S) {
	if (!noKnockInMatrix) solve();
	std::vector<std::pair<double, double>> res;
	for (auto iter = S.begin(); iter != S.end(); ++iter) {
		int x = int(S0[0] / dx);
		int y = int(*iter / dy);
		double currentPrice = (*noKnockInMatrix)[x][y];
		double upperPrice = (*noKnockInMatrix)[x][y+1];
		double delta = (upperPrice - currentPrice) / dy;
		res.push_back(std::pair<double, double>(*iter, delta));
	}
	return res;
}

std::vector<std::pair<double, double>> ELStwo::gammaX(Vector S) {
	if (!noKnockInMatrix) solve();
	std::vector<std::pair<double, double>> res;
	for (auto iter = S.begin(); iter != S.end(); ++iter) {
		int x = int(*iter / dx);
		int y = int(S0[1] / dy);
		double currentPrice = (*noKnockInMatrix)[x][y];
		double upperPrice = (*noKnockInMatrix)[x + 1][y];
		double lowerPrice = (*noKnockInMatrix)[x - 1][y];
		double upperDelta = (upperPrice - currentPrice) / dx;
		double lowerDelta = (currentPrice - lowerPrice) / dx;
		double gamma = (upperDelta - lowerDelta) / (2 * dx);
		res.push_back(std::pair<double, double>(*iter, gamma));
	}
	return res;
}

std::vector<std::pair<double, double>> ELStwo::gammaY(Vector S) {
	if (!noKnockInMatrix) solve();
	std::vector<std::pair<double, double>> res;
	for (auto iter = S.begin(); iter != S.end(); ++iter) {
		int x = int(S0[0] / dx);
		int y = int(*iter / dy);
		double currentPrice = (*noKnockInMatrix)[x][y];
		double upperPrice = (*noKnockInMatrix)[x][y+1];
		double lowerPrice = (*noKnockInMatrix)[x][y-1];
		double upperDelta = (upperPrice - currentPrice) / dy;
		double lowerDelta = (currentPrice - lowerPrice) / dy;
		double gamma = (upperDelta - lowerDelta) / (2 * dy);
		res.push_back(std::pair<double, double>(*iter, gamma));
	}
	return res;
}

std::vector<std::pair<double, double>> ELStwo::vegaX(Vector vol) {
	std::vector<std::pair<double, double>> res;
	double originalSigma = sigma1;
	double dv = sigma1 * 0.001;
	for (auto iter = vol.begin(); iter != vol.end(); ++iter) {
		int x = int(S0[0] / dx);
		int y = int(S0[1] / dy);
		setVol1(*iter);
		updatePrice();
		double currentPrice = (*noKnockInMatrix)[x][y];
		setVol1(*iter + dv);
		updatePrice();
		double upperPrice = (*noKnockInMatrix)[x][y];
		double vega = (upperPrice - currentPrice) / dv;
		res.push_back(std::pair<double, double>(*iter, vega));
	}
	sigma1 = originalSigma;
	reset();
	return res;
}

std::vector<std::pair<double, double>> ELStwo::vegaY(Vector vol) {
	std::vector<std::pair<double, double>> res;
	double originalSigma = sigma2;
	double dv = sigma2 * 0.001;
	for (auto iter = vol.begin(); iter != vol.end(); ++iter) {
		int x = int(S0[0] / dx);
		int y = int(S0[1] / dy);
		setVol2(*iter);
		updatePrice();
		double currentPrice = (*noKnockInMatrix)[x][y];
		setVol2(*iter + dv);
		updatePrice();
		double upperPrice = (*noKnockInMatrix)[x][y];
		double vega = (upperPrice - currentPrice) / (dv);
		res.push_back(std::pair<double, double>(*iter, vega));
	}
	sigma2 = originalSigma;
	reset();
	return res;
}

void ELStwo::updatePrice() {
	reset();
	solve();
}

void ELStwo::solve(){
	if (noKnockInMatrix) return;
	cout << "ELStwo::solve() is called" << endl;
	// calculate boundary(bottom) condition of ELS
	bool knockIn = true;
	std::deque<std::shared_ptr<Matrix>> bcK = initPayoff(knockIn); // boundary condition from maturity to present time
	std::deque<std::shared_ptr<Matrix>> bc = initPayoff(!knockIn); // boundary condition from maturity to present time
	std::shared_ptr<Matrix> vLowK = bcK.front();
	std::shared_ptr<Matrix> vLow = bc.front();
	bcK.pop_front();
	bc.pop_front();
	std::shared_ptr<Matrix> vMidK = createMatrix(Nx + 1, Ny + 1);
	std::shared_ptr<Matrix> vUpK = createMatrix(Nx + 1, Ny + 1);
	std::shared_ptr<Matrix> vMid = createMatrix(Nx + 1, Ny + 1);
	std::shared_ptr<Matrix> vUp = createMatrix(Nx + 1, Ny + 1);
	// make tridiagonal matrix
	// m->m* = A1   // m* -> m+1 = A2
	std::shared_ptr<Matrix> A1 = createCoefMatrix(Nx, sigma1, div1);
	std::shared_ptr<Matrix> A2 = createCoefMatrix(Ny, sigma2, div2);

	Date currDate = redemDate.back();
	// start step
	int kiX = S0[0] * KI / dx;
	int kiY = S0[1] * KI / dy;
	while (currDate >= evalDate) {
		// First step : implicit scheme in x only
		// KNOCK-IN
		for (int y = 1; y <= Ny - 1; ++y) {
			std::shared_ptr<Vector> d = createFirstConstVector(y, vLowK);
			ThomasAlgorithm ta(A1, d);
			std::shared_ptr<Vector> tempVec = ta.solve();
			for (int x = 1; x <= Nx - 1; ++x) (*vMidK)[x][y] = (*tempVec)[x - 1];	// copy to upper matrix
		}
		writeBoundaryCondition(vMidK);

		// NO KNOCK-IN
		for (int y = 1; y <= Ny - 1; ++y) {
			for (int i = 0; i < kiX; ++i) 
				for (int j = 0; j <= Ny; ++j) 
					(*vLow)[i][j] = (*vLowK)[i][j];
			for (int i = kiX; i <= Nx; ++i)
				for (int j = 0; j < kiY; ++j)
					(*vLow)[i][j] = (*vLowK)[i][j];
			std::shared_ptr<Vector> d = createFirstConstVector(y, vLow);
			ThomasAlgorithm ta(A1, d);
			std::shared_ptr<Vector> tempVec = ta.solve();
			for (int x = 1; x <= Nx - 1; ++x) (*vMid)[x][y] = (*tempVec)[x - 1];	// copy to upper matrix
		}
		writeBoundaryCondition(vMid);

		// Second step : implicit scheme in y only
		// KNOCK - IN
		for (int x = 1; x <= Nx - 1; ++x) {
			std::shared_ptr<Vector> d = createSecondConstVector(x, vMidK);
			ThomasAlgorithm ta(A2, d);
			std::shared_ptr<Vector> tempVec = ta.solve();
			std::copy(tempVec->begin(), tempVec->end(), (vUpK->begin() + x)->begin() + 1);	// copy to upper matrix
		}
		writeBoundaryCondition(vUpK);

		// NO KNOCK - IN
		for (int x = 1; x <= Nx - 1; ++x) {
			for (int i = 0; i < kiX; ++i)
				for (int j = 0; j <= Ny; ++j)
					(*vMid)[i][j] = (*vMidK)[i][j];
			for (int i = kiX; i <= Nx; ++i)
				for (int j = 0; j < kiY; ++j)
					(*vMid)[i][j] = (*vMidK)[i][j];
			std::shared_ptr<Vector> d = createSecondConstVector(x, vMid);
			ThomasAlgorithm ta(A2, d);
			std::shared_ptr<Vector> tempVec = ta.solve();
			std::copy(tempVec->begin(), tempVec->end(), (vUp->begin() + x)->begin() + 1);	// copy to upper matrix
		}
		writeBoundaryCondition(vUp);

		// KNOCK - IN
		vLowK = vUpK;
		vMidK = createMatrix(Nx + 1, Ny + 1);
		vUpK = createMatrix(Nx + 1, Ny + 1);
		// NO KNOCK - IN
		vLow = vUp;
		vMid = createMatrix(Nx + 1, Ny + 1);
		vUp = createMatrix(Nx + 1, Ny + 1);

		--currDate;
		if (std::find(redemDate.begin(), redemDate.end(), currDate) != redemDate.end()) {
			printDegreeOfProcess(currDate);
			updatePayoff(vLow, bc); // vLowÀ§¿¡ bc.front()ÀÇ 0¾Æ´Ñ ºÎºÐÀ» µ¤¾î¾º¿ò (¿©±â ±ò²ûÇÏ°ÔÁ» °íÃÄ¾ß ÇÒµí)
			updatePayoff(vLowK, bcK);
		}
	}
	//for (int k = 0; k < K.size(); ++k) {
	//	for (int i = 0; i < pp; ++i) {
	//		// First step : implicit scheme in x only
	//		// KNOCK - IN
	//		for (int y = 1; y <= Ny - 1; ++y) {
	//			std::shared_ptr<Vector> d = createFirstConstVector(y, vLowK);
	//			ThomasAlgorithm ta(A1, d);
	//			std::shared_ptr<Vector> tempVec = ta.solve();
	//			for (int x = 1; x <= Nx - 1; ++x) (*vMidK)[x][y] = (*tempVec)[x - 1];	// copy to upper matrix
	//		}
	//		writeBoundaryCondition(vMidK);

	//		// NO KNOCK - IN
	//		for (int y = 1; y <= Ny - 1; ++y) {
	//			int kiX = S0[0] * KI / dx;
	//			int kiY = S0[1] * KI / dy;
	//			for (int i = 0; i < kiX; ++i) 
	//				for (int j = 0; j <= Ny; ++j) 
	//					(*vLow)[i][j] = (*vLowK)[i][j];
	//			for (int i = kiX; i <= Nx; ++i)
	//				for (int j = 0; j < kiY; ++j)
	//					(*vLow)[i][j] = (*vLowK)[i][j];
	//			std::shared_ptr<Vector> d = createFirstConstVector(y, vLow);
	//			ThomasAlgorithm ta(A1, d);
	//			std::shared_ptr<Vector> tempVec = ta.solve();
	//			for (int x = 1; x <= Nx - 1; ++x) (*vMid)[x][y] = (*tempVec)[x - 1];	// copy to upper matrix
	//		}
	//		writeBoundaryCondition(vMid);

	//		// Second step : implicit scheme in y only
	//		// KNOCK - IN
	//		for (int x = 1; x <= Nx - 1; ++x) {
	//			std::shared_ptr<Vector> d = createSecondConstVector(x, vMidK);
	//			ThomasAlgorithm ta(A2, d);
	//			std::shared_ptr<Vector> tempVec = ta.solve();
	//			std::copy(tempVec->begin(), tempVec->end(), (vUpK->begin() + x)->begin() + 1);	// copy to upper matrix
	//		}
	//		writeBoundaryCondition(vUpK);

	//		// NO KNOCK - IN
	//		for (int x = 1; x <= Nx - 1; ++x) {
	//			int kiX = S0[0] * KI / dx;
	//			int kiY = S0[1] * KI / dy;
	//			for (int i = 0; i < kiX; ++i)
	//				for (int j = 0; j <= Ny; ++j)
	//					(*vMid)[i][j] = (*vMidK)[i][j];
	//			for (int i = kiX; i <= Nx; ++i)
	//				for (int j = 0; j < kiY; ++j)
	//					(*vMid)[i][j] = (*vMidK)[i][j];
	//			std::shared_ptr<Vector> d = createSecondConstVector(x, vMid);
	//			ThomasAlgorithm ta(A2, d);
	//			std::shared_ptr<Vector> tempVec = ta.solve();
	//			std::copy(tempVec->begin(), tempVec->end(), (vUp->begin() + x)->begin() + 1);	// copy to upper matrix
	//		}
	//		writeBoundaryCondition(vUp);

	//		// KNOCK - IN
	//		vLowK = vUpK;
	//		vMidK = createMatrix(Nx + 1, Ny + 1);
	//		vUpK = createMatrix(Nx + 1, Ny + 1);

	//		// NO KNOCK - IN
	//		vLow = vUp;
	//		vMid = createMatrix(Nx + 1, Ny + 1);
	//		vUp = createMatrix(Nx + 1, Ny + 1);
	//	}
	//	printDegreeOfProcess(k, pp);
	//	updatePayoff(vLow, bc); // vLowÀ§¿¡ bc.front()ÀÇ 0¾Æ´Ñ ºÎºÐÀ» µ¤¾î¾º¿ò (¿©±â ±ò²ûÇÏ°ÔÁ» °íÃÄ¾ß ÇÒµí)
	//	updatePayoff(vLowK, bcK);
	//}

	knockInMatrix = vLowK;
	noKnockInMatrix = vLow;
}

void ELStwo::updatePayoff(std::shared_ptr<Matrix>& vLow, std::deque<std::shared_ptr<Matrix>>& bc) {
	int xThreshold = -1, yThreshold = -1;
	bool isThresholdSetted = false;
	if (!bc.empty()) {
		for (int i = 0; i < bc.front()->size(); ++i) {
			if (isThresholdSetted) break;
			for (int j = 0; j < (*bc.front())[i].size(); ++j)
				if ((*bc.front())[i][j] != 0) {
					xThreshold = i;
					yThreshold = j;
					isThresholdSetted = true;
					break;
				}
		}
		for (int i = xThreshold; i < bc.front()->size(); ++i)
			for (int j = yThreshold; j < (*bc.front())[i].size(); ++j)
				(*vLow)[i][j] = (*bc.front())[i][j];
		bc.pop_front();
	}
}

void ELStwo::printDegreeOfProcess(Date currDate) {
	cout << "matrix from " << redemDate.back().to_str() << " to " << currDate.to_str() << " is completed" << endl;
}

// calculate explicit boundary condition (payoff)
std::deque<std::shared_ptr<Matrix>> ELStwo::initPayoff(bool knockIn){
	double F = 100;
	std::deque<std::shared_ptr<Matrix>> bc; // boundary condition from maturity to present time
	auto iter = std::find_if(redemDate.begin(), redemDate.end(), [this](Date d) {return d >= this->evalDate; });
	int bcNum = std::distance(iter, redemDate.end());
	for (int i = 0; i < bcNum; ++i) bc.push_back(createMatrix(Nx + 1, Ny + 1));
	for (int i = 0; i < bc[0]->size(); ++i)
		for (int j = 0; j < (*bc[0])[i].size(); ++j) {
			double X = i * dx, Y = j * dy, min = MIN(X / S0[0], Y / S0[1]);
			if ((min < K.back() && knockIn) || (min < KI && !knockIn)) (*bc[0])[i][j] = min * F;
			else (*bc[0])[i][j] = (1 + cr.back()) * F;
		}
	for (int k = 1; k < bcNum; ++k)
		for (int i = 0; i < bc[k]->size(); ++i)
			for (int j = 0; j < (*bc[k])[i].size(); ++j) {
				double X = i * dx, Y = j * dy, min = MIN(X / S0[0], Y / S0[1]);
				if (min >= K[K.size() - 1 - k]) (*bc[k])[i][j] = (1 + cr[cr.size() - 1 - k]) * F;
			}
	return bc;
}

void ELStwo::writeBoundaryCondition(const std::shared_ptr<Matrix>& v) {
	// ¸ð¼­¸®´Â ÀÏ´Ü ÁÂ¿ì·Î ¶¯±è
	for (int x = 1; x <= Nx-1; ++x) {
		(*v)[x][0] = 2 * (*v)[x][1] - (*v)[x][2];
		(*v)[x][Ny] = 2 * (*v)[x][Ny - 1] - (*v)[x][Ny - 2];
	}
	for (int y = 0; y <= Ny; ++y) {
		(*v)[0][y] = 2 * (*v)[1][y] - (*v)[2][y];
		(*v)[Nx][y] = 2 * (*v)[Nx - 1][y] - (*v)[Nx - 2][y];
	}
}

std::shared_ptr<Vector> ELStwo::createFirstConstVector(int y, const std::shared_ptr<Matrix>& vLow){
	std::shared_ptr<Vector> d = createVector(Nx - 1);
	double commonTerm = rho * sigma1 * sigma2 * 0.125 * dt;
	for (int x = 1; x <= Nx - 1; ++x) {
		(*d)[x - 1] = (*vLow)[x][y] + commonTerm *
			x * y * ((*vLow)[x + 1][y + 1] + (*vLow)[x - 1][y - 1] - (*vLow)[x - 1][y + 1] - (*vLow)[x + 1][y - 1]);
	}

	return d;
}

std::shared_ptr<Vector> ELStwo::createSecondConstVector(int x, const std::shared_ptr<Matrix>& vMid) {
	std::shared_ptr<Vector> d = createVector(Ny - 1);
	double commonTerm = rho * sigma1 * sigma2 * 0.125 * dt;
	for (int y = 1; y <= Ny - 1; ++y)
		(*d)[y - 1] = (*vMid)[x][y] + commonTerm *
		x * y * ((*vMid)[x + 1][y + 1] + (*vMid)[x - 1][y - 1] - (*vMid)[x - 1][y + 1] - (*vMid)[x + 1][y - 1]);
	return d;
}

std::shared_ptr<Matrix> ELStwo::createCoefMatrix(int N, double sigma, double div){
	std::shared_ptr<Vector> a = createVector(N - 1);
	std::shared_ptr<Vector> b = createVector(N - 1);
	std::shared_ptr<Vector> c = createVector(N - 1);
	for (int i = 1; i <= N - 1; ++i) {
		(*a)[i-1] = 1 + dt * (0.5 * r + (r-div) * i + sigma * sigma * i * i);
		(*b)[i-1] = -(sigma * sigma * i * i) * dt * 0.5;
		(*c)[i-1] = (*b)[i-1] - (r-div) * i * dt;
	}
	(*a)[0] += 2 * (*b)[0];
	(*a)[N - 2] += 2 * (*c)[N - 2];
	(*c)[0] -= (*b)[0];
	(*b)[N - 2] -= (*c)[N - 2];
	c->erase(c->end() - 1);
	b->erase(b->begin());
	return createTridiagonalMatrix(a, b, c);
}


ELStwo::ELStwo(Vector S0, Vector S, double sigma1, double sigma2, 
	double rho, double r, double div1, double div2,
	Date evalDate, std::vector<Date> redemDate,
	double TTM, Vector cr, Vector K, double KI, int pp, int Nx, int Ny,
	double S1max, double S2max, double S1min, double S2min)
	: S0(S0), S(S), sigma1(sigma1), sigma2(sigma2), rho(rho), r(r), div1(div1), div2(div2),
	evalDate(evalDate), redemDate(redemDate),
	TTM(TTM), cr(cr), K(K), KI(KI), pp(pp), Nx(Nx), Ny(Ny), M(pp*K.size()),
	S1max(S1max), S2max(S2max), S1min(S1min), S2min(S2min),
	dx((S1max-S1min)/Nx), dy((S2max-S2min)/Ny), 
	dt(TTM / daysBetween(evalDate, redemDate.back())),
	knockInMatrix(nullptr), noKnockInMatrix(nullptr){
}

void ELStwo::reset() {
	knockInMatrix.reset();
	noKnockInMatrix.reset();
}

void ELStwo::setVol1(double vol1) {
	sigma1 = vol1;
	reset();
}

void ELStwo::setVol2(double vol2) {
	sigma2 = vol2;
	reset();
}

void ELStwo::setInterestRate(double r) {
	this->r = r;
	reset();
}

void ELStwo::setCorrelation(double rho) {
	this->rho = rho;
	reset();
}

void ELStwo::setS1(double S1) {
	this->S0[0] = S1;
	reset();
}
void ELStwo::setS2(double S2) {
	this->S0[1] = S2;
	reset();
}
void ELStwo::setDiv1(double div1) {
	this->div1 = div1;
	reset();
}
void ELStwo::setDiv2(double div2) {
	this->div2 = div2;
	reset();
}

Matrix ELStwo::getPriceMatrix() {
	if (!knockInMatrix) solve();
	return *noKnockInMatrix;
}

ELStwo::ELStwo()
{
}

ELStwo::ELStwo(const ELStwo& els)
{
}