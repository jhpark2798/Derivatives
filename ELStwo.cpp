#include "ThomasAlgorithm.h"
#include "FDMutils.h"
#include "ELStwo.h"

#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))

#include <time.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <memory>
#include <deque>
#include <algorithm>

using std::cout;
using std::endl;

void testELStwo()
{
	Vector S0{ 3500, 3000 }; // eurostoxx and sp500
	double sigma1 = 0.25, sigma2 = 0.2;
	double rho = 0.2;
	double r = 0.02;
	double div1 = 0, div2 = 0;
	double TTM = 3;
	Vector cr{ 0.0285, 0.057, 0.0855, 0.114, 0.1425, 0.171 };
	Vector K{ 0.85, 0.85, 0.85, 0.85, 0.85, 0.85 };
	double KI = 0.55;
	int pp = 100;
	int Nx = 300, Ny = 300;
	double S1max = S0[0] *3, S1min = 0;
	double S2max = S0[1] *3, S2min = 0;

	ELStwo els(S0, sigma1, sigma2, rho, r, div1, div2, TTM, cr, K, KI, pp, Nx, Ny,
		S1max, S2max, S1min, S2min);

	clock_t start, end;
	start = clock();
	cout << "price = " << els.price(S0) << endl;
	end = clock();
	cout << "time = " << double(end - start) << "ms" << endl;
}

double ELStwo::price(Vector S) {
	std::unique_ptr<Matrix> knockInMatrix = solve(true);
	std::unique_ptr<Matrix> noKnockInMatrix = solve(false);
	int x = int(S[0] / dx);
	int y = int(S[1] / dy);
	double p = 0.255; // 시뮬레이션 결과 나오는 값으로 변경해야 함
	return p * (*knockInMatrix)[x][y] + (1 - p) * (*noKnockInMatrix)[x][y];
}

std::unique_ptr<Matrix> ELStwo::solve(bool knockIn){
	cout << "ELStwo::solve() is called" << endl;
	// calculate boundary(bottom) condition of ELS
	std::deque<std::unique_ptr<Matrix>> bc = initPayoff(knockIn); // boundary condition from maturity to present time
	std::unique_ptr<Matrix> vLow = std::move(bc.front());
	bc.pop_front();
	std::unique_ptr<Matrix> vMid = createMatrix(Nx + 1, Ny + 1);
	std::unique_ptr<Matrix> vUp = createMatrix(Nx + 1, Ny + 1);
	
	// make tridiagonal matrix
	// m->m* = A1   // m* -> m+1 = A2
	std::unique_ptr<Matrix> A1 = createCoefMatrix(Nx, sigma1);
	std::unique_ptr<Matrix> A2 = createCoefMatrix(Ny, sigma2);

	// start step
	for (int k = 0; k < K.size(); ++k) {
		for (int i = 0; i < pp; ++i) {
			// First step : implicit scheme in x only
			for (int y = 1; y <= Ny - 1; ++y) {
				std::unique_ptr<Vector> d = createFirstConstVector(y, vLow);
				ThomasAlgorithm ta(A1, d);
				std::unique_ptr<Vector> tempVec = ta.solve();
				for (int x = 1; x <= Nx - 1; ++x) (*vMid)[x][y] = (*tempVec)[x - 1];	// copy to upper matrix
			}
			writeBoundaryCondition(vMid);

			// Second step : implicit scheme in y only
			for (int x = 1; x <= Nx - 1; ++x) {
				std::unique_ptr<Vector> d = createSecondConstVector(x, vMid);
				ThomasAlgorithm ta(A2, d);
				std::unique_ptr<Vector> tempVec = ta.solve();
				std::copy(tempVec->begin(), tempVec->end(), (vUp->begin() + x)->begin() + 1);	// copy to upper matrix
			}
			writeBoundaryCondition(vUp);

			vLow = std::move(vUp);
			vUp.reset(); vMid.reset();
			vMid = createMatrix(Nx + 1, Ny + 1);
			vUp = createMatrix(Nx + 1, Ny + 1);
		}
		printDegreeOfProcess(k, pp);
		updatePayoff(vLow, bc); // vLow위에 bc.front()의 0아닌 부분을 덮어씌움 (여기 깔끔하게좀 고쳐야 할듯)
	}
	return vLow;
}

void ELStwo::updatePayoff(std::unique_ptr<Matrix>& vLow, std::deque<std::unique_ptr<Matrix>>& bc) {
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

void ELStwo::printDegreeOfProcess(int k, int pp) {
	cout << (k + 1) * pp << "th step is completed" << endl;
}

std::deque<std::unique_ptr<Matrix>> ELStwo::initPayoff(bool knockIn){
	std::deque<std::unique_ptr<Matrix>> bc; // boundary condition from maturity to present time
	for (int i = 0; i < K.size(); ++i) bc.push_back(createMatrix(Nx + 1, Ny + 1));
	for (int i = 0; i < bc[0]->size(); ++i)
		for (int j = 0; j < (*bc[0])[i].size(); ++j) {
			double X = i * dx, Y = j * dy, min = MIN(X / S0[0], Y / S0[1]);
			if ((min < K.back() && knockIn) || (min < KI && !knockIn)) (*bc[0])[i][j] = min;
			else (*bc[0])[i][j] = 1 + cr.back();
		}
	for (int k = 1; k < bc.size(); ++k)
		for (int i = 0; i < bc[k]->size(); ++i)
			for (int j = 0; j < (*bc[k])[i].size(); ++j) {
				double X = i * dx, Y = j * dy, min = MIN(X / S0[0], Y / S0[1]);
				if (min > K[K.size() - 1 - k]) (*bc[k])[i][j] = 1 + cr[cr.size() - 1 - k];
			}
	return bc;
}

void ELStwo::writeBoundaryCondition(const std::unique_ptr<Matrix>& v) {
	// 모서리는 일단 좌우로 땡김
	for (int x = 1; x <= Nx-1; ++x) {
		(*v)[x][0] = 2 * (*v)[x][1] - (*v)[x][2];
		(*v)[x][Ny] = 2 * (*v)[x][Ny - 1] - (*v)[x][Ny - 2];
	}
	for (int y = 0; y <= Ny; ++y) {
		(*v)[0][y] = 2 * (*v)[1][y] - (*v)[2][y];
		(*v)[Nx][y] = 2 * (*v)[Nx - 1][y] - (*v)[Nx - 2][y];
	}
}

std::unique_ptr<Vector> ELStwo::createFirstConstVector(int y, const std::unique_ptr<Matrix>& vLow){
	std::unique_ptr<Vector> d = createVector(Nx - 1);
	double commonTerm = rho * sigma1 * sigma2 * 0.125 * dt;
	for (int x = 1; x <= Nx - 1; ++x)
		(*d)[x - 1] = (*vLow)[x][y] + 	commonTerm * 
		x * y * ((*vLow)[x + 1][y + 1] + (*vLow)[x - 1][y - 1] - (*vLow)[x - 1][y + 1] - (*vLow)[x + 1][y - 1]);
	return d;
}

std::unique_ptr<Vector> ELStwo::createSecondConstVector(int x, const std::unique_ptr<Matrix>& vMid) {
	std::unique_ptr<Vector> d = createVector(Ny - 1);
	double commonTerm = rho * sigma1 * sigma2 * 0.125 * dt;
	for (int y = 1; y <= Ny - 1; ++y)
		(*d)[y - 1] = (*vMid)[x][y] + commonTerm *
		x * y * ((*vMid)[x + 1][y + 1] + (*vMid)[x - 1][y - 1] - (*vMid)[x - 1][y + 1] - (*vMid)[x + 1][y - 1]);
	return d;
}

std::unique_ptr<Matrix> ELStwo::createCoefMatrix(int N, double sigma){
	std::unique_ptr<Vector> a = createVector(N - 1);
	std::unique_ptr<Vector> b = createVector(N - 1);
	std::unique_ptr<Vector> c = createVector(N - 1);
	for (int i = 1; i <= N - 1; ++i) {
		(*a)[i-1] = 1 + dt * (0.5 * r + r * i + sigma * sigma * i * i);
		(*b)[i-1] = -(sigma * sigma * i * i) * dt * 0.5;
		(*c)[i-1] = (*b)[i-1] - r * i * dt;
	}
	(*a)[0] += 2 * (*b)[0];
	(*a)[N - 2] += 2 * (*c)[N - 2];
	(*c)[0] -= (*b)[0];
	(*b)[N - 2] -= (*c)[N - 2];
	c->erase(c->end() - 1);
	b->erase(b->begin());
	return createTridiagonalMatrix(a, b, c);
}


ELStwo::ELStwo(Vector S0, double sigma1, double sigma2, double rho, double r, double div1, double div2,
	double TTM, Vector cr, Vector K, double KI, int pp, int Nx, int Ny,
	double S1max, double S2max, double S1min, double S2min)
	: S0(S0), sigma1(sigma1), sigma2(sigma2), rho(rho), r(r), div1(div1), div2(div2),
	TTM(TTM), cr(cr), K(K), KI(KI), pp(pp), Nx(Nx), Ny(Ny), M(pp*K.size()),
	dx((S1max-S1min)/Nx), dy((S2max-S2min)/Ny), dt(TTM/(double(pp)*K.size())),
	S1max(S1max), S2max(S2max), S1min(S1min), S2min(S2min)
{
}

ELStwo::ELStwo()
{
}
