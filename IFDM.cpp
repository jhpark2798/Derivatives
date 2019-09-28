#include "IFDM.h"
#include "ThomasAlgorithm.h"
#include "FDMutils.h"
	
#define MAX(x,y) (((x)>(y)) ? (x) : (y))

#include <time.h>
#include <cmath>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

void testIFDM()
{
	int K = 10;
	double r = 0.05;
	double sigma = 0.2;
	double TTM = 1. / 2.;
	double q = 0;
	int N = 1000;
	int M = 10000;
	double S_max = 20;
	
	clock_t start, end;
	start = clock();

	IFDM ifdm(K, r, sigma, TTM, q, N, M, S_max);
	std::unique_ptr<Vector> price = ifdm.solve();
	
	end = clock();
	cout << "option price is : " << (*price)[N / 2-1] << endl;
	cout << "time = " << double(end-start) << "ms" << endl;
}

IFDM::IFDM()
	: K(0), r(0), sigma(0), TTM(0), q(0), delta_tau(0), delta_x(0), S_max(0), N(0), M(0)
{
}

IFDM::IFDM(int K, double r, double sigma, double TTM, double q,
	int N, int M, double S_max)
	: K(K), r(r), sigma(sigma), TTM(TTM), q(q), S_max(S_max), N(N), M(M), 
	delta_x(S_max/N), delta_tau(TTM/M)
{
}

IFDM::IFDM(const IFDM& ifdm)
	: K(ifdm.K), r(ifdm.r), sigma(ifdm.sigma), TTM(ifdm.TTM), q(ifdm.q), 
	delta_tau(ifdm.delta_tau), delta_x(ifdm.delta_x), S_max(ifdm.S_max), N(ifdm.N), M(ifdm.M)
{
}

IFDM& IFDM::operator=(const IFDM& ifdm)
{
	if (this != &ifdm) {
		K = ifdm.K;
		r = ifdm.r;
		sigma = ifdm.sigma;
		TTM = ifdm.TTM;
		q = ifdm.q;
		delta_tau = ifdm.delta_tau;
		delta_x = ifdm.delta_x;
		S_max = ifdm.S_max;
		N = ifdm.N;
		M = ifdm.M;
	}
	return *this;
}

IFDM::~IFDM()
{
}

std::unique_ptr<Vector> IFDM::solve()
{
	cout << "IFDM::solve() is called" << endl;
	// calculate boundary condition of european put option
	std::unique_ptr<Matrix> v = createMatrix(N + 1, M + 1);
	for (int n = 0; n <= N; ++n) {
		(*v)[n].push_back(MAX(K - (N - n) * delta_x, 0)); // 격자 왼쪽 벽 위부터 채우기
	}
	for (int m = 1; m <= M; ++m) {
		(*v)[0].push_back(0);		// 격자 위쪽 벽 왼쪽부터 채우기
		(*v)[N].push_back(K * std::exp(-r * delta_tau * m) - 0);		// 격자 아래쪽 벽 왼쪽부터 채우기
	}

	// Ax = d를 풀기 위해 A,d를 셋팅. a,b,c,는 A의 대각 벡터들
	std::unique_ptr<Vector> a = createVector(N - 1), b = createVector(N-1), c = createVector(N-1);
	for (int i = 1; i <= N - 1; ++i) {
		a->push_back(1 + delta_tau * (r + sigma * sigma * i * i));
		//b->push_back(delta_tau * 0.5 * ((r - q) * i - sigma * sigma * i * i));
		//c->push_back(delta_tau * 0.5 * (-(r - q) * i - sigma * sigma * i * i));
	}
	for (int i = 1; i <= N - 2; ++i) {
		b->push_back(delta_tau * 0.5 * ((r - q) * (i + 1) - sigma * sigma * (i + 1) * (i + 1)));
		c->push_back(delta_tau * 0.5 * (-(r - q) * i - sigma * sigma * i * i));
	}
	std::unique_ptr<Matrix> A = createTridiagonalMatrix(a, b, c);
	
	c->push_back(delta_tau * 0.5 * (-(r - q) * (N - 1) - sigma * sigma * (N - 1) * (N - 1)));
	b = createVector(N - 1);
	for (int i = 1; i <= N - 1; ++i) {
		b->push_back(delta_tau * 0.5 * ((r - q) * (i + 1) - sigma * sigma * (i + 1) * (i + 1)));
	}

	// d 만들어야 댐. a4에 써논거 필요함.
	std::unique_ptr<Vector> d = createVector(N - 1);
	d->push_back((*v)[N - 1][0] - (*b)[0]*(*v)[N][1]);
	for (int i = 2; i <= N-2; ++i)
		d->push_back((*v)[N - i][0]);
	d->push_back((*v)[1][0] - (*c)[N - 2] * (*v)[0][1]);
	
	for (int m = 1; m <= M; ++m) {
		if (m % 500 == 0) cout << m << "th step is completed" << endl;
		ThomasAlgorithm ta(A, d);
		d = ta.solve();

		for (int i = 1; i <= N - 1; ++i) {
			(*v)[N - i].push_back((*d)[i - 1]);
		}
		if (m != M) {
			(*d)[0] = (*v)[N - 1][m] - (*b)[0] * (*v)[N][m + 1];
			(*d)[N - 2] = (*v)[1][m] - (*c)[N - 2] * (*v)[0][m + 1];
		}
		// printMatrix(v);
	}

	return d;
}

