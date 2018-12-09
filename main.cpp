#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <cassert>
#include <numeric>
#include <random>

#include "include/matplotlibcpp.h"

namespace plt = matplotlibcpp;

enum Alpha_method
{
	LIPSCHITZ = 0,
	DIVISION = 1,
	END = 2
};

double getLipschitzConstant(double c0, double c1, double eps1 = 0.1, double eps2 = 0.1)
{
	double lip = std::sqrt(2.0 * c0 * c1);
	std::default_random_engine gen;
	std::uniform_real_distribution<double> dist(eps1, 2.0 / (lip + 2.0 * eps2));
	return dist(gen);
}

std::vector<double> solve(const std::vector<double> &alphas,
						  const std::vector<double> &bettas,
						  const std::vector<double> &gammas,
						  const std::vector<double> &b)
{
	unsigned n = b.size();
	std::vector<double> result(n), tmp1(n), tmp2(n);

	// Go forward
	double y = bettas[0];
	tmp1[0] = -gammas[0] / y;
	tmp2[0] = b[0] / y;

	for (unsigned i = 1; i < n; i++)
	{
		y = bettas[i] + alphas[i] * tmp1[i - 1];
		tmp1[i] = -gammas[i] / y;
		tmp2[i] = (b[i] - alphas[i] * tmp2[i - 1]) / y;
	}

	// Go backward
	result[n - 1] = (b[n - 1] - alphas[n - 1] * tmp2[n - 2]) / (bettas[n - 1] + alphas[n - 1] * tmp1[n - 2]);
	for (int i = n - 2; i >= 0; i--)
	{
		result[i] = tmp1[i] * result[i + 1] + tmp2[i];
	}

	return result;
}

bool test_1()
{
	auto solution = solve({0, 3, 1, 1}, {5, 6, 4, -3}, {3, 1, -2, 0}, {8, 10, 3, -2});
	for (unsigned i = 0; i < solution.size(); i++)
	{
		if (solution[i] - 1.0 > 0.0001)
			return false;
	}
	return true;
}

bool test_2(Alpha_method method)
{
	// Size of grid along time axis
	const unsigned TIME_GRID_SIZE = 1000;

	// Size of grid along rod length axis
	const unsigned POS_GRID_SIZE = 200;

	// Time in seconds
	const double T = 1.0;

	// Length of rod in meters
	const double L = 1.0;

	// Time step
	const double tau = T / (TIME_GRID_SIZE - 1);

	// Length step
	const double h = L / (POS_GRID_SIZE - 1);

	// Coeff
	const double nu = 1.0;

	const double a = 1.0;

	const double p_min = -100;
	const double p_max = 100;

	//const double R = 50.0;

	const double J_EPS = 0.0001;

	const double C0 = std::max((std::pow(a, 4.0) * nu * nu + 2.0 * L) / (a * a * nu), 2.0 * L / (a * a));
	const double C1 = std::max(a * a * nu / 0.1, 1.0 / (a * a * 0.1));

	// Density of heat sources in rod
	auto f = [](double s, double t) -> double {
		return 2.0;
	};

	// Outside temperature
	auto p = [](double t) -> double {
		return 3.0 + t;
	};

	// Descrete form of p
	std::vector<double> p_discrete(TIME_GRID_SIZE);
	for (unsigned i = 0; i < TIME_GRID_SIZE; i++)
	{
		p_discrete[i] = p(i * tau);
	}

	// phi(s) = x(s, 0), t == 0
	auto phi = [](double s) -> double {
		return std::exp( s );
	};

	// Required distribution of heat in rod
	auto y = [](double s) -> double {
		return s * s - std::cos( s );
	};

	std::vector<std::vector<double>> f_discrete(TIME_GRID_SIZE, std::vector<double>(POS_GRID_SIZE));
	for (unsigned i = 0; i < TIME_GRID_SIZE; i++)
	{
		for (unsigned j = 0; j < POS_GRID_SIZE; j++)
		{
			f_discrete[i][j] = f(h * j, tau * i);
		}
	}

	// Error functional
	auto J = [POS_GRID_SIZE, h](const std::vector<double> &x, const std::function<double(double)> &y) -> double {
		double sum = 0;
		for (unsigned i = 0; i < POS_GRID_SIZE; i++)
		{
			sum += (x[i] - y(i * h)) * (x[i] - y(i * h));
		}

		return sum * h;
	};

	std::vector<double> s_values(POS_GRID_SIZE);
	for (unsigned i = 0; i < s_values.size(); i++)
		s_values[i] = i * h;

	std::vector<double> x(POS_GRID_SIZE, 0.0);

	// Here we have SLE with trigiagonal matrix. bettas contains elements from main diagonal.
	// alphas - from bottom diagonal, gammas - from upper diagonal.
	std::vector<double> alphas(POS_GRID_SIZE);
	alphas[0] = 0.0;
	for (unsigned j = 1; j < POS_GRID_SIZE - 1; j++)
	{
		alphas[j] = -a * a * tau / (h * h);
	}
	alphas[POS_GRID_SIZE - 1] = -1.0;

	std::vector<double> bettas(POS_GRID_SIZE);
	bettas[0] = 1.0 + h * h / (2.0 * a * a * tau);
	for (unsigned j = 1; j < POS_GRID_SIZE - 1; j++)
	{
		bettas[j] = 1.0 + 2.0 * a * a * tau / (h * h);
	}
	bettas[POS_GRID_SIZE - 1] = 1.0 + h * h / (2.0 * a * a * tau) + h * nu;

	std::vector<double> gammas(POS_GRID_SIZE);
	gammas[0] = -1.0;
	for (unsigned j = 1; j < POS_GRID_SIZE - 1; j++)
	{
		gammas[j] = -a * a * tau / (h * h);
	}
	gammas[POS_GRID_SIZE - 1] = 0.0;

	std::vector<double> b(POS_GRID_SIZE);

	double alpha = getLipschitzConstant(C0, C1);

	// Iterate while J( x, y ) > J_EPS
	std::vector<double> x_old = x;
	double alpha_div = 1;
	unsigned iters = 0;
	do
	{
		x_old = x;

		alphas[0] = 0.0;
		alphas[POS_GRID_SIZE - 1] = -1.0;

		bettas[0] = 1.0 + h * h / (2.0 * a * a * tau);
		bettas[POS_GRID_SIZE - 1] = 1.0 + h * h / (2.0 * a * a * tau) + h * nu;

		gammas[0] = -1.0;
		gammas[POS_GRID_SIZE - 1] = 0.0;

		// Step #1. Solve diff eq for each time slice:
		// dx/dt = a^2 * d2x/ds2 + f(s, t)
		//  dx/dt(s = 0) = 0
		//  dx/ds(s = l) = nu( p(t) - x(l) )
		//  x(s, 0) = φ(s)

		// Apply last requirement for x(s, 0)
		for (unsigned i = 0; i < POS_GRID_SIZE; i++)
		{
			x[i] = phi(i * h);
		}

		for (unsigned i = 0; i < TIME_GRID_SIZE; i++)
		{
			b[0] = (f_discrete[i][0] + x[0] / tau) * h * h / (2.0 * a * a);
			for (unsigned j = 1; j < POS_GRID_SIZE - 1; j++)
			{
				b[j] = x[j] + f_discrete[i][j] * tau;
			}
			b[POS_GRID_SIZE - 1] = (f_discrete[i][POS_GRID_SIZE - 1] + x[POS_GRID_SIZE - 1] / tau) * h * h / (2.0 * a * a) + nu * p_discrete[i] * h;

			// Use tridiagonal matrix algorithm (Thomas algorithm) to solve SLE
			x = solve(alphas, bettas, gammas, b);
		}

		// Step #2. Solve another diff eq for each time slice:
		// dψ/dt = -a^2 * d2ψ/ds2
		//  dψ/ds(s = 0, t) = 0
		//  dψ/ds(s = l, t) = -nu * ψ(l, t)
		//  ψ(s, t = T) = 2( x(s, T) - y(s) )
		alphas[0] = 0.0;
		alphas[POS_GRID_SIZE - 1] = 1.0;

		bettas[0] = -1.0 + h * h / (2.0 * a * a * tau);
		bettas[POS_GRID_SIZE - 1] = -1.0 - nu * h + h * h / (2.0 * a * a * tau);

		gammas[0] = 1.0;
		gammas[POS_GRID_SIZE - 1] = 0.0;

		std::vector<double> u0(x.size());
		for (unsigned i = 0; i < x.size(); i++)
			u0[i] = 2.0 * (x[i] - y(i * h));

		std::vector<std::vector<double>> KSI(TIME_GRID_SIZE, std::vector<double>(POS_GRID_SIZE));
		for (int i = TIME_GRID_SIZE - 1; i >= 0; i--)
		{
			KSI[i] = u0;

			b = u0;
			b[0] = u0[0] * h * h / (2.0 * a * a * tau);
			b[POS_GRID_SIZE - 1] = u0[POS_GRID_SIZE - 1] * h * h / (2.0 * a * a * tau);

			u0 = solve(alphas, bettas, gammas, b);
		}

		// Step #3. Make a step of the method of conditional gradient
		alpha_div = J(x, y) > J(x_old, y) && iters > 0 ? alpha_div / 1.1 : alpha_div;
		for (unsigned i = 0; i < p_discrete.size(); i++)
		{
			// KSI[i][POS_GRID_SIZE - 1] = ψ(l, t)
			double tmp = KSI[i][POS_GRID_SIZE - 1] > 0 ? p_min : p_max;
			double alpha_lip = std::min(1.0, alpha * a * a * nu * std::abs(a * a * nu * KSI[i][POS_GRID_SIZE - 1] * (tmp - p_discrete[i])) / ((tmp - p_discrete[i]) * (tmp - p_discrete[i])));
			double alpha_k = 0;
			switch (method)
			{
			case LIPSCHITZ:
				alpha_k = alpha_lip;
				break;
			case DIVISION:
				alpha_k = alpha_div;
				break;
			default:
				assert(false);
			}
			p_discrete[i] = p_discrete[i] + alpha_k * (tmp - p_discrete[i]);
		}

		// Step #4. Make a step of gradient projection by f
		// TODO: do later
		std::cout << "J(x, y) = " << J(x, y) << std::endl;
		iters++;
	} while (J(x, y) > J_EPS);

	std::cout << "Iters: " << iters << ", J_EPS: " << J_EPS << std::endl;
	switch (method)
	{
	case LIPSCHITZ:
		std::cout << "Constant of gradient descent was choosed by Lipschitz method" << std::endl;
		break;
	case DIVISION:
		std::cout << "Constant of gradient descent was choosed by division method" << std::endl;
		break;
	default:
		assert(false);
	}

	plt::plot(s_values, x, "r-", s_values, y, "k-");
	plt::save("xy.jpg");
	plt::clf();

	std::vector<double> t_discrete( TIME_GRID_SIZE );
	for( unsigned i = 0; i < TIME_GRID_SIZE; i++ )
	{
		t_discrete[i] = i * tau;
	}
	plt::plot(t_discrete, p_discrete, "r-");
	plt::save("p.jpg");

	return J(x, y) < J_EPS;
}

int main(int argc, char **argv)
{
	Alpha_method method = DIVISION;

	assert(test_1());
	assert(test_2(method));

	return 0;
}