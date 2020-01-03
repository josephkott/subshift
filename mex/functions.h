#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <cmath>

#include "Point.h"
using namespace NR;
using namespace std;

// NLS with piecewise constant nonlinearity with absence of linear potential.
// Equation for nonlinear stationary modes:
// u_{xx} + \mu u + \sigma(x) u^3 = 0
// where \sigma(x) = -1, x \in [0, L_1) & \sigma = +1, x \in [L_1, L_2),
// L = L_1 + L_2 -- period of nonlinear potential.
// Parameters: [\mu L_1 L_2]
Point<2> f(double x, Point<2> u, double* param)
{
	Point<2> du;
	double mu = param[0], L1 = param[1], L2 = param[2];
    double L = L1 + L2;
    
    double n = floor(x / L);
    double x0 = x - n * L;
    double sigma = (x0 < L1 ? -1 : +1);
    
	du[0] = u[1];
	du[1] = -mu * u[0] - sigma * pow(u[0], 3);
	return du;
}

#endif