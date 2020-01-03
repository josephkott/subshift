#ifndef SOLVER_H_
#define SOLVER_H_

#include <cstring>
#include <iostream>
// #include <initializer_list>

#include "Point.h"
#include "Solution.h"
#include "OdeProblem.h"

namespace NR
{

//const double UPPER_BOUND = 10000000000000.0;
const double UPPER_BOUND = 100000.0;
	
template <int dim>
class Solver
{
public:
	Solution<dim> Solve(OdeProblem<dim> OP, const int intCount)
	//Solution<dim> Solve(OdeProblem<dim> OP, const double step)
	{
		Solution<dim> sol;
		sol.nodesCount = intCount + 1;
		/*
		int intCount = floor((OP.tspan[1] - OP.tspan[0]) / step);
		sol.nodesCount = intCount + 1;
		sol.step = (OP.tspan[1] - OP.tspan[0]) / intCount; 
		*/
		sol.step = (OP.tspan[1] - OP.tspan[0]) / intCount;
		sol.points = new Point<dim>[sol.nodesCount];
		sol.isInfinit = false;
		sol.points[0] = OP.Init;
		sol.tspan = new double[2];
		sol.tspan[0] = OP.tspan[0];
		sol.tspan[1] = OP.tspan[1];

		double t = sol.tspan[0];
		double h = sol.step;
		Point<dim> k[4];
		for (int i = 0; i < intCount; i++)
		{
			k[0] = OP.f(t, sol.points[i]) * h;
			k[1] = OP.f(t + 0.5 * h, sol.points[i] + k[0] * 0.5) * h;
			k[2] = OP.f(t + 0.5 * h, sol.points[i] + k[1] * 0.5) * h;
			k[3] = OP.f(t + h, sol.points[i] + k[2]) * h;

			sol.points[i + 1] = sol.points[i] + (k[0] + k[1] * 2.0 + k[2] * 2.0 + k[3]) * (1.0 / 6.0);

			if (sol.points[i + 1].abs() > UPPER_BOUND)
			{
				sol.isInfinit = true;
				sol.tspan[1] = t;
				sol.nodesCount = i;
				return sol;
			}
			t += h;
		}
		return sol;
	}
};

//class Timer
/*{
public:
	Timer()
	{
		t = clock();
		std::cout<<"Start"<<std::endl;
	}
	~Timer()
	{
		std::cout<<"Done!"<<std::endl;
		std::cout<<"Elapsed time:"<< (double)(clock() - t) / CLOCKS_PER_SEC<<std::endl;
	}
private:
	clock_t t;
};*/

}
#endif
