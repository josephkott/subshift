#include "mex.h"
#include "functions.h"
#include "Solver.h"
using namespace NR;

#define PI 3.14159265358979323846

// MATLAB CALL:
// [xgrid, solution] =  f_solve(params, xspan, y0, intCount)

// INPUT:
// :params:   [\omega \alpha]
// :xspan:	  spatial interval of solution
// :y0:		  initial values for Cauchy problem
// :intCount: number of intervals

// OUTPUT:
// :xgrid:	   massive of spatial coordinates
// :solutions: n-dimentioonal massive of solution values on xgrid 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//Unpacking
	double *i_param = mxGetPr(prhs[0]),
		*i_xspan = mxGetPr(prhs[1]),
		*i_U0 = mxGetPr(prhs[2]);

	const double defIntCount = 1024;
	const double* pIntCount;

	if (nrhs == 4)
		pIntCount = mxGetPr(prhs[3]);
	else
		pIntCount = &defIntCount;

	//Preparing
	Point<2> U0(i_U0);
	int intCount = (int)(*pIntCount);

	//Solving
	OdeProblem<2> OP(f, i_param, U0, i_xspan);
	Solver<2> Slv;
	Solution<2> Sol = Slv.Solve(OP, intCount);

	//Packing
	int M = Sol.nodesCount;
	plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(M, 2, mxREAL);
	plhs[2] = mxCreateDoubleScalar((double)Sol.isInfinit);
	double *X = mxGetPr(plhs[0]), *U = mxGetPr(plhs[1]);
	const Point<2>* Data = Sol.GetData();
	for (int i = 0; i < M; i++)
	{
		X[i] = Sol.tspan[0] + i*Sol.step;
		U[i] = Sol[i][0];
		U[i + M] = Sol[i][1];
	}
	return;
}