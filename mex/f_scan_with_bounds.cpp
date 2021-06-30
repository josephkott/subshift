#include "mex.h"
#include "functions.h"
#include "Solver.h"
using namespace NR;

#define PI 3.14159265358979323846

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//Unpacking
	double *i_param = mxGetPr(prhs[0]),
		*i_tspan = mxGetPr(prhs[1]),
		*i_xspan = mxGetPr(prhs[2]),
		*i_yspan = mxGetPr(prhs[3]),
		*i_M = mxGetPr(prhs[4]),
		*i_N = mxGetPr(prhs[5]),
		*i_xBound = mxGetPr(prhs[6]),
		*i_yBound = mxGetPr(prhs[7]);

	const double defIntCount = 32;
	const double* pIntCount;

	if (nrhs == 9)
		pIntCount = mxGetPr(prhs[8]);
	else
		pIntCount = &defIntCount;

	//Preparing
	int intCount = (int)(*pIntCount);
	int M = (int)(*i_M), N = (int)(*i_N);
	double xBound = (double)(*i_xBound);
	double yBound = (double)(*i_yBound);
	double xstep = (i_xspan[1] - i_xspan[0]) / (M-1);
	double ystep = (i_yspan[1] - i_yspan[0]) / (N-1);
	Solver<2> Slv;

	//Scaning
	plhs[0] = mxCreateDoubleMatrix(N, M, mxREAL);
	double *o_Plane = mxGetPr(plhs[0]);

	double x, y;
	double *i_Y0 = new double[2];
	
	for (int i = 0; i < N; i++)
	{
		y = i_yspan[0] + i*ystep;
		for (int j = 0; j < M; j++)
		{
			x = i_xspan[0] + j*xstep;
			
			i_Y0[0] = x;
			i_Y0[1] = y;
			Point<2> Y0(i_Y0);
			
			OdeProblem<2> OP(f, i_param, Y0, i_tspan);
			Solution<2> Sol = Slv.SolveWithBounds(OP, intCount, xBound, yBound);
			o_Plane[i + N*j] = (double)Sol.isInfinit;
		}
	}
	
	delete[] i_Y0;
	return;
}