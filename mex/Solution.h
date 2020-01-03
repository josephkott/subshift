#ifndef SOLUTION_H_
#define SOLUTION_H_

#include <iostream>

namespace NR {

template <int dim>
class Solution
{
public:
	int nodesCount;
	double step;
	Point<dim>* points;
	bool isInfinit;
	double* tspan;
	int* pCopyCounter;

	Solution(const Solution& s)
		: nodesCount(s.nodesCount)
		, step(s.step)
		, points(s.points)
		, isInfinit(s.isInfinit)
		, tspan(s.tspan)
		, pCopyCounter(s.pCopyCounter)
	{
		(*pCopyCounter)++;
	}

	~Solution()
	{
		(*pCopyCounter)--;
		if (*pCopyCounter == 0)
		{
			delete pCopyCounter;
			delete[] points;
		}
	}


	const Point<dim>& operator[](int i)
	{
		return points[i];
	}

	const Point<dim>* GetData() const
	{
		return points;
	}
private:
	Solution()
	{
		pCopyCounter = new int;
		*pCopyCounter = 1;
	}

	Solution& operator=(const Solution& s)
	{
		return *this;
	}

	template <int d1>
	friend class Solver;

	template<int d>
	friend std::ostream& operator<<(std::ostream& stream, Solution<d> sol);
};

template <int d>
std::ostream &operator<<(std::ostream& stream, Solution<d> sol)
{
	for (int i = 0; i < sol.nodesCount; i++)
	{
		stream << sol.points[i] << std::endl;
	}
	return stream;
};

}

#endif