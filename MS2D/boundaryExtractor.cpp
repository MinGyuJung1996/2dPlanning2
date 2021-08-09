#include "boundaryExtractor.hpp"
#include "MS2D.h"

/*
DeF:
	struct for saving intersection-point info
*/
class ipt
{
public:
	double t; // [0, 1]
	bool   f; // true iff forward is inner region

	inline bool operator<(const ipt& rhs) const
	{
		return this->t < rhs.t;
	}
};



/*
Def:
	do the actual boundary extraction process
Assume:
	all loops are ccw (inner-region on left side of tangent)
	no triple intersection pt of arcs;
Summary:
	1. Find all intersection btw arcs
	2. Subdivide/get rid of those which fall in the inner region
	3. Find neighbor relationship
	4. Find those which forms loop
*/
void bex::calculate()
{
	// 0. alias
	auto& in = *_in;
	auto& out = *_out;

	vector<set<ipt>> interPts(in.size());
	// 1. find intersections & reduce
	//for (all pair)
	for (int i = 0  ; i < in.size(); i++)
	for (int j = i+1; j < in.size(); j++)
	{
		// 1-1. get intersection
		auto& arc0 = in[i];
		auto& arc1 = in[j];
	}

}

/*
Def:
	given two circularArcs which are from a globally-ccw loop, cut off parts that doe not form the boundary
Warning:
	changes input param.

*/
void reduceArcPair(CircularArc& _in_arc0, CircularArc& _in_arc1)
{
	// find intersections
	int intersectionType; 
		// 0: no 
		// 1: one intersection btw arcs / two inter btw circles
		// 2: two intersection
		// 3: one inter btw 
}

/*
Def:
	return truee iff two arcs have common theta interval
*/
bool isArcPairCommonTheta(CircularArc& _in_arc0, CircularArc& _in_arc1)
{
	auto par0 = _in_arc0.param();
	auto par1 = _in_arc1.param();

	auto swap = [](double& a, double& b) -> void
	{
		auto temp = a;
		a = b;
		b = temp;
	};

	// make ccw
	if (par0.first > par0.second)
		swap(par0.first, par0.second);
	if (par1.first > par1.second)
		swap(par1.first, par1.second);

	// confine interval and test
	// 1.
	{
		auto t = par1.first;
		while (t < par0.first)
			t += PI2;
		while (t > par0.first + PI2)
			t -= PI2;

		if (t < par0.second)
			return true;
	}
	{
		auto t = par0.first;
		while (t < par1.first)
			t += PI2;
		while (t > par1.first + PI2)
			t -= PI2;

		if (t < par1.second)
			return true;
	}

	return false;


}

/*
Def: 
	find ipts
Assume:
	r is not 0
Ret:
	_out_special_case:
		0: no problem
		1: the two arcs: have same cc(), cr() & intersect
		2: r = 0
*/
void findIpt(CircularArc& _in_arc0, CircularArc& _in_arc1, vector<ipt>& _out_ipt0, vector<ipt>& _out_ipt1, int& _out_special_case)
{
	constexpr double _h_eq_pt_dist = 1e-8;


	auto& a0 = _in_arc0;
	auto& a1 = _in_arc1;

	_out_special_case = 0;
	
	// 1. check circle intersection
	auto v = a1.cc() - a0.cc();
	double centerDist = (v).length1();
	{
		// i. too far
		if (centerDist >= a0.cr() + a1.cr())
			return;

		// ii. degenerate (same center)
		if (centerDist < _h_eq_pt_dist)
		{
			// if(same rad)
			if (abs(a0.cr() - a1.cr()))
			{
				if (isArcPairCommonTheta(a0, a1))
				{
					if (BOUNDARY_EXTRACTOR_DO_ERROR_CHECK)
					{
						std::cerr << "special case 1 in findIpt" << std::endl;
					}

					_out_special_case = 1;
					return;
				}
			}
			else
				return;

		}

		// iii. r same
		if (BOUNDARY_EXTRACTOR_DO_ERROR_CHECK)
		{
			if (a0.cr() < 1e-100 || a1.cr() < 1e-100)
				std::cerr << "special case 2 in findIpt" << std::endl;
		}

	}

	// 2. find circle intersection pts; (exists)
	Point p0, p1; // Warning, order of intersection is not p0 than p1 for each arc
	double theta0, theta1;
	{
		double theta;

		// 2-1. use cos-law
		auto a = a1.cr();
		auto b = centerDist;
		auto c = a0.cr();

		double cosTheta = (b*b + c*c - a*a) / (2*b*c);
		theta = acos(cosTheta);
		theta0 = theta1 = atan2(v.y(), v.x());
		theta0 += theta;
		theta1 -= theta;

		p0 = a0.xt(theta0);
		p1 = a0.xt(theta1);
	}

	// 3. get arc concavity;
	// as loop(global) is ccw...
	bool a0convex = a0.lccw();
	bool a1convex = a1.lccw();

	// 4. check corresponding param
	double
		arc0_t0, // t-value of p0 of arc0
		arc0_t1, // t-value of p1 of arc0
		arc1_t0,
		arc1_t1;
	auto par0 = a0.param();
	auto par1 = a1.param();
	{
		//auto getT = [](std::pair<double, double>& par, )
		//
		//if (a0.lccw())
		//{
		//	while
		//}
		//else
		//{
		//
		//}
	}

	// 5. check whether p0/p1 is valid for both theta
	bool valid0, valid1;
}


