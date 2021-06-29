#include "refine.hpp"
#include "MS2D.h"


double Refiner::_h_g1_error = 1e-4;
bool   Refiner::_be_verbose = true;
double Refiner::_h_t0_begin = 0.95;
double Refiner::_h_t1_end   = 0.10;


/*
Def:
	refine model to make it g1
Assume:
	_in is g0 continuous loop (:= _in[0].x1 = _in[1].x0)
	single loop
Desc:
	find touching circle near non-g1 points
		always exists, since voronoi diagram will be there, and VD will be circle center
*/
void Refiner::g1(vector<CircularArc>& _in, vector<CircularArc>& _out)
{
	auto& in  = _in;
	auto& out = _out;

	auto err = 1.0 - g1Error();

	std::vector<int> idx; // idx to take care later
	for (int i = 0; i < in.size(); i++)
	{
		// 1. alias
		int i2 = i + 1;
		if (i2 == in.size())
			i2 = 0;

		auto& arc0 = in[i];
		auto& arc1 = in[i2];

		// 2. check g1
		bool g1;
		auto dot = arc0.n1() * arc1.n0();
		if (abs(dot) < err)
			g1 = false;
		else
			g1 = true;

		// 3. make list
		out.push_back(in[i]);

		// if(not g1) {place temp arc / save idx;}
		if (!g1)
		{
			idx.push_back(out.size());
			out.emplace_back();
		}
	}

	// now take care of those that need change
	vector<int> eraseIdx;
	for (auto& i : idx)
	{
		// 4. alias
		auto iPrev = i - 1; //need no err check for iPrev
		auto iNext = i + 1;
		if (iNext == out.size())
			iNext = 0;

		auto& arc0 = out[iPrev];
		auto& arc1 = out[iNext];
		auto& arcNew = out[i];

		// 5. find theta
		// Note : iff (st > 0) normal changes ccw
		auto st = arc0.tan1() ^ arc1.tan0();
		//auto theta = asin(theta);

		Point 
			p0,		//point on arc0 
			p1,		//point on arc2 
			cen;
		auto
			par0 = arc0.param(),
			par1 = arc1.param();
		double t = _h_t0_begin;
		double t1;
		double rad;
		for (int i = 0; i < 10; i++)
		{
			// 6. find point-normal
			auto theta0 = (1 - t) * par0.first + (t)*par0.second;
			p0 = arc0.xt(theta0);
			Point t0 = arc0.tan(theta0);
			Point n0 = st > 0 ? t0.rotate() : -t0.rotate();

			// 7. find touching point
			rad = arc1.maxTouchRad(p0, n0, t1);

			// 8. check end condition
			if (t1 < _h_t1_end || rad < 0.0)
				break;
			else
				t = 1 - (1 - t) * 0.5;

		}

		// error-check when verbose
		if (_be_verbose)
		{
			if(t1 > _h_t1_end)
				std::cerr<<"!!! WARNING: refiner::g1 : t1 is " << t1 << std::endl;
		}
		if (rad < 0.0)
		{
			if (_be_verbose)
			{
				std::cerr << "!!! WARNING: refiner::g1 : rad is " << rad << std::endl
					<< "		rad0, rad1 : " << arc0.cr() << " , " << arc1.cr() << std::endl;
			}
			eraseIdx.push_back(i);
			continue;
		}
		
		// 9. set stuffs

		// 9-1. first arc
		auto theta0 = (1 - t) * par0.first + (t)*par0.second;
		p0 = arc0.xt(theta0);
		arc0.x1() = p0;
		arc0.n1() = (arc0.x1() - arc0.cc()).normalize();

		// 9-2. second arc
		auto theta1 = (1 - t1) * par1.first + (t1)*par1.second;
		p1 = arc1.xt(theta1);
		arc1.x0() = p1;
		arc1.n0() = (arc1.x0() - arc1.cc()).normalize();

		// 9-3. mid  arc
		Point t0 = arc0.tan(theta0);
		Point n0 = st > 0 ? t0.rotate() : -t0.rotate();
		cen = p0 + rad * n0;
		arcNew.cc() = cen;
		arcNew.x0() = arc0.x1();
		arcNew.x1() = arc1.x0();
		arcNew.cr() = (arcNew.cc() - arcNew.x0()).length1();
		arcNew.n0() = (arcNew.x0() - arcNew.cc()) / arcNew.cr();
		arcNew.n1() = (arcNew.x1() - arcNew.cc()) / arcNew.cr();
		arcNew.ccw  = (st>0);
	}

	// 10. erase stuff
	for (int i = eraseIdx.size() - 1; i >= 0; i--)
	{
		auto ei = eraseIdx[i];
		out.erase(out.begin() + ei);
	}
}