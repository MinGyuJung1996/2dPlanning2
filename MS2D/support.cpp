#include "support.hpp"
#include "MS2D.h"
#include "omp.h"
#include <iomanip>
#include "collision detection.hpp"
#include "refine.hpp"

#define derr cerr
extern vector<double> debugBlock;

/*
Def: 
	Simple constructor
*/
Support::Support(CircularArc& arc)
{
	build(arc);
}

/*
Def:
	Build internal data from arc.
*/
void
Support::build(CircularArc& arc)
{
	// 1. orig
	_ori = &arc;

	// 2. Assume ccw of arc : needs to change cw arcs
	Point p0, p1;
	if (arc.lccw())
	{
		auto param = arc.param();
		_t0 = param.first;
		_t1 = param.second;
		p0 = arc.x0();
		p1 = arc.x1();
		_flipped = false;
	}
	else
	{
		auto param = arc.param();
		_t0 = param.second;
		_t1 = param.first;
		p0 = x1();
		p1 = x0();
		_flipped = true;
	}

	// set tX
	{
		_t3 = _t0 + PI2;
		_t2 = (_t3 + _t1) / 2;
	}

	// set x, y, r
	{
		_x0 = arc.cx();
		_y0 = arc.cy();
		_r0 = arc.cr();

		_x1 = _x0 + _r0 * cos(_t1);
		_y1 = _y0 + _r0 * sin(_t1);

		_x2 = _x0 + _r0 * cos(_t0);
		_y2 = _y0 + _r0 * sin(_t0);
	}
}


/*
Def: get support value at theta
Param:
	theta : theta
	c : (opt) cos(theta)
	s : (opt) sin(theta)
	type : (opt) type of best support point
		0 : point on arc
		1 : endpoint1
		2 : endpoint0
Ret:
	Support val
*/
double
Support::eval(double theta, double c, double s)
{
	double ret;

	// 1. find type
	auto t = type(theta);

	// 2. get cos/sin, if not given;
	if (c > 1.1)
	{
		c = cos(theta);
		s = sin(theta);
	}

	// 3. switch
	switch (t)
	{
	case 0:
		return f0(c, s);
	case 1:
		return f1(c, s);
	case 2:
		return f2(c, s);
	}
}

/*
Def: returns the type of the best supporting point
Ret:
	0: point on arc
	1: ep1;
	2: ep0;
*/
int
Support::type(double theta)
{
	// 1. confine theta between t0 and (t0+2pi)
	while (theta < _t0)
		theta += PI2;
	while (theta > _t3)
		theta -= PI2;

	// 2. compare interval
	if (theta < _t1)
		return 0;
	if (theta < _t2)
		return 1;
	else
		return 2;
}

/*
Def: evaluate func val for theta : t0 < theta < t1
Param: cos & sin value at theta
*/
double 
Support::f0(double c, double s)
{
	return _x0 * c + _y0 + s + _r0;
}

/*
Def: evaluate func val for theta : t1 < theta < t2
Param: cos & sin value at theta
*/
double 
Support::f1(double c, double s)
{
	return _x1 * c + _y1 + s;
}

/*
Def: evaluate func val for theta : t2 < theta < t3
Param: cos & sin value at theta
*/
double Support::f2(double c, double s)
{
	return _x2 * c + _y2 + s;
}

void Support::print()
{
	using namespace std;
	const int s = 10;

	cout << "*** Support info ***" << endl;
	cout << "a : "
		<< setw(s) << x0() << " "
		<< setw(s) << x1() << " "
		<< setw(s) << x2() << " " << endl;
	cout << "b : "
		<< setw(s) << y0() << " "
		<< setw(s) << y1() << " "
		<< setw(s) << y2() << " " << endl;
	cout << "r : "
		<< setw(s) << r0() << endl;
	cout << "t : "
		<< setw(s) << t0() << " "
		<< setw(s) << t1() << " "
		<< setw(s) << t2() << " "
		<< setw(s) << t3() << " "
		<< endl;

}

/*****************************************************************************************************************************************************
**																																					**
**															Common Tangent Calculator																**
**																																					**
*****************************************************************************************************************************************************/

/*
Def : find common tangents between two circular arc, given its Sup
Desc:
	Consider endpoints as circles/arcs with radius 0, then 4 case
		ep  x ep
		arc x ep
		ep  x arc
		arc x arc
Ret :
	ret[2n], ret[2n+1] forms the common tangent line.
*/
vector<Point> 
CommonTangentCalculator::eval(Sup& a, Sup& b)
{

}

vector<Point> 
CommonTangentCalculator::eval(CircularArc& a, CircularArc& b)
{
	Sup sa(a);
	Sup sb(b);
	return eval(sa, sb);
}

/*
Def: Set input0
*/
void 
CommonTangentCalculator::setInput0(vector<CircularArc>& loop0)
{
	_loop0 = &loop0;
}

/*
Def: Set input1
*/
void 
CommonTangentCalculator::setInput1(vector<CircularArc>& loop1)
{
	_loop1 = &loop1;
}

/*
Def: Set output
*/
void 
CommonTangentCalculator::setOutput(vector<Point>& out)
{
	_out = &out;
}

/*
Def: calculate lineSegments that are common tangents of loop0 and loop1
*/
void 
CommonTangentCalculator::calc()
{
	// 1.
	vector<Sup>
		sup0,
		sup1;
	{
		for (auto& a : *_loop0)
			sup0.emplace_back(a);
		for (auto& a : *_loop1)
			sup1.emplace_back(a);
	}

	// 2.
	vector<SupCom>
		sc0,
		sc1;
	{
		for (auto& s : sup0)
			sc0.emplace_back(s);
		for (auto& s : sup1)
			sc1.emplace_back(s);
	}

	// 3. env
	SupCom
		env0,
		env1;
	{
		env0 = sc0[0];
		env1 = sc1[0];

		for (int i = 1; i < sc0.size(); i++)
			env0 = env0.upperEnvelope(sc0[i]);
		for (int i = 1; i < sc1.size(); i++)
			env1 = env1.upperEnvelope(sc1[i]);
	}

	// 4. case 1;
	auto inter = env0.intersection(env1);
	{
		for (auto t : inter)
		{
			_out->push_back(env0.point(t));
			_out->push_back(env1.point(t));
		}
	}

	// 5. case2
	auto env2 = env1;
	env1.flip();
	inter = env0.intersection(env1);
	{
		for (auto t : inter)
		{
			_out->push_back(env0.point(t));
			_out->push_back(env2.point(t+PI));
		}
	}
}


/*****************************************************************************************************************************************************
**																																					**
**															  Convex Hull Calculator																**
**																																					**
*****************************************************************************************************************************************************/

double CvxHullCalc::_lineRad = 1000.0;
int CvxHullCalc::_nSam = 360;

/*
Def: Register vector<CA> to read
*/
void CvxHullCalc::setInput(vector<CircularArc>& inputArcs)
{
	_iArcs = &inputArcs;
}

/*
Def: Register vector<CA> to write
*/
void CvxHullCalc::setOutput(vector<CircularArc>& outputArcs)
{
	_oArcs = &outputArcs;
}

/*
Def: find cvxHull using sampling
*/
void CvxHullCalc::calcSam()
{
	// alias;
	auto& in  = *_iArcs;
	auto& out = *_oArcs;

	// 1. find cos/sin first
	vector<double> cosArr(_nSam);
	vector<double> sinArr(_nSam);

	for (int i = 0; i < _nSam; i++)
	{
		double theta = double(i) / _nSam * PI2;
		cosArr[i] = cos(theta);
		sinArr[i] = sin(theta);
	}

	// 2. find support for each arc
	vector<Sup> sups(in.size());

	for (int i = 0; i < in.size(); i++)
	{
		sups[i].build(in[i]);
	}

	// 3. find Max support;
	vector<double> maxSup(_nSam, -1e+100);
	vector<int> maxSupIdx(_nSam);

	{
		int len = in.size();
		#pragma omp parallel for
		for (int i = 0; i < _nSam; i++)
		{
			auto theta = double(i) / _nSam * PI2;
			auto c = cos(theta);
			auto s = sin(theta);
			for (int j = 0; j < len; j++)
			{
				auto& arc = in[j];
				auto sup = sups[j].eval(theta, c, s);
				if (sup > maxSup[i])
				{
					maxSup[i] = sup;
					maxSupIdx[i] = j;
				}
			}
		}

	}

	// 4. find type of maxSup
	vector<int> maxSupTyp(_nSam);
	{
		for (int i = 0; i < _nSam; i++)
		{
			maxSupTyp[i] = sups[maxSupIdx[i]].type(double(i) / _nSam * PI2);
		}
	}
	// 4. from max supports, build convex hull
	/* 
	If maxSupIdx[i] != maxSupIdx[i+1];
		simply connect line between it...
		or there maybe no primitive, when g1 cont arcs
	If maxSupIdx[i] == maxSupIdx[i+1];
		it could be
			a point
			a line
			an arc
	Need to merge last&first arc if 
	*/
	for (int i = 0; i < _nSam; i++)
	{
		// 4-1. find arc's idx and support type;
		int idx0, idx1;
		int typ0, typ1;
		{
			idx0 = maxSupIdx[i];
			typ0 = maxSupTyp[i];
			if (i != _nSam - 1)
			{
				idx1 = maxSupIdx[i + 1];
				typ1 = maxSupTyp[i + 1];
			}
			else
			{
				idx1 = maxSupIdx[0];
				typ1 = maxSupTyp[0];
			}
		}

		// 4-2. divide cases
		//if (same arc)
		if (idx0 == idx1)
		{

			// if (same type)
			if (typ0 == typ1)
			{
			//what if (_t3 - _t1) <<<< 360.0 / _nSample ???

			}
			// if (diff type)
			else
			{

			}
		} // ~if (idx0 == idx1)
		else
		{
			if (1)
				1;
		} // ~if (idx0 == idx1)
	}


}

void ConvexHullCalculator::calcAlg()
{
	// 1. find envelope
	vector<Support> sup;
	for (auto& a : *_iArcs)
		sup.emplace_back(a);

	vector<SupportComplex> supc;
	for (auto& s : sup)
		supc.emplace_back(s);

	auto env = supc.front();
	for (auto& sc : supc)
		env = env.upperEnvelope(sc);

	// 2. gen primitive (ccw) circular arcs
	vector<CircularArc> arcsOnly;

	auto& t = env.t();
	auto& iArcs = env.arc();
	for (int i = 0; i < t.size() - 1; i++)
	{
		// 2-1. check param
		auto t0 = t[i + 0];
		auto t1 = t[i + 1];
		//{
		//	// if t0, t1 is both inside original arc's param => just make arc
		//	double theta0, theta1;
		//	auto param = iArcs[i][0].param();
		//	if (iArcs[i][0].lccw())
		//	{
		//		theta0 = param.first;
		//		theta1 = param.second;
		//	}
		//	else
		//	{
		//		theta0 = param.second;
		//		theta1 = param.first;
		//	}
		//
		//}

		// 2-2-1.
		// if (this is arc)
		if (env.r()[i] > 1e-300)
		{
			CircularArc temp = iArcs[i][0];
			temp.ccw = true;
			temp.n0() = Point(cos(t0), sin(t0));
			temp.n1() = Point(cos(t1), sin(t1));
			temp.x0() = temp.cc() + temp.cr() * temp.n0();
			temp.x1() = temp.cc() + temp.cr() * temp.n1();
			arcsOnly.push_back(temp);
		}
		// 2-2-2.
		//else (this is point)
		else
		{
			// !!!! Note : this is not actually an valid arc => need removal at end;
			CircularArc temp;

			temp.cr() = -1.0; // just to note that this is not real arc

			//find which point it is
			Point p;
			{
				auto& arc = iArcs[i][0];
				auto t = (t0+t1)/2.0;
				Point temp = arc.cr() * Point(cos(t), sin(t)) + arc.cc();
				if ((temp - arc.x0()).length2() < (temp - arc.x1()).length2())
					p = arc.x0();
				else
					p = arc.x1();
			}

			temp.x0() = p;
			temp.x1() = p;
			arcsOnly.push_back(temp);
		}

	}

	// 3. check last/first arc is same
	bool lfSame = false;
	{
		auto& a0 = arcsOnly.front();
		auto& a1 = arcsOnly.back();

		bool cond0;
		cond0 = (a1.x1() == a0.x0());

		bool cond1;
		cond1 =	(a0.cc() == a1.cc());

		bool cond2;
		cond2 = fabs(a0.cr() - a1.cr()) < 1e-8;

		lfSame == (cond0 && cond1 && cond2);
	}
	if (lfSame)
	{
		arcsOnly.front().x0() = arcsOnly.back().x0();
		arcsOnly.front().n0() = arcsOnly.back().n0();
		arcsOnly.pop_back();
	}

	// 4. add lines
	auto& out = *_oArcs;
	for (int i = 0; i < arcsOnly.size(); i++)
	{
		// 4-1.
		if(arcsOnly[i].cr() > 1e-10) // if (not eps or error-arc from 2-2-2)
			out.push_back(arcsOnly[i]);

		// 4-2.
		int iNext = i + 1;
		if (iNext == arcsOnly.size())
			iNext = 0;

		// 4-3.
		auto& a0 = arcsOnly[i];
		auto& a1 = arcsOnly[iNext];

		auto& p0 = a0.x1();
		auto& p1 = a1.x0();

		if ((p0 - p1).length2() < 1e-8) // if p0 and p1 are too close, construction of line-circularArc becomes erroneous
			continue;

		// 4-4.
		CircularArc arc = cd::constructArc(p0, p1, lineRad(), true, true);
		out.push_back(arc);
	}
}

/*
Def: offset convex hull
Assume:
	"in" is properly-ordered convex hull
	all arcs are convex arcs
Desc:
	offset radius for original arcs,

*/
void ConvexHullCalculator::offset(vector<CircularArc>& in, vector<CircularArc>& out, double offRad)
{
	for (int i = 0; i < in.size(); i++)
	{
		// 1. find next arc
		int iNext = i + 1;
		if (iNext == in.size())
			//iNext == 0; // bug-fixed : bi2.png
			iNext = 0;

		// 2. offset original Arc
		auto arc = in[i];

		arc.cr() = arc.cr() + offRad;
		arc.x0() = arc.n0() * arc.cr() + arc.cc();
		arc.x1() = arc.n1() * arc.cr() + arc.cc();

		//if (arc.cr() < lineRad())
		//{
		//	// when this is normal arc => change r
		//	arc.cr() = arc.cr() + offRad;
		//	arc.x0() = arc.n0() * arc.cr() + arc.cc();
		//	arc.x1() = arc.n1() * arc.cr() + arc.cc();
		//}
		//else
		//{
		//	//when the arc represents a straight line => translate
		//	Point trans = arc.x1() - arc.x0();
		//	trans.normalize();
		//	trans = -trans.rotate(); // as ccw;
		//
		//	arc.cc() = arc.cc() + trans * offRad;
		//	arc.x0() = arc.x0() + trans * offRad;
		//	arc.x1() = arc.x1() + trans * offRad;
		//}

		out.push_back(arc);

		// 3. check g1
		const double errG1 = 1e-8;
		bool g1;
		{
			auto dot = arc.n1() * in[iNext].n0();
			if (fabs(dot) > 1 - errG1) // note : theoeretically, fabs is not needed.
				g1 = true;
			else
				g1 = false;
		}

		// 4. if(not g1, add intermediate arc)
		if (!g1)
		{
			CircularArc arc;
			arc.cc() = in[iNext].x0();
			arc.cr() = offRad;
			arc.n0() = in[i].n1();
			arc.n1() = in[iNext].n0();
			arc.x0() = arc.n0() * arc.cr() + arc.cc();	
			arc.x1() = arc.n1() * arc.cr() + arc.cc();

			out.push_back(arc);
		}

	}

	// do refinement
	vector<CircularArc> temp;
	Refiner::g1(out, temp);
	out.swap(temp);
}


/*****************************************************************************************************************************************************
**																																					**
**																Support(Complex)																	**
**																																					**
*****************************************************************************************************************************************************/
int SupportComplex::_drawSampleNo = 400;

SupportComplex::SupportComplex(Support& sup)
{
	build(sup);
}

/*
Def: From Support class, build internal data.
*/
void SupportComplex::build(Support& sup)
{
	// alias
	auto& t0 = sup.t0();
	auto& t1 = sup.t1();
	auto& t2 = sup.t2();
	auto& t3 = sup.t3();

	// 1. find where 2n*pi is
	double zero = -2 * PI2;
	while (zero < sup.t0())
		zero += PI2;

	// 2. 3cases
	//	0 := t0 < 2npi < t1 
	//	1 := t1 < 2npi < t2
	//	2 := t2 < 2npi < t3
	int c = 2;
	if (zero < t1)
		c = 0;
	else if (zero < t2)
		c = 1;

	// 3. switch
	switch (c)
	{
	case 0:
		{
		_interval.push_back(0.0);
		_interval.push_back(t1 - zero);
		_interval.push_back(t2 - zero);
		_interval.push_back(t3 - zero);
		_interval.push_back(PI2);

		_a.push_back(sup.x0());
		_a.push_back(sup.x1());
		_a.push_back(sup.x2());
		_a.push_back(sup.x0());

		_b.push_back(sup.y0());
		_b.push_back(sup.y1());
		_b.push_back(sup.y2());
		_b.push_back(sup.y0());

		_r.push_back(sup.r0());
		_r.push_back(0.0);
		_r.push_back(0.0);
		_r.push_back(sup.r0());

		_arcIdx.push_back(&sup.ori());
		_arcIdx.push_back(&sup.ori());
		_arcIdx.push_back(&sup.ori());
		_arcIdx.push_back(&sup.ori());
		}
		break;
	case 1:
		{
		_interval.push_back(0.0);
		_interval.push_back(t2 - zero);
		_interval.push_back(t3 - zero);
		_interval.push_back(t1 - zero + PI2);
		_interval.push_back(PI2);

		_a.push_back(sup.x1());
		_a.push_back(sup.x2());
		_a.push_back(sup.x0());
		_a.push_back(sup.x1());

		_b.push_back(sup.y1());
		_b.push_back(sup.y2());
		_b.push_back(sup.y0());
		_b.push_back(sup.y1());

		_r.push_back(0.0);
		_r.push_back(0.0);
		_r.push_back(sup.r0());
		_r.push_back(0.0);

		_arcIdx.push_back(&sup.ori());
		_arcIdx.push_back(&sup.ori());
		_arcIdx.push_back(&sup.ori());
		_arcIdx.push_back(&sup.ori());
		}
		break;
	case 2:
	default:
	{
		_interval.push_back(0.0);
		_interval.push_back(t3 - zero);
		_interval.push_back(t1 - zero + PI2);
		_interval.push_back(t2 - zero + PI2);
		_interval.push_back(PI2);

		_a.push_back(sup.x2());
		_a.push_back(sup.x0());
		_a.push_back(sup.x1());
		_a.push_back(sup.x2());

		_b.push_back(sup.y2());
		_b.push_back(sup.y0());
		_b.push_back(sup.y1());
		_b.push_back(sup.y2());

		_r.push_back(0.0);
		_r.push_back(sup.r0());
		_r.push_back(0.0);
		_r.push_back(0.0);

		_arcIdx.push_back(&sup.ori());
		_arcIdx.push_back(&sup.ori());
		_arcIdx.push_back(&sup.ori());
		_arcIdx.push_back(&sup.ori());
	}
		break;
	}

	// check redundancy
	{
		// (fabs < 1e-8) ??
		if (_interval.back() == _interval[_interval.size() - 2])
		{
			_interval.pop_back();
			_a.pop_back();
			_b.pop_back();
			_r.pop_back();
			_arcIdx.pop_back();
		}
	}

}

/*
Def: draw graph
*/
void SupportComplex::draw()
{
	// 1. build p
	vector<Point> p(_drawSampleNo+ 1);

	// interval idx, where t0 < t < t1 holds
	int intIdx = 0; 
	double t0 = _interval[intIdx + 0];
	double t1 = _interval[intIdx + 1];

	// for i = [0, dsn)
	for (int i = 0; i < _drawSampleNo; i++)
	{
		// t as in C(t)
		double t = double(i) / _drawSampleNo * PI2;

		// check whether t0 < t < t1 holds
		while (t > t1)
		{
			intIdx++;
			t0 = _interval[intIdx + 0];
			t1 = _interval[intIdx + 1];
		}

		// evaluate point
		p[i].x() = t - PI;
		p[i].y() = _a[intIdx] * cos(t) + _b[intIdx] * sin(t) + _r[intIdx];

	}

	// for i = [dsn, dsn]
	{
		p[_drawSampleNo].x() = PI2 - PI;
		p[_drawSampleNo].y() = _a.back() + _r.back(); //as c = 1, s = 0
	}


	// 2. rendering
	{
		glBegin(GL_LINE_STRIP);
		for (auto& xy : p)
			glVertex2dv(xy.P);
		glEnd();
	}
}

/*
Def: 
	given two SupCom, compute a new SupCom, which has its upper envelope
Hyper:
	eraseEps 
Date:
	210618
*/
SupportComplex SupportComplex::upperEnvelope(SupportComplex& rhs)
{
	// 0. Alias
	auto& lhs = *this;
	SupCom ret, temp;

	// 0.5 quick test
	{
		// min max comparison
	}

	// 1. build interval
	set<double> interv;
	for (auto i : lhs._interval)
		interv.insert(i);
	for (auto i : rhs._interval)
		interv.insert(i);

	// 2. loop variables
	if (interv.size() == 0)
	{
		std::cerr << "!!! ERR in upperEnvelope() : interv.size() == 0" << endl;
	}
	auto it0 = interv.begin();
	auto it1 = it0;
	it1++;

	int idx0 = -1;
	int idx1 = -1;

	// 3. compare curves in interval
	// find intersections and btw intersections, compare size
	while (it1 != interv.end())
	{
		// i. interval is [t0, t1)
		auto t0 = *it0;
		auto t1 = *it1;

		// ii. set idx_
		// if(t0 is inside lhs.interval) idx0++;
		//if (find(lhs._interval.begin(), lhs._interval.end(), t0) != lhs._interval.end())
		//	idx0++;
		//if (find(rhs._interval.begin(), rhs._interval.end(), t0) != rhs._interval.end())
		//	idx1++;
		if (lhs._interval[idx0 + 1] == t0)
			idx0++;
		if (rhs._interval[idx1 + 1] == t0)
			idx1++;
		
		// flag for no intersection : if (true) simply pick a larger one
		bool noInter = false;

		// iii. find coeff
		double 
			&a0 = lhs._a[idx0], 
			&b0 = lhs._b[idx0], 
			&r0 = lhs._r[idx0],
			&a1 = rhs._a[idx1],
			&b1 = rhs._b[idx1],
			&r1 = rhs._r[idx1];

		double
			da = a0 - a1,
			db = b0 - b1,
			dr = r0 - r1;

		

		// iv. solve da(cost) + db(sint) + dr = 0;
		// change it to cos(t - phi) = -dr/sqrt(...)
		double len = sqrt(da * da + db * db);
		double con = 2.0, conTheta;
		{
			if (len == 0.0)
				//only r is different => infinite intersection is though of as no inter
				noInter = true;
			else
				con = -dr / len;

			if (fabs(con) > 1.0)
				noInter = true;
			else
				conTheta = acos(con); // [0, pi]
		}

		// v. add intervals to temp
		auto s0 = temp._interval.size();
		temp._interval.push_back(t0);

		if (!noInter)
		{
			double phi = atan2(db, da); // [-pi,pi]
			double
				theta0 = phi - conTheta, // [-2pi, pi]
				theta1 = phi + conTheta; // [-pi, 2pi]

			//dbg_out
			{
				
			}

			// turn them to [0, 2pi]
			if (theta0 < 0.0)
				theta0 += PI2;
			if (theta1 < 0.0)
				theta1 += PI2;

			// check whether they are valid and insert
			// bug-fixed : no guarantee that order is : 1.theta0 2. theta1
			if (theta0 < theta1)
			{
				if (t0 < theta0 && theta0 < t1)
					temp._interval.push_back(theta0);
				if (t0 < theta1 && theta1 < t1)
					temp._interval.push_back(theta1);
			}
			else
			{
				if (t0 < theta1 && theta1 < t1)
					temp._interval.push_back(theta1);
				if (t0 < theta0 && theta0 < t1)
					temp._interval.push_back(theta0);
			}
		}
		s0 = temp._interval.size() - s0; //number of newly added elements

		// vi. find a b r for intrval
		temp._interval.push_back(t1); // just to make the code better
		for (int i = s0; i > 0; i--)
		{
			// vi-1. get t
			auto t0 = temp._interval[temp._interval.size() - i - 1];
			auto t1 = temp._interval[temp._interval.size() - i];
			auto t = (t0 + t1) / 2;

			// vi-2. compar who is bigger at t
			auto f0 = a0 * cos(t) + b0 * sin(t) + r0;
			auto f1 = a1 * cos(t) + b1 * sin(t) + r1;

			if (f0 > f1)
			{
				temp._a.push_back(a0);
				temp._b.push_back(b0);
				temp._r.push_back(r0);
				temp._arcIdx.push_back(lhs._arcIdx[idx0]);
			}
			else
			{
				temp._a.push_back(a1);
				temp._b.push_back(b1);
				temp._r.push_back(r1);
				temp._arcIdx.push_back(rhs._arcIdx[idx1]);
			}

		}
		temp._interval.pop_back(); // as this element will be inserted twice, otherwise



		it0++;
		it1++;
	}
	temp._interval.push_back(PI2);

	//if (debugBlock.size() != 0)
	//{
	//	cout << "idx0 idx1" << endl;
	//	cout << idx0 << "  " << lhs._interval.size() << endl;
	//	cout << idx1 << "  " << rhs._interval.size() << endl;
	//}

	ret = _simplify(temp);

	// 6. set max min.
	{

	}

	return ret;
}

void SupportComplex::print()
{
	using namespace std;
	const int s = 10;

	cout << "*** SupportComplex info ***" << endl;
	cout << "a : ";
	for (auto a : _a)
		cout << setw(s) << a << " ";
	cout << endl;
	
	cout << "b : ";
	for (auto a : _b)
		cout << setw(s) << a << " ";
	cout << endl;

	cout << "a : ";
	for (auto a : _r)
		cout << setw(s) << a << " ";
	cout << endl;

	cout << "t : ";
	for (auto a : _interval)
		cout << setw(s) << a/PI << "p ";
	cout << endl;

}

/*
Def:
	erase redudant intervals / too small intervals
	and save it to *this
Assume:
	
*/
SupportComplex SupportComplex::_simplify(SupportComplex& temp)
{
	SupportComplex ret;

	// 4. merge intervals that are the same
	//	Same := same a, b, r
	// simply insert the first a,b,r,t 
	{
		ret._interval.push_back(temp._interval.front());
		ret._a.push_back(temp._a.front());
		ret._b.push_back(temp._b.front());
		ret._r.push_back(temp._r.front());
		ret._arcIdx.push_back(temp._arcIdx.front());

		for (int i = 1; i < temp._a.size(); i++)
		{

			// if(same interval content) skip
			if (ret._a.back() == temp._a[i] &&
				ret._b.back() == temp._b[i] &&
				ret._r.back() == temp._r[i] &&
				ret._arcIdx.back() == temp._arcIdx[i])
				continue;

			ret._a.push_back(temp._a[i]);
			ret._b.push_back(temp._b[i]);
			ret._r.push_back(temp._r[i]);
			ret._arcIdx.push_back(temp._arcIdx[i]);
			ret._interval.push_back(temp._interval[i]);
		}

		ret._interval.push_back(PI2);
	}

	// 5. erase intervals with 0 length
	// no need : set is used
	// but if diff is too small -> problematic in codes
	// ==> erase intervals with length < eps
	double eraseEps = 1e-4;
	//vector<int> keepIdx;
	//keepIdx.push_back(0);
	{
		for (int i = ret._interval.size() - 2; i > 0; i--)
		{
			if (fabs(ret._interval[i] - ret._interval[i + 1]) < eraseEps)
			{
				//delete the one in front
				ret._interval.erase(ret._interval.begin() + i);
				ret._a.erase(ret._a.begin() + i);
				ret._b.erase(ret._b.begin() + i);
				ret._r.erase(ret._r.begin() + i);
				ret._arcIdx.erase(ret._arcIdx.begin() + i);
			}
		}
	}

	return ret;
}

/*
Def:
translate 180 degrees and flip sign
*/
void SupportComplex::flip()
{
	// 0. Alias
	double trans = PI;
	SupportComplex temp;

	// 1. Where is trans?
	int idx = 0;
	while (t()[idx + 1] < trans)
		idx++;

	// now: t()[i] < PI < t()[i+1]

	// 2. add
	temp.t().push_back(0.0);
	for (int i = idx + 1; i < t().size(); i++)
	{
		temp.a().push_back(a()[i - 1]);
		temp.b().push_back(b()[i - 1]);
		temp.r().push_back(r()[i - 1]);
		temp.arc().push_back(arc()[i - 1]);

		temp.t().push_back(t()[i] - trans);
	}
	for (int i = 0; i <= idx; i++)
	{
		temp.t().push_back(t()[i + 1] + PI2 - trans);
		
		temp.a().push_back(a()[i]);
		temp.b().push_back(b()[i]);
		temp.r().push_back(r()[i]);
		temp.arc().push_back(arc()[i]);
	}
	temp.t().back() = PI2; // temp.t().back() is wrong, fixing is done.

	// 3. simplify
	*this = _simplify(temp);

	// 4. flip a, b, r;
	//for (auto& a : a())
	//	a = -a;
	//for (auto& a : b())
	//	a = -a;
	for (auto& a : r())
		a = -a;


}

/*
Def: find t-value with intersection
Return:
	t \in [0, 2pi), ordered
*/
vector<double> SupportComplex::intersection(SupportComplex& rhs)
{
	vector<double> ret;

	// 0. Alias
	auto& lhs = *this;

	// 1. build interval
	set<double> interv;
	for (auto i : lhs._interval)
		interv.insert(i);
	for (auto i : rhs._interval)
		interv.insert(i);

	// 2. loop variables
	if (interv.size() == 0)
	{
		std::cerr << "!!! ERR in intersection() : interv.size() == 0" << endl;
	}
	auto it0 = interv.begin();
	auto it1 = it0;
	it1++;

	int idx0 = -1;
	int idx1 = -1;

	// 3. compare curves in interval
	// find intersections and btw intersections, compare size
	while (it1 != interv.end())
	{
		// i. interval is [t0, t1)
		auto t0 = *it0;
		auto t1 = *it1;

		// ii. set idx_
		// if(t0 is inside lhs.interval) idx0++;
		//if (find(lhs._interval.begin(), lhs._interval.end(), t0) != lhs._interval.end())
		//	idx0++;
		//if (find(rhs._interval.begin(), rhs._interval.end(), t0) != rhs._interval.end())
		//	idx1++;
		if (lhs._interval[idx0 + 1] == t0)
			idx0++;
		if (rhs._interval[idx1 + 1] == t0)
			idx1++;

		// flag for no intersection : if (true) simply pick a larger one
		bool noInter = false;

		// iii. find coeff
		double
			& a0 = lhs._a[idx0],
			& b0 = lhs._b[idx0],
			& r0 = lhs._r[idx0],
			& a1 = rhs._a[idx1],
			& b1 = rhs._b[idx1],
			& r1 = rhs._r[idx1];

		double
			da = a0 - a1,
			db = b0 - b1,
			dr = r0 - r1;

		// iv. solve da(cost) + db(sint) + dr = 0;
		// change it to cos(t - phi) = -dr/sqrt(...)
		double len = sqrt(da * da + db * db);
		double con = 2.0, conTheta;
		{
			if (len == 0.0)
				//only r is different => infinite intersection is though of as no inter
				noInter = true;
			else
				con = -dr / len;

			if (fabs(con) > 1.0)
				noInter = true;
			else
				conTheta = acos(con); // [0, pi]
		}

		if (!noInter)
		{
			double phi = atan2(db, da); // [-pi,pi]
			double
				theta0 = phi - conTheta, // [-2pi, pi]
				theta1 = phi + conTheta; // [-pi, 2pi]

			// turn them to [0, 2pi]
			if (theta0 < 0.0)
				theta0 += PI2;
			if (theta1 < 0.0)
				theta1 += PI2;

			// check whether they are valid and insert
			// bug-fixed : no guarantee that order is : 1.theta0 2. theta1
			if (theta0 < theta1)
			{
				if (t0 < theta0 && theta0 < t1)
					ret.push_back(theta0);
				if (t0 < theta1 && theta1 < t1)
					ret.push_back(theta1);
			}
			else
			{
				if (t0 < theta1 && theta1 < t1)
					ret.push_back(theta1);
				if (t0 < theta0 && theta0 < t1)
					ret.push_back(theta0);
			}
		}

		it0++;
		it1++;
	}

	return ret;
}

/*
Def:
	given t-value, find position in 2d space.
*/
Point SupportComplex::point(double t)
{
	Point ret;

	// confine t
	while (t < 0.0)
		t += PI2;
	while (t >= PI2)
		t -= PI2;

	// find interval
	int idx = 0;
	while (_interval[idx + 1] <= t)
		idx++;

	// find type
	int type;
	auto& arc = _arcIdx[idx][0];
	double
		t0, t1, t2, t3;
	{
		// find tX
		if (arc.lccw())
		{
			auto par = arc.param();
			t0 = par.first;
			t1 = par.second;
			t3 = t0 + PI2;
			t2 = (t3 + t1) * 0.5;
		}
		else
		{
			auto par = arc.param();
			t0 = par.second;
			t1 = par.first;
			t3 = t0 + PI2;
			t2 = (t3 + t1) * 0.5;
		}

		// confine t
		double theta = t;
		while (theta < t0)
			theta += PI2;
		while (theta > t3)
			theta -= PI2;

		if (theta < t1)
			type = 0;
		else if (theta < t2)
			type = 1;
		else
			type = 2;
	}

	// eval point
	switch (type)
	{
	case 2:
		return arc.cc() + arc.cr() * Point(cos(t0), sin(t0));
	case 1:
		return arc.cc() + arc.cr() * Point(cos(t1), sin(t1));
	case 0:
		return arc.cc() + arc.cr() * Point(cos(t), sin(t));
	}

}
