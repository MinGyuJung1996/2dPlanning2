#pragma once

/*
******************************************************************************************************************************************************
**																																					**
**															Support Functions																		**
**																																					**
******************************************************************************************************************************************************

A support given point p is a function of theta
	f(theta) := dot((cosTheta, sinTheta), p)
A suuport given set P is a funciton of theta
	F(theta) := max{p \in P} f(theta)
Here, the sets are circular arcs or a set of circular arcs.

Set of points that 

*/

#include <vector>
#include <set>
#include <iostream>
#include <fstream>

#pragma region Class/Funcs from other files
namespace ms
{
	class CircularArc;
	class Point;
}
//using namespace ms;
using ms::CircularArc;
using ms::Point;
using std::vector;
using std::set;


//template<class T>
//using vector = std::vector<T>;
//template<class T>
//using set = std::set<T>;

using std::cout;
using std::cerr;
using std::endl;
using std::find;

#pragma endregion

/*
Notes About Support-function of circular arc
	Given 
		Arc(x,y,r,t0,t1) (t1 > t0)
		t2 = (t1 + (t0+2pi)) / 2	
		t3 = t0 + 2pi
	support(arc, theta) :=
		if(t0 < theta < t1) x*cos(t) + y*sin(t) + r 
		if(t1 < theta < t2) (x+r*cos(t1), y+r*sin(t1)) * (cost, sint)
		if(t2 < theta < t3) (x+r*cos(t0), y+r*sin(t0)) * (cost, sint)
		Above funcs are C0

Notice that if one model is completely inside another, no commonTangent
*/


/*****************************************************************************************************************************************************
**																																					**
**															Support Class																			**
**																																					**
*****************************************************************************************************************************************************/

/*
Def: Support of an arc

API:
	construct
	evaluate: theta --> support-val
*/
struct Support
{
public:
	// Const/Dest
	Support() = default;
	~Support() = default;
	Support(CircularArc& arc);
	
	// funcs
	void 
		build(CircularArc& arc);
	double
		eval(double theta, double c = 2.0, double s = 2.0);
	int
		type(double theta);
	double
		f0(double c, double s),
		f1(double c, double s),
		f2(double c, double s);

	void
		print();
	//double
	//	eval(double cosTheta, double sinTheta);
	
	// accessor
	inline CircularArc& ori() { return *_ori; }
	inline double& t0(){ return _t0; }
	inline double& t1(){ return _t1; }
	inline double& t2(){ return _t2; }
	inline double& t3(){ return _t3; }
	inline double& x0(){ return _x0; }
	inline double& y0(){ return _y0; }
	inline double& r0(){ return _r0; }
	inline double& x1(){ return _x1; }
	inline double& y1(){ return _y1; }
	inline double& x2(){ return _x2; }
	inline double& y2(){ return _y2; }

private:
	CircularArc* _ori;
	double
		_t0, _t1, _t2, _t3;
	double
		_x0, _y0, _r0,
		_x1, _y1,
		_x2, _y2;
	int 
		_flipped;
};
using Sup = Support;


/*****************************************************************************************************************************************************
**																																					**
**															Common Tangent Calculator																**
**																																					**
*****************************************************************************************************************************************************/
struct CommonTangentCalculator
{
public:
	static vector<Point> eval(Sup& a, Sup& b); //TODO
	static vector<Point> eval(CircularArc& a, CircularArc& b); //TODO

	void setInput0(vector<CircularArc>& loop0);
	void setInput1(vector<CircularArc>& loop1);
	void setOutput(vector<Point>& out);
	void calc();

private:
	vector<CircularArc>
		* _loop0,
		* _loop1;
	vector<Point>
		* _out;
};
using ComTanCalc = CommonTangentCalculator;


/*****************************************************************************************************************************************************
**																																					**
**															  Convex Hull Calculator																**
**																																					**
*****************************************************************************************************************************************************/

/*
Def:
	Convex Hull given lots of loops.
Assume:
	The order of the circular arcs doesn't matter, but they should be part of some loop.
	Sampling is dense enough
	??? G1 of model. ? G0 of model
*/
struct ConvexHullCalculator
{
public:
	void setInput(vector<CircularArc>& inputArcs);
	void setOutput(vector<CircularArc>& outputArcs);
	void calcSam(); // TODO
	void calcAlg();

	static void offset(vector<CircularArc>& in, vector<CircularArc>& out, double offRad);

	static inline double& lineRad() { return _lineRad; }

private:
	static double _lineRad;
	static int _nSam;

	vector<CircularArc>* _iArcs;
	vector<CircularArc>* _oArcs;
};
using CvxHullCalc = ConvexHullCalculator;

/*
Def: Carries pieces of g0-continuous {a*cos(theta) + b*sin(theta) + r)

Member:
	interval : size() = (No. of pieces + 1)
		how [0, 2pi) is divided;

*/
struct SupportComplex
{
public:
	SupportComplex() = default;
	~SupportComplex() = default;
	SupportComplex(Support& sup);
	
	void build(Support& sup);
	void draw();
	void print();
	
	void flip();

	SupportComplex upperEnvelope(SupportComplex& rhs);
	vector<double> intersection(SupportComplex& rhs);
	Point point(double t);


	inline vector<double>& a() { return _a; }
	inline vector<double>& b() { return _b; }
	inline vector<double>& r() { return _r; }
	inline vector<double>& t() { return _interval; }
	inline vector<CircularArc*>& arc() { return _arcIdx; }
private:
	static SupportComplex _simplify(SupportComplex& temp);

	vector<double>
		_a,				// coeff
		_b,				// coeff
		_r;				// coeff
	vector<double>
		_interval;		// interval endpoints, where curve equation changes. (always contain 0, 2pi at each end)
	vector<CircularArc*>
		_arcIdx;		// arc ptr, where corresponding interval comes from

	double _max, _min;

	static int _drawSampleNo;
};
using SupCom = SupportComplex;
