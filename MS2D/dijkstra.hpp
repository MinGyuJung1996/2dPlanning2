#pragma once
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include "xyt-cs.hpp"

namespace ms
{
	class Point;
}
using ms::Point;

namespace planning
{
	namespace output_to_file
	{
		struct v_edge;
	}
}

class xyt;
using std::vector;

namespace graphSearch
{

#pragma region Alias 

	/*
	lineSeg:
	API used:
		v0, v1, clr0, clr1
	*/
	using lseg = planning::output_to_file::v_edge;
	using namespace std;
	using namespace boost;
	//using ms::Point;

	using
		Weight =
		double;
	using
		Edge =
		std::pair<int, int>;
	using
		Graph =
		adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, Weight > >;
	using
		Vdesc =
		graph_traits < Graph >::vertex_descriptor;
	using
		Edesc =
		graph_traits < Graph >::edge_descriptor;

#pragma endregion

	enum EdgeType
	{
		etVor = 0, // voronoi
		etOff = 1, // offset curve
		etOcn = 2, // connection: vor - off
		etTan = 3, // common tangnet
		etTcn = 4, // connection : tan - (others)
		etIsc = 5  // InterSliceConnection
	};

	/*
	Def: vert used for dijkstra search
	*/
	struct Vert
	{
		inline double& x() {return _x;}
		inline double& y() {return _y;}
		inline double& z() {return _z;}
		
		//inline EdgeType& type() { return _t; }
		//bool operator<(const Vert& rhs) const;
	private:
		double _x, _y, _z;
		
		//EdgeType _t;
		//static double eqPtEps;
	};

	/*
	Not used
	*/
	struct vertData
	{

	};

	struct PointCompare
	{
	public:
		static double ptEqDist;
		bool operator()(const Point& lhs, const Point& rhs) const;
	};
}
namespace gs = graphSearch;

#define djkUseInterSliceBound true
#define djkInterSliceBound 0.01
//#define ptEqdist // currently static double of PointCompare
class djkCalc
{
	using Graph	 = gs::Graph;
	using Edge	 = gs::Edge;
	using Weight = gs::Weight;
	using lseg   = gs::lseg;

	using Vdesc  = gs::Vdesc;
	using Edesc  = gs::Edesc;

private:

	Graph			_graph;

	vector<Edge>	_edge;
	vector<Weight>	_weight;

	vector<xyt>		_vert;
	vector<double>  _clearance;
	vector<int>		_vertIdx; // [vertIdx[0], vertIdx[1]) : interval of idx of verts at slice 0

public:
	vector<vector<Point>>
		_vMat;
	vector<vector<Edge>>
		_eMat,
		_eMat2;
	vector<vector<Weight>>
		_wMat;

public:
	// con/destructor
	 djkCalc() = default;
	~djkCalc() = default;

	// member-access
	inline decltype(_graph)&
		getGraph() { return _graph; }
	inline decltype(_edge)&
		getEdge() { return _edge; }
	inline decltype(_weight)&
		getWeight() { return _weight; }
	inline decltype(_vert)&
		getVert() { return _vert; }
	inline decltype(_clearance)&
		getClr() { return _clearance; }

	void
	buildGraph
		(
			int _in_nSlice,
			vector<vector< lseg >>& _in_vor,		// in: voronoi Edges
			vector<vector< lseg >>& _in_off,		// in: offset
			vector<vector< lseg >>& _in_offCon,	// in: offset Connector
			vector<vector< lseg >>& _in_tan,		// in: comTan
			vector<vector< lseg >>& _in_tanCon,	// in: comTan connector
			Point& _in_robot_forward_direction
		);

	int
	searchGraph
		(
			xyt& _in_start,
			xyt& _in_end,
			vector<int>& _out_path
		);

private:
	//inline Weight wgtVor(Point& p, Point& q, Point& forward) { return (p - q).length1(); }
	//inline Weight wgtOff(Point& p, Point& q, Point& forward) { return (p - q).length1(); }
	//inline Weight wgtOcn(Point& p, Point& q, Point& forward) { return (p - q).length1(); }
	//inline Weight wgtTan(Point& p, Point& q, Point& forward) { return (p - q).length1(); }
	//inline Weight wgtTcn(Point& p, Point& q, Point& forward) { return (p - q).length1(); }
	//
	//using weightFunction = decltype(&wgtVor);
	//
	//void 
	//geometryToGraph
	//	(
	//	vector<lseg>&						_in_geo,
	//	//weightFunction						_in_wgtFunc,
	//	Point&								_in_forward,
	//	gs::EdgeType						_in_et,
	//
	//	map<Point, int, gs::PointCompare>&	_out_vm,
	//	vector<Point>&						_out_pl,
	//	vector<double>&						_out_cl,
	//	vector<Edge>&						_out_el,
	//	vector<Weight>&						_out_ew,
	//	vector<gs::EdgeType>&				_out_et
	//	);



};

using dijkstraCalculator = djkCalc;