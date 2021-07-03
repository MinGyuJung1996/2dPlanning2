#include "dijkstra.hpp"
#include "voronoi.hpp"
#include "xyt-cs.hpp"

#define COUT(x) #x " : " << x

namespace graphSearch
{


#pragma region PointCompare
	double PointCompare::ptEqDist = 1e-4;
	/*
	Def: comparator for points,
		but if they are in a box with length=ptEqDist, they are the same.
	Assume all points are in same slice
	*/
	bool PointCompare::operator()(const Point& lhs, const Point& rhs) const
	{

		if (lhs.xc() < rhs.xc() - ptEqDist)
			return true;
		if (lhs.xc() > rhs.xc() + ptEqDist)
			return false;

		if (lhs.yc() < rhs.yc() - ptEqDist)
			return true;
		else
			return false;
	}
#pragma endregion

	

	void test();
}

/*********************************************************************************************************************************************
**																																			**
**																Dijkstra Calculator															**
**																																			**
*********************************************************************************************************************************************/



/*
Def:
	build graph structure from geometries
Assume:
	Clearance data is built.
		Currently, all vertices will have their clr set. But not all edges. 
		(This means: at least one of the edges that contains a certain vertex will have clr value set, else edges have -1.0 as clr)
Desc:

Summary:

*/
void
djkCalc::buildGraph
	(
		int _in_nSlice,
		vector<vector<lseg>>& _in_vor,		// in: voronoi Edges
		vector<vector<lseg>>& _in_off,		// in: offset
		vector<vector<lseg>>& _in_offCon,	// in: offset Connector
		vector<vector<lseg>>& _in_tan,		// in: comTan
		vector<vector<lseg>>& _in_tanCon,	// in: comTan connector
		Point& _in_robot_forward_direction
	)
{
	using namespace graphSearch;

	// 0.  Alias
	auto& ns = _in_nSlice;

	auto& vor = _in_vor;
	auto& off = _in_off;
	auto& ocn = _in_offCon;
	auto& tan = _in_tan;
	auto& tcn = _in_tanCon;

	auto& frw = _in_robot_forward_direction;

	// 1. input processing
	vector<map<Point, int, PointCompare>> 
		viMaps		(ns);	// viMaps	  [sliceNo][pt]    -> ptIdx (start from 1, not zero)
	vector<vector<Point>>
		ptLists		(ns);	// ptList	  [sliceNo][ptIdx] -> pt
	vector<vector<double>>
		clrLists	(ns);	// clrList	  [sliceNo][ptIdx] -> clearance
	vector<vector<Edge>>
		edgeLists	(ns);	// edgeLists  [sliceNo][eIdx] -> edge
	vector<vector<Weight>>
		edgeWeights	(ns);	// edgeWeights[sliceNo][eIdx] -> edge weight
	vector<vector<EdgeType>>
		edgeTypes	(ns);	// edgeTypes  [sliceNo][eIdx] -> edge type
	{
		for (int sn = 0; sn < ns; sn++)
		{
			// 1-0. alias (for this slice)
			auto& vm = viMaps[sn];
			auto& pl = ptLists[sn];
			auto& cl = clrLists[sn];
			auto& el = edgeLists[sn];
			auto& ew = edgeWeights[sn];
			auto& et = edgeTypes[sn];

			// push dummy point -> now ptIdx corresponds to pt;
			pl.push_back(Point(0, 0));
			cl.push_back(-1.0);

			// !!!: later: use lambda

			// 1-1. voronoiEdges
			for (auto& e : vor[sn])
			{
				// i. vm & pl & cl (new)
				auto& idx0 = vm[e.v0];
				{
					// if(e.v0 not inside map) new idx
					if (idx0 == 0)
					{
						idx0 = vm.size();
						pl.push_back(e.v0);
						cl.push_back(e.clr0());
					}
				}

				auto& idx1 = vm[e.v1];
				{
					// if(e.v1 not inside map) new idx
					if (idx1 == 0)
					{
						idx1 = vm.size();
						pl.push_back(e.v1);
						cl.push_back(e.clr1());
					}
				}

				// ii. cl (update)
				{
					// if(cl not set)
					if (cl[idx0] < -1e-100)
						cl[idx0] = e.clr0();

					// if(cl not set)
					if (cl[idx1] < -1e-100)
						cl[idx1] = e.clr1();
				}

				if (idx0 != idx1)
				{
					// iii. el
					el.emplace_back(idx0, idx1);

					// iv. ew
					Weight w;
					{
						//
						w = (e.v0 - e.v1).length1();
					}
					ew.push_back(w);

					// v. et
					et.push_back(etVor);
				}
			}

			// 1-2. offset curves
			for (auto& e : off[sn])
			{
				// i. vm & pl & cl (new)
				auto& idx0 = vm[e.v0];
				{
					// if(e.v0 not inside map) new idx
					if (idx0 == 0)
					{
						idx0 = vm.size();
						pl.push_back(e.v0);
						cl.push_back(e.clr0());
					}
				}

				auto& idx1 = vm[e.v1];
				{
					// if(e.v1 not inside map) new idx
					if (idx1 == 0)
					{
						idx1 = vm.size();
						pl.push_back(e.v1);
						cl.push_back(e.clr1());
					}
				}

				// ii. cl (update)
				{
					// if(cl not set)
					if (cl[idx0] < -1e-100)
						cl[idx0] = e.clr0();

					// if(cl not set)
					if (cl[idx1] < -1e-100)
						cl[idx1] = e.clr1();
				}

				if (idx0 != idx1)
				{
					// iii. el
					el.emplace_back(idx0, idx1);

					// iv. ew
					Weight w;
					{
						//
						w = (e.v0 - e.v1).length1();
					}
					ew.push_back(w);

					// v. et
					et.push_back(etVor);
				}
			}

			// 1-3. offset connector
			for (auto& e : ocn[sn])
			{
				// i. vm & pl & cl (new)
				auto& idx0 = vm[e.v0];
				{
					// if(e.v0 not inside map) new idx
					if (idx0 == 0)
					{
						idx0 = vm.size();
						pl.push_back(e.v0);
						cl.push_back(e.clr0());
					}
				}

				auto& idx1 = vm[e.v1];
				{
					// if(e.v1 not inside map) new idx
					if (idx1 == 0)
					{
						idx1 = vm.size();
						pl.push_back(e.v1);
						cl.push_back(e.clr1());
					}
				}

				// ii. cl (update)
				{
					// if(cl not set)
					if (cl[idx0] < -1e-100)
						cl[idx0] = e.clr0();

					// if(cl not set)
					if (cl[idx1] < -1e-100)
						cl[idx1] = e.clr1();
				}

				if (idx0 != idx1)
				{
					// iii. el
					el.emplace_back(idx0, idx1);

					// iv. ew
					Weight w;
					{
						//
						w = (e.v0 - e.v1).length1();
					}
					ew.push_back(w);

					// v. et
					et.push_back(etVor);
				}
			}

			// 1-4. common tangent
			for (auto& e : tan[sn])
			{
				// i. vm & pl & cl (new)
				auto& idx0 = vm[e.v0];
				{
					// if(e.v0 not inside map) new idx
					if (idx0 == 0)
					{
						idx0 = vm.size();
						pl.push_back(e.v0);
						cl.push_back(e.clr0());
					}
				}

				auto& idx1 = vm[e.v1];
				{
					// if(e.v1 not inside map) new idx
					if (idx1 == 0)
					{
						idx1 = vm.size();
						pl.push_back(e.v1);
						cl.push_back(e.clr1());
					}
				}

				// ii. cl (update)
				{
					// if(cl not set)
					if (cl[idx0] < -1e-100)
						cl[idx0] = e.clr0();

					// if(cl not set)
					if (cl[idx1] < -1e-100)
						cl[idx1] = e.clr1();
				}

				if (idx0 != idx1)
				{
					// iii. el
					el.emplace_back(idx0, idx1);

					// iv. ew
					Weight w;
					{
						//
						w = (e.v0 - e.v1).length1();
					}
					ew.push_back(w);

					// v. et
					et.push_back(etVor);
				}
			}

			// 1-5. common tangent connectors
			for (auto& e : tcn[sn])
			{
				// i. vm & pl & cl (new)
				auto& idx0 = vm[e.v0];
				{
					// if(e.v0 not inside map) new idx
					if (idx0 == 0)
					{
						idx0 = vm.size();
						pl.push_back(e.v0);
						cl.push_back(e.clr0());
					}
				}

				auto& idx1 = vm[e.v1];
				{
					// if(e.v1 not inside map) new idx
					if (idx1 == 0)
					{
						idx1 = vm.size();
						pl.push_back(e.v1);
						cl.push_back(e.clr1());
					}
				}

				// ii. cl (update)
				{
					// if(cl not set)
					if (cl[idx0] < -1e-100)
						cl[idx0] = e.clr0();

					// if(cl not set)
					if (cl[idx1] < -1e-100)
						cl[idx1] = e.clr1();
				}

				if (idx0 != idx1)
				{
					// iii. el
					el.emplace_back(idx0, idx1);

					// iv. ew
					Weight w;
					{
						//
						w = (e.v0 - e.v1).length1();
					}
					ew.push_back(w);

					// v. et
					et.push_back(etVor);
				}
			}
		}
	}

	// 2. connect slices
	vector<vector<Edge>> 
		edgBtwSlc(ns); // edge between slices
	vector<vector<Weight>> 
		wgtBtwSlc(ns); // weight between slices
	// !!! TODO: Need optimization : Grid?
	{
		for (int sn = 0; sn < ns; sn++)
		{
			// 2-0. Alias
			auto snn = sn + 1;
			if (snn == ns)
				snn = 0;

			//cout << COUT(sn) << COUT(snn) << endl;

			auto& edg = edgBtwSlc[sn];
			auto& wgt = wgtBtwSlc[sn];

			auto& v0 = ptLists[sn];
			auto& v1 = ptLists[snn];
			auto& c0 = clrLists[sn];
			auto& c1 = clrLists[snn];

			// dbg
			int dbgVar = 1;

			// for(each pt pair)
			for(int i = 1; i < v0.size(); i++)
			for(int j = 1; j < v1.size(); j++)
			{
				auto d = (v0[i] - v1[j]).length1();
				if (djkUseInterSliceBound && d > djkInterSliceBound)
					continue;

				//if (two pts are inside clearance)
				if (d < c0[i] || d < c1[j])
				{
					edg.push_back(make_pair(i, j));
					
					Weight w = d;
					wgt.push_back(w);

					////dbg_out
					//if (dbgVar == 1)
					//{
					//	cout << COUT(sn) << endl;
					//	cout << COUT(i) << endl;
					//	cout << COUT(j) << endl;
					//	cout << COUT(v0[i]) << endl;
					//	cout << COUT(v1[j]) << endl;
					//
					//	dbgVar = 0;
					//}
				}

			}


		}
		
	}

	// 3. build index translation (needed for section 4)
	vector<int> idxTrn(ns); // index translation := change local idx in a slice to global idx // (sliceNo, vertNo) -> idxTrn[sliceNo] + vertNo;
	int totVrtCnt; // := number of all vertices
	{
		idxTrn[0] = -1; //since there is dummy node;

		for (int sn = 1; sn < ns; sn++)
		{
			idxTrn[sn] = idxTrn[sn - 1] + (ptLists[sn - 1].size() - 1);
			// idxTrn[sn] = (sum of number of enrties in [0, sn)) - 1;
		}

		totVrtCnt = idxTrn[ns - 1] + ptLists[ns - 1].size() - 1 + 1;
		
	}
	{
		// set _vertIdx
		_vertIdx.resize(ns + 1);
		for(int i = 0; i < _vertIdx.size(); i++)
		{
			_vertIdx[i] = idxTrn[i] + 1;
		}
		_vertIdx[_vertIdx.size()] = totVrtCnt;
	}
	
	// 4. Merge edges
	//map<Edge, int> edgeInvMap;
	vector<Edge>&	edge	= _edge;
	vector<Weight>&	weight	= _weight;
	//vector<EdgeType> etype;
	{
		edge.clear();
		weight.clear();

		// 4-1. single slice edges
		for (int sn = 0; sn < ns; sn++)
		{
			// Alias
			auto& el = edgeLists[sn];
			auto& ew = edgeWeights[sn];
			auto& et = edgeTypes[sn];

			// find translation
			auto& trn = idxTrn[sn];
			for (int i = 0; i < el.size(); i++)
			{
				edge.emplace_back(trn + el[i].first, trn + el[i].second);
				weight.push_back(ew[i]);
				//etype.push_back(et[i]);

				////dbg_out : e with idx -1
				//if (trn + el[i].first == -1)
				//{
				//	cout << "trn: " << trn << endl
				//		<< COUT(el[i].first) << endl;
				//}
				//if (trn + el[i].second == -1)
				//{
				//	cout << "trn: " << trn << endl
				//		<< COUT(el[i].second) << endl;
				//}
			}
		}

		// 4-2. inter-slice edges
		{
			for (int sn = 0; sn < edgBtwSlc.size(); sn++)
			{
				// Alias
				auto& el = edgBtwSlc[sn];
				auto& wl = wgtBtwSlc[sn];

				int snn = sn + 1;
				if (snn == edgBtwSlc.size())
					snn = 0;

				// find translation
				int trn0 = idxTrn[sn];
				int trn1 = idxTrn[snn];

				for (int i = 0; i < el.size(); i++)
				{
					edge.emplace_back(trn0 + el[i].first, trn1 + el[i].second);
					weight.push_back(wl[i]);

					////dbg_out : e with idx -1
					//if (trn0 + el[i].first == -1)
					//{
					//	cout << COUT(trn0) << endl
					//		<< COUT(el[i].first) << endl;
					//}
					//if (trn1 + el[i].second == -1)
					//{
					//	cout << COUT(trn1) << endl
					//		<< COUT(el[i].second) << endl;
					//}
				}
			}
		}
	}

	// 5. Merge verts
	vector<xyt>& vert = _vert;
	vector<double>& clearance = _clearance;
	{
		for (int sn = 0; sn < ns; sn++)
		{
			auto& vl = ptLists[sn];
			auto& cl = clrLists[sn];

			auto t = double(sn) / ns * 360.0;

			// for (all vertice, except dummy)
			for (int i = 1; i < vl.size(); i++)
			{
				xyt temp;
				temp.x() = vl[i].x();
				temp.y() = vl[i].y();
				temp.t() = t;

				vert.push_back(temp);
				clearance.push_back(cl[i]);
			}
		}
	}
	
	// 5.5 test validity
	{
		int minIdx = 1e8;
		int maxIdx = -1e8;
		for (auto& e : edge)
		{
			auto e0 = e.first;
			if (e0 < minIdx)
				minIdx = e0;
			if (e0 > maxIdx)
				maxIdx = e0;

			e0 = e.second;
			if (e0 < minIdx)
				minIdx = e0;
			if (e0 > maxIdx)
				maxIdx = e0;
		}

		
		cout << "dijkstra info" << endl;
		cout << "edge.size(): " << edge.size() << endl;
		cout << "minIdx: " << minIdx << endl;
		cout << "maxIdx: " << (maxIdx) << endl;
		cout << "vert.size()" << _vert.size() << endl;

		int sum = 0;
		for (auto& ptl : ptLists)
		{
			sum += ptl.size();
		}
		cout << "size sum : " << sum << endl;
		cout << "totvrtcnt: " << totVrtCnt << endl;

	}

	// 6. build Graph & swap;
	//Graph g(edge.begin(), edge.end(), &(weight.front()), totVrtCnt);
	Graph g(edge.begin(), edge.end(), weight.begin(), totVrtCnt);
	_graph.swap(g);

	// 7. swap (deubg)
	{
		_vMat.swap(ptLists);
		_eMat.swap(edgeLists);
		_eMat2.swap(edgBtwSlc);
		_wMat.swap(wgtBtwSlc);
	}


}


/*

Assume:
	start/end is each in some slice (that was calculated)
	start != end
Return:
	0 : no path
	1 : path exists
*/
int
djkCalc::searchGraph
	(
		xyt& _in_start,
		xyt& _in_end,
		vector<int>& _out_path
	)
{
	// alias
	auto& start = _in_start;
	auto& end	= _in_end;
	auto& path	= _out_path;

	int ns = _vertIdx.size() - 1;		// No. of slice
	int sz = start.t() / 360.0 * ns;	// sliceNo of start
	int ez = end.t() / 360.0 * ns;		// sliceNo of end

	// 1. find closest ver to start
	
	int vis, vie; //vertex index start/end 
	{
		xyt p = start;
		double minDist = 1.0e10;
		double minVert = -1;

		// for (all vert in same slice) find closest & clear vert
		for (int i = _vertIdx[sz]; i < _vertIdx[sz + 1]; i++)
		{
			auto& v = _vert[i];
			auto& c = _clearance[i];

			auto dx = p.x() - v.x();
			auto dy = p.y() - v.y();
			auto dist = sqrt(dx * dx + dy * dy);

			// if(p is inside clearance of v && dist updatable)
			if (dist < c && dist < minDist)
			{
				minVert = i;
				minDist = dist;
			}
		}

		if (minVert == -1)
			return false;
		else
			vis = minVert;
	}
	{
		xyt p = end;
		double minDist = 1.0e+10;
		double minVert = -1;

		// for (all vert in same slice) find closest & clear vert
		for (int i = _vertIdx[ez]; i < _vertIdx[ez + 1]; i++)
		{
			auto& v = _vert[i];
			auto& c = _clearance[i];

			auto dx = p.x() - v.x();
			auto dy = p.y() - v.y();
			auto dist = sqrt(dx * dx + dy * dy);

			// if(p is inside clearance of v && dist updatable)
			if (dist < c && dist < minDist)
			{
				minVert = i;
				minDist = dist;
			}
		}

		if (minVert == -1)
			return false;
		else
			vie = minVert;
	}

	// 2. do graph search
	Vdesc			sv = boost::vertex(vis, _graph);
	vector<Weight>	dist(boost::num_vertices(_graph));
	vector<Vdesc>	pred(boost::num_vertices(_graph));

	boost::dijkstra_shortest_paths(_graph, sv, boost::predecessor_map(&pred[0]).distance_map(&dist[0]));

	// 3. return 
	//if(end point not reachable)
	if (pred[vie] == vie)
		return false;

	// follow pred and build path
	int cur = vie;
	while (cur != vis)
	{
		path.push_back(cur);
		cur = pred[cur];
	}
	path.push_back(vis);

	reverse(path.begin(), path.end());

	//debug
	if(1)
	{
		using namespace boost;

		double d0 = dist[vie];
		double d1 = 0.0;
		for (int i = 0; i < path.size() - 1; i++)
		{
			int in = i + 1;
			int idx0 = path[i];  //vertIdx;
			int idx1 = path[in]; //vertIdx;

			// weight in graph
			Vdesc vd0 = vertex(idx0, _graph);
			Vdesc vd1 = vertex(idx1, _graph);
			auto edt = edge(vd0, vd1, _graph);
			Edesc ed = edt.first;
			auto weight_g = get(edge_weight, _graph, ed);

			// weight before making graph
			int eidxOrig = -1;
			for (int i = 0; i < _edge.size(); i++)
			{
				if (_edge[i].first == idx0 && _edge[i].second == idx1)
				{
					eidxOrig = i;
					break;
				}
				if (_edge[i].first == idx1 && _edge[i].second == idx0)
				{
					eidxOrig = i;
					break;
				}
			}
			auto weight_orig = _weight[eidxOrig];

			auto& v0 = _vert[idx0];
			auto& v1 = _vert[idx1];
			auto dv = v0 - v1;
			
			auto weight_calc = sqrt((dv.x() * dv.x()) + (dv.y() * dv.y()));
			d1 += weight_calc;

			cout << i << ": wo wg wc : " << weight_orig << "   " << weight_g << " " << weight_calc << endl;
			// conc : weight is wrongly set at build graph
		}

		// _weight (right before graph)
		// boost (current graph)
		// sqrt(dP)
		// 

		cout << COUT(d0) << endl;
		cout << COUT(d1) << endl;
	}

	return true;
}

//void djkCalc::geometryToGraph
//(
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
//)
//{
//	// 0. Alias
//	auto& geo = _in_geo;
//
//	auto& vm = _out_vm;
//	auto& pl = _out_pl;
//	auto& cl = _out_cl;
//	auto& el = _out_el;
//	auto& ew = _out_ew;
//	auto& et = _out_et;
//
//	for (auto& e : geo)
//	{
//		// i. vm & pl & cl (new)
//		auto& idx0 = vm[e.v0];
//		{
//			// if(e.v0 not inside map) new idx
//			if (idx0 == 0)
//			{
//				idx0 = vm.size();
//				pl.push_back(e.v0);
//				cl.push_back(e.clr0());
//			}
//		}
//
//		auto& idx1 = vm[e.v1];
//		{
//			// if(e.v1 not inside map) new idx
//			if (idx1 == 0)
//			{
//				idx1 = vm.size();
//				pl.push_back(e.v1);
//				cl.push_back(e.clr1());
//			}
//		}
//
//		// ii. cl (update)
//		{
//			// if(cl not set)
//			if (cl[idx0] < -1e-100)
//				cl[idx0] = e.clr0();
//
//			// if(cl not set)
//			if (cl[idx1] < -1e-100)
//				cl[idx1] = e.clr1();
//		}
//
//		if (idx0 != idx1)
//		{
//			// iii. el
//			el.emplace_back(idx0, idx1);
//
//			// iv. ew
//			Weight w = _in_wgtFunc(e.v0, e.v1, _in_forward);
//			ew.push_back(w);
//
//			// v. et
//			et.push_back(_in_et);
//		}
//	}
//}
