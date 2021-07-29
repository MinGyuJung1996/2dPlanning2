#include "dijkstra.hpp"
#include "voronoi.hpp"
#include "xyt-cs.hpp"

#define COUT(x) #x " : " << x
#define ONLY_USE_VORONOI_IN_GRAPH true

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

			auto degree = double(sn) / ns * 360.0;

			// push dummy point -> now ptIdx corresponds to pt;
			pl.push_back(Point(0, 0));
			cl.push_back(-1.0);

			// !!!: later: use lambda

			// 1-1. voronoiEdges
			for (auto& e : vor[sn])
			{
				//// i. vm & pl & cl (new)
				//auto& idx0 = vm[e.v0];
				//{
				//	// if(e.v0 not inside map) new idx
				//	if (idx0 == 0)
				//	{
				//		idx0 = vm.size();
				//		pl.push_back(e.v0);
				//		cl.push_back(e.clr0());
				//	}
				//}
				//
				//auto& idx1 = vm[e.v1];
				//{
				//	// if(e.v1 not inside map) new idx
				//	if (idx1 == 0)
				//	{
				//		idx1 = vm.size();
				//		pl.push_back(e.v1);
				//		cl.push_back(e.clr1());
				//	}
				//}

				// i. vm & pl & cl (new)
				auto& pair0 = vm.insert({ e.v0, 0 });
				auto& v0 = pair0.first->first;
				auto& idx0 = pair0.first->second;
				{
					// if(e.v0 not inside map) new idx
					if (idx0 == 0)
					{
						idx0 = vm.size();
						pl.push_back(e.v0);
						cl.push_back(e.clr0());
					}
				}

				auto& pair1 = vm.insert({ e.v1, 0 });
				auto& v1 = pair1.first->first;
				auto& idx1 = pair1.first->second;
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
					Weight w = findWeight(v0, v1, frw, degree);
					//Weight w = findWeight(e.v0, e.v1, frw);
					//Weight w;
					//{
					//	//
					//	w = (e.v0 - e.v1).length1();
					//}
					ew.push_back(w);

					// v. et
					et.push_back(etVor);
				}
			}

			// 1-2. offset curves
			if (!ONLY_USE_VORONOI_IN_GRAPH)
			for (auto& e : off[sn])
			{
				//// i. vm & pl & cl (new)
				//auto& idx0 = vm[e.v0];
				//{
				//	// if(e.v0 not inside map) new idx
				//	if (idx0 == 0)
				//	{
				//		idx0 = vm.size();
				//		pl.push_back(e.v0);
				//		cl.push_back(e.clr0());
				//	}
				//}
				//
				//auto& idx1 = vm[e.v1];
				//{
				//	// if(e.v1 not inside map) new idx
				//	if (idx1 == 0)
				//	{
				//		idx1 = vm.size();
				//		pl.push_back(e.v1);
				//		cl.push_back(e.clr1());
				//	}
				//}

				// i. vm & pl & cl (new)
				auto& pair0 = vm.insert({ e.v0, 0 });
				auto& v0 = pair0.first->first;
				auto& idx0 = pair0.first->second;
				{
					// if(e.v0 not inside map) new idx
					if (idx0 == 0)
					{
						idx0 = vm.size();
						pl.push_back(e.v0);
						cl.push_back(e.clr0());
					}
				}

				auto& pair1 = vm.insert({ e.v1, 0 });
				auto& v1 = pair1.first->first;
				auto& idx1 = pair1.first->second;
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
					Weight w = findWeight(v0, v1, frw, degree);
					//Weight w = findWeight(e.v0, e.v1, frw);
					//Weight w;
					//{
					//	//
					//	w = (e.v0 - e.v1).length1();
					//}
					ew.push_back(w);

					// v. et
					et.push_back(etOff);
				}
			}

			// 1-3. offset connector
			if (!ONLY_USE_VORONOI_IN_GRAPH)
			for (auto& e : ocn[sn])
			{
				//// i. vm & pl & cl (new)
				//auto& idx0 = vm[e.v0];
				//{
				//	// if(e.v0 not inside map) new idx
				//	if (idx0 == 0)
				//	{
				//		idx0 = vm.size();
				//		pl.push_back(e.v0);
				//		cl.push_back(e.clr0());
				//	}
				//}
				//
				//auto& idx1 = vm[e.v1];
				//{
				//	// if(e.v1 not inside map) new idx
				//	if (idx1 == 0)
				//	{
				//		idx1 = vm.size();
				//		pl.push_back(e.v1);
				//		cl.push_back(e.clr1());
				//	}
				//}

				// i. vm & pl & cl (new)
				auto& pair0 = vm.insert({ e.v0, 0 });
				auto& v0 = pair0.first->first;
				auto& idx0 = pair0.first->second;
				{
					// if(e.v0 not inside map) new idx
					if (idx0 == 0)
					{
						idx0 = vm.size();
						pl.push_back(e.v0);
						cl.push_back(e.clr0());
					}
				}

				auto& pair1 = vm.insert({ e.v1, 0 });
				auto& v1 = pair1.first->first;
				auto& idx1 = pair1.first->second;
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
					Weight w = findWeight(v0, v1, frw, degree);
					//Weight w = findWeight(e.v0, e.v1, frw);
					//Weight w;
					//{
					//	//
					//	w = (e.v0 - e.v1).length1();
					//}
					ew.push_back(w);

					// v. et
					et.push_back(etOcn);
				}
			}

			// 1-4. common tangent
			if(!ONLY_USE_VORONOI_IN_GRAPH)
			for (auto& e : tan[sn])
			{
				//// i. vm & pl & cl (new)
				//auto& idx0 = vm[e.v0];
				//{
				//	// if(e.v0 not inside map) new idx
				//	if (idx0 == 0)
				//	{
				//		idx0 = vm.size();
				//		pl.push_back(e.v0);
				//		cl.push_back(e.clr0());
				//	}
				//}
				//
				//auto& idx1 = vm[e.v1];
				//{
				//	// if(e.v1 not inside map) new idx
				//	if (idx1 == 0)
				//	{
				//		idx1 = vm.size();
				//		pl.push_back(e.v1);
				//		cl.push_back(e.clr1());
				//	}
				//}

				// i. vm & pl & cl (new)
				auto& pair0 = vm.insert({ e.v0, 0 });
				auto& v0 = pair0.first->first;
				auto& idx0 = pair0.first->second;
				{
					// if(e.v0 not inside map) new idx
					if (idx0 == 0)
					{
						idx0 = vm.size();
						pl.push_back(e.v0);
						cl.push_back(e.clr0());
					}
				}

				auto& pair1 = vm.insert({ e.v1, 0 });
				auto& v1 = pair1.first->first;
				auto& idx1 = pair1.first->second;
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
					Weight w = findWeight(v0, v1, frw, degree);
					//Weight w = findWeight(e.v0, e.v1, frw);
					//Weight w;
					//{
					//	//
					//	w = (e.v0 - e.v1).length1();
					//}
					ew.push_back(w);

					// v. et
					et.push_back(etTan);
				}
			}

			// 1-5. common tangent connectors
			if(!ONLY_USE_VORONOI_IN_GRAPH)
			for (auto& e : tcn[sn])
			{
				//// i. vm & pl & cl (new)
				//auto& idx0 = vm[e.v0];
				//{
				//	// if(e.v0 not inside map) new idx
				//	if (idx0 == 0)
				//	{
				//		idx0 = vm.size();
				//		pl.push_back(e.v0);
				//		cl.push_back(e.clr0());
				//	}
				//}
				//
				//auto& idx1 = vm[e.v1];
				//{
				//	// if(e.v1 not inside map) new idx
				//	if (idx1 == 0)
				//	{
				//		idx1 = vm.size();
				//		pl.push_back(e.v1);
				//		cl.push_back(e.clr1());
				//	}
				//}

				// i. vm & pl & cl (new)
				auto& pair0 = vm.insert({ e.v0, 0 });
				auto& v0 = pair0.first->first;
				auto& idx0 = pair0.first->second;
				{
					// if(e.v0 not inside map) new idx
					if (idx0 == 0)
					{
						idx0 = vm.size();
						pl.push_back(e.v0);
						cl.push_back(e.clr0());
					}
				}

				auto& pair1 = vm.insert({ e.v1, 0 });
				auto& v1 = pair1.first->first;
				auto& idx1 = pair1.first->second;
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
					Weight w = findWeight(v0, v1, frw, degree);
					//Weight w = findWeight(e.v0, e.v1, frw);
					//Weight w;
					//{
					//	//
					//	w = (e.v0 - e.v1).length1();
					//}
					ew.push_back(w);

					// v. et
					et.push_back(etTcn);
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

			auto degree = double(2 * sn + 1) / 2.0 / ns * 360.0;

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
					
					//Weight w = d;
					Weight w = findWeight(v0[i], v1[j], frw, degree) * djkInterSliceWeightMultiplier;
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


			////dbg_out
			//std::cout << "Connect btw: " << v0.size() << " x " << v1.size() << " = " << edg.size() << endl;
		}
		
	}

	// 3. build index translation (needed for section 4)
	vector<int>& idxTrn = _idxTrn; // index translation := change local idx in a slice to global idx // (sliceNo, vertNo) -> idxTrn[sliceNo] + vertNo;
	idxTrn.resize(ns);
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
	vector<EdgeType>& etype = _etype;
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
				etype.push_back(et[i]);

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
					etype.push_back(gs::etIsc);

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
	if(false)
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

	// 7. build invEdge;
	{
		for(int i = 0; i < edge.size(); i++)
		{
			_invEdge[edge[i]] = i;
		}

		cout << COUT(edge.size()) << endl;
		cout << COUT(_invEdge.size()) << endl;
	}

	// 99. swap (deubg)
	{
		_vMat.swap(ptLists);
		_cMat.swap(clrLists);

		_eMat.swap(edgeLists);
		_eMat2.swap(edgBtwSlc);
		_wMat.swap(wgtBtwSlc);
	}


}


/*
Def:
	1. finds closest point to start/end in graph
	2. do dijk from those points.
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

	//debug dbg_out
	if(1)
	{
		using namespace boost;

		cout << COUT(_edge.size()) << endl;
		cout << COUT(_weight.size()) << endl;
		cout << COUT(_etype.size()) << endl;

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

			bool allEQ = false;
			if (abs(weight_orig - weight_g) < 1e-4 && abs(weight_orig - weight_calc) < 1e-4)
				allEQ = true;
			cout << i << ": wo wg wc : "<<allEQ <<" type:" << _etype[eidxOrig] << "   " << weight_orig << "   " << weight_g << " " << weight_calc << endl;
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

#include "collision detection.hpp"
namespace graphSearch
{
	extern std::vector<cd::lineSegmentCollisionTester> testers;
	extern std::vector<cd::pointCollisionTester> testers2;
}
/*
Def:
	similar to searchGraph, but the robot rotates near start/end
Assume:
	start/end is each in some slice (that was calculated)
	start != end
Return:
	0 : no path
	1 : path exists
*/
int
djkCalc::searchGraph2
(
	xyt& _in_start,
	xyt& _in_end,
	vector<xyt>& _out_path,
	vector<double>& _out_clr,
	vector<gs::EdgeType>& _out_et,
	Point& _in_robotForward
)
{
	auto stamp0 = clock();

	// hyper
	double distBound = 0.3;
	double distBoundSq = distBound * distBound;
	constexpr bool alterDjkCalc = false; //TODO
	constexpr bool useLineTest = false;
	auto cpyGraph = _graph;

	// alias
	auto& start = _in_start;
	auto& end = _in_end;
	//auto& path = _out_path;

	int ns = _vertIdx.size() - 1;		// No. of slice
	int sz = start.t() / 360.0 * ns;	// sliceNo of start
	int ez = end.t() / 360.0 * ns;		// sliceNo of end

	// 1. find rotatable angle
	auto stamp1 = clock();
	vector<int> rotatableS, rotatableE; // idx [0, number_of_slice)
	vector<double> rotatableSclr, rotatableEclr;
	int sRotIdx, eRotIdx;
	{
		Point p(start.x(), start.y());
		double z = sz;
		auto& res = rotatableS;
		auto& res2 = rotatableSclr;
		auto& rotIdx = sRotIdx;

		// 1-1. test orig
		Point closest;
		if (gs::testers2[z].testPrecise(p, closest))
			return false;
		auto clr0 = (p - closest).length1();

		// 1-2. rising index;
		vector<int> temp;
		vector<double> tempClr;
		for (int i = z + 1; i < z + ns; i++)
		{
			auto idx = i;
			if (i >= ns)
				idx -= ns;
			if (gs::testers2[idx].testPrecise(p, closest))
				break;

			temp.push_back(i);
			tempClr.push_back((p - closest).length1());
		}

		// 1-3. decreasing
		vector<int> temp2;
		vector<double> tempClr2;
		for (int i = z - 1; i > z + temp.size() - ns; i--)
		{
			auto idx = i;
			if (i < 0)
				idx += ns;
			if (gs::testers2[idx].testPrecise(p, closest))
				break;

			temp2.push_back(i);
			tempClr2.push_back((p - closest).length1());
		}

		// build res
		res.insert(res.end(), temp2.rbegin(), temp2.rend());
		res.push_back(z);
		res.insert(res.end(), temp.begin(), temp.end());

		res2.insert(res2.end(), tempClr2.rbegin(), tempClr2.rend());
		res2.push_back(clr0);
		res2.insert(res2.end(), tempClr.begin(), tempClr.end());

		rotIdx = temp2.size();
	}
	{
		// 1-0. changed part
		Point p(end.x(), end.y());
		double z = ez;
		auto& res = rotatableE;
		auto& res2 = rotatableEclr;
		auto& rotIdx = eRotIdx;

		// 1-1. test orig
		Point closest;
		if (gs::testers2[z].testPrecise(p, closest))
			return false;
		auto clr0 = (p - closest).length1();

		// 1-2. rising index;
		vector<int> temp;
		vector<double> tempClr;
		for (int i = z + 1; i < z + ns; i++)
		{
			auto idx = i;
			if (i >= ns)
				idx -= ns;
			if (gs::testers2[idx].testPrecise(p, closest))
				break;

			temp.push_back(i);
			tempClr.push_back((p - closest).length1());
		}

		// 1-3. decreasing
		vector<int> temp2;
		vector<double> tempClr2;
		for (int i = z - 1; i > z + temp.size() - ns; i--)
		{
			auto idx = i;
			if (i < 0)
				idx += ns;
			if (gs::testers2[idx].testPrecise(p, closest))
				break;

			temp2.push_back(i);
			tempClr2.push_back((p - closest).length1());
		}

		// build res
		res.insert(res.end(), temp2.rbegin(), temp2.rend());
		res.push_back(z);
		res.insert(res.end(), temp.begin(), temp.end());

		res2.insert(res2.end(), tempClr2.rbegin(), tempClr2.rend());
		res2.push_back(clr0);
		res2.insert(res2.end(), tempClr.begin(), tempClr.end());

		rotIdx = temp2.size();
	}

	// 2. find connectable
	auto stamp2 = clock();
	vector<vector<int>> connectableS(rotatableS.size()), connectableE(rotatableE.size());
	vector<vector<Weight>> conWeightS(rotatableS.size()), conWeightE(rotatableE.size());
	{
		auto& r = rotatableS;
		auto& clrs = rotatableSclr;
		auto& res = connectableS;
		auto& res2 = conWeightS;
		Point p(start.x(), start.y());
		
		for (int i = 0; i < r.size(); i++)
		{
			// alias
			auto& c = res[i];
			auto& clr0 = clrs[i];
			auto& w = res2[i];
			auto& idx = r[i];
			if (idx < 0)
				idx += ns;
			if (idx >= ns)
				idx -= ns;

			auto trn = _idxTrn[idx];
			auto deg = double(idx) / ns * 360.0;
			
			for (int j = 1; j < _vMat[idx].size(); j++)
			{
				auto& v = _vMat[idx][j];

				auto dsq = (p - v).length2();
				if (useLineTest)
				{
					if (dsq < distBoundSq)
					{
						if (gs::testers[idx].test(p, v))
							continue;
						else
						{
							c.push_back(trn + j);
							w.push_back(findWeight(p, v, _in_robotForward, deg));
						}

					}
				}
				else
				{
					auto d = sqrt(dsq);
					auto clr1 = _cMat[idx][j];

					if (d < clr1 || d < clr0)
					{
						c.push_back(trn + j);
						w.push_back(findWeight(p, v, _in_robotForward, deg));
					}
				}
			}
		}
	}
	{
		auto& r = rotatableE;
		auto& clrs = rotatableEclr;
		auto& res = connectableE;
		auto& res2 = conWeightE;
		Point p(end.x(), end.y());

		for (int i = 0; i < r.size(); i++)
		{
			// alias
			auto& c = res[i];
			auto& clr0 = clrs[i];
			auto& w = res2[i];
			auto& idx = r[i];
			if (idx < 0)
				idx += ns;
			if (idx >= ns)
				idx -= ns;

			auto trn = _idxTrn[idx];
			auto deg = double(idx) / ns * 360.0;

			for (int j = 1; j < _vMat[idx].size(); j++)
			{
				auto& v = _vMat[idx][j];

				auto dsq = (p - v).length2();
				if (useLineTest)
				{
					if (dsq < distBoundSq)
					{
						if (gs::testers[idx].test(p, v))
							continue;
						else
						{
							c.push_back(trn + j);
							w.push_back(findWeight(p, v, _in_robotForward, deg));
						}

					}
				}
				else
				{
					auto d = sqrt(dsq);
					auto clr1 = _cMat[idx][j];

					if (d < clr1 || d < clr0)
					{
						c.push_back(trn + j);
						w.push_back(findWeight(p, v, _in_robotForward, deg));
					}
				}
			}
		}
	}

	// 3. add vertices
	auto stamp3 = clock();
	int sTrn = boost::num_vertices(_graph);
	vector<Vdesc> svd;
	{
		for (int i = 0; i < rotatableS.size(); i++)
		{
			svd.push_back(add_vertex(_graph));

			//change original graph
			{
				xyt temp = start;
				auto deg = double(rotatableS[i]) / ns * 360.0;
				if (deg < 0.0)
					deg += 360.0;
				if (deg >= 360.0)
					deg -= 360.0;
				temp.t() = deg;

				if (alterDjkCalc)
				{
					_vert.push_back(temp);
					_clearance.push_back(rotatableSclr[i]);
				}
			}
		}
	}
	int eTrn = boost::num_vertices(_graph);
	vector<Vdesc> evd;
	{
		for (int i = 0; i < rotatableE.size(); i++)
		{
			evd.push_back(add_vertex(_graph));

			//change original graph
			{
				xyt temp = end;
				auto deg = double(rotatableE[i]) / ns * 360.0;
				if (deg < 0.0)
					deg += 360.0;
				if (deg >= 360.0)
					deg -= 360.0;
				temp.t() = deg;

				if (alterDjkCalc)
				{
					_vert.push_back(temp);
					_clearance.push_back(rotatableEclr[i]);
				}
			}
		}
	}

	// 4. add edges (horizontal)
	auto stamp4 = clock();
	vector<int> eIdxList;
	eIdxList.reserve(5);
	eIdxList.push_back(boost::num_edges(_graph));
	vector<Edesc> sed;
	{
		auto& c = connectableS;
		auto& w = conWeightS;
		auto& edv = sed;
		auto& vdv = svd;

		for (int i = 0; i < c.size(); i++)
		{
			for (int j = 0; j < c[i].size(); j++)
			{
				// boost
				auto v = boost::vertex(c[i][j], _graph);
				//boost::property<boost::edge_weight_t, Weight > prop()
				auto ed = add_edge(vdv[i], v, _graph).first;
				boost::put(boost::edge_weight, _graph, ed, w[i][j]);
				edv.push_back(ed);

				// djkCalc

				// TODO
			}
		}
	}
	eIdxList.push_back(boost::num_edges(_graph));
	vector<Edesc> eed;
	{
		auto& c = connectableE;
		auto& w = conWeightE;
		auto& edv = eed;
		auto& vdv = evd;

		for (int i = 0; i < c.size(); i++)
		{
			for (int j = 0; j < c[i].size(); j++)
			{
				auto v = boost::vertex(c[i][j], _graph);
				//boost::property<boost::edge_weight_t, Weight > prop()
				auto ed = add_edge(vdv[i], v, _graph).first;
				boost::put(boost::edge_weight, _graph, ed, w[i][j]);
				edv.push_back(ed);
			}
		}

	}
	eIdxList.push_back(boost::num_edges(_graph));

	// 5. add edges (vertical)
	auto stamp5 = clock();
	{
		auto& vdv = svd;
		auto& edv = sed;

		for (int i = 0; i < vdv.size() - 1; i++)
		{
			auto ed = add_edge(vdv[i], vdv[i+1], _graph).first;
			boost::put(boost::edge_weight, _graph, ed, 1e-12);
			edv.push_back(ed);
		}

		// if(robot can rotate full 360 degree) connect last-first vert
		if (vdv.size() == ns)
		{
			cout << "robot can rotate fully at start" << endl;
			auto ed = add_edge(vdv.front(), vdv.back(), _graph).first;
			boost::put(boost::edge_weight, _graph, ed, 1e-12);
			edv.push_back(ed);
		}
	}
	eIdxList.push_back(boost::num_edges(_graph));
	{
		auto& vdv = evd;
		auto& edv = eed;

		for (int i = 0; i < vdv.size() - 1; i++)
		{
			auto ed = add_edge(vdv[i], vdv[i + 1], _graph).first;
			boost::put(boost::edge_weight, _graph, ed, 1e-10);
			edv.push_back(ed);
		}

		// if(robot can rotate full 360 degree) connect last-first vert
		if (vdv.size() == ns)
		{
			cout << "robot can rotate fully at end" << endl;
			auto ed = add_edge(vdv.front(), vdv.back(), _graph).first;
			boost::put(boost::edge_weight, _graph, ed, 1e-12);
			edv.push_back(ed);
		}
	}
	eIdxList.push_back(boost::num_edges(_graph));


	// 6. do graph search
	auto stamp6 = clock();
	Vdesc			sv = svd[sRotIdx];
	vector<Weight>	dist(boost::num_vertices(_graph));
	vector<Vdesc>	pred(boost::num_vertices(_graph));

	boost::dijkstra_shortest_paths(_graph, sv, boost::predecessor_map(&pred[0]).distance_map(&dist[0]));

	// 7. evaluate result;
	auto stamp7 = clock();
	auto vie = eTrn + eRotIdx;
	auto vis = sTrn + sRotIdx;

	// 7-1. failure
	if (pred[vie] == vie)
		return false;

	// 7-2. path of vertIdx;
	vector<int> path;
	int cur = vie;
	while (cur != vis)
	{
		path.push_back(cur);
		cur = pred[cur];
	}
	path.push_back(vis);

	reverse(path.begin(), path.end());

	// 7-3. path xyt
	auto& out = _out_path;
	for (auto i : path)
	{
		if (i < sTrn)
		{
			// v from original graph
			out.push_back(getVert()[i]);
		}
		else if (i < eTrn)
		{
			// v from start-connection
			auto localIdx = i - sTrn;
			auto z = rotatableS[localIdx];
			if (z < 0)
				z += ns;
			if (z >= ns)
				z -= ns;
			double t = double(z) / ns * 360.0;
			out.emplace_back(start.x(), start.y(), t);

			// dbg_out
			{
				if (t < 0.0 || t >= 360.0)
					cout << "t is strange: " << t << endl;
			}
		}
		else
		{
			// v from end-connection
			auto localIdx = i - eTrn;
			auto z = rotatableE[localIdx];
			if (z < 0)
				z += ns;
			if (z >= ns)
				z -= ns;
			double t = double(z) / ns * 360.0;
			out.emplace_back(end.x(), end.y(), t);

			// dbg_out
			{
				if (t < 0.0 || t >= 360.0)
					cout << "t is strange: " << t << endl;
			}
		}
	}

	// 8. clr list build;
	{
		auto& outclr = _out_clr;

		for (auto i : path)
		{
			if (i < sTrn)
			{
				outclr.push_back(getClr()[i]);
			}
			else if (i < eTrn)
			{
				// v from start-connection
				auto localIdx = i - sTrn;
				outclr.push_back(rotatableSclr[localIdx]);
			}
			else
			{
				auto localIdx = i - eTrn;
				outclr.push_back(rotatableEclr[localIdx]);
			}
		}
	}

	// 9. edge type build
	{
		auto& outet = _out_et;
		for (int i = 0; i < path.size() - 1; i++)
		{
			// get ed;
			//int i1 = i + 1;
			auto v0 = path[i + 0];
			auto v1 = path[i + 1];
			//auto vd0 = boost::vertex(v0, _graph);
			//auto vd1 = boost::vertex(v1, _graph);
			//auto ed = boost::edge(vd0, vd1, _graph).first;

			// find edgeType
			if (v0 < sTrn)
			{
				if (v1 < sTrn)
				{
					// [0, 6)
					auto query = _invEdge.find(make_pair(v0, v1));
					if (query == _invEdge.end())
					{
						auto query2 = _invEdge.find(make_pair(v1, v0));
						if (query2 == _invEdge.end())
						{
							cout << " Err at searchGraph2 : 9 : _invEdge not found " << endl;
							outet.push_back(gs::etUnknow);
							continue;
						}

						auto idx = query2->second;
						outet.push_back(_etype[idx]);
						continue;
					}

					auto idx = query->second;
					outet.push_back(_etype[idx]);
				}
				else if (v1 < eTrn)
				{
					// 7
					outet.push_back(gs::etStrHor);
				}
				else
				{
					// 9
					outet.push_back(gs::etEndHor);
				}
			}
			else if (v0 < eTrn)
			{
				if (v1 < sTrn)
				{
					// 7
					outet.push_back(gs::etStrHor);
				}
				else if (v1 < eTrn)
				{
					// 6
					outet.push_back(gs::etStrVer);
				}
				else
				{
					// 10
					outet.push_back(gs::etStrEnd);
				}
			}
			else
			{
				if (v1 < sTrn)
				{
					// 9
					outet.push_back(gs::etEndHor);
				}
				else if (v1 < eTrn)
				{
					// 10
					outet.push_back(gs::etStrEnd);
				}
				else
				{
					// 8
					outet.push_back(gs::etEndVer);
				}
			}
		}
	}

	// 99. remove vert/edge
	auto stamp8 = clock();
	{
		//No?
		if (!alterDjkCalc) 
		{
			//// remove edge
			//for (auto& vd : svd)
			//	boost::clear_vertex(vd, _graph);
			//for (auto& vd : evd)
			//	boost::clear_vertex(vd, _graph);
			//
			//// remove vert
			//int s = sTrn;
			//int e = boost::num_vertices(_graph);
			//for (int i = e - 1; i >= s; i--)
			//{
			//	auto vd = boost::vertex(i, _graph);
			//	boost::remove_vertex(vd, _graph);
			//}

			//_graph.swap(cpyGraph); // not used anymore: this is time consuming
			_graph = (cpyGraph);
		}
	}
	auto stamp9 = clock();
	cout << "time: init:   " << double(stamp1 - stamp0) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	cout << "time: rotAng: " << double(stamp2 - stamp1) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	cout << "time: connec: " << double(stamp3 - stamp2) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	cout << "time: addVer: " << double(stamp4 - stamp3) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	cout << "time: addEho: " << double(stamp5 - stamp4) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	cout << "time: addEvr: " << double(stamp6 - stamp5) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	cout << "time: search: " << double(stamp7 - stamp6) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	cout << "time: wrtPth: " << double(stamp8 - stamp7) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	cout << "time: rmvGrp: " << double(stamp9 - stamp8) / CLOCKS_PER_SEC * 1000 << "ms" << endl;

	return true;
}

/*
Def:
	multiplier constant representing non-holo
Assume:
	p != q
*/
gs::Weight djkCalc::nhMultiplier(const Point& p, const Point& q, Point& forward, double sliceDegree)
{
	auto dot = (p - q).normalize() * forward.rotate(sliceDegree / 360.0 * PI2);
	//auto theta = acos(dot);
	//auto degree = theta * 180.0 / PI;
	//auto flr = floor(degree / 10.0);
	//return flr + 1.0;

	//// if (angle < 25 deg)
	//if (abs(dot) > 0.90)
	//	return 1.0;
	//else
	//	return 1e100;

	return floor(((1 - abs(dot)) * 20.0)) + 1.0;
	return floor(((1 - abs(dot)) * 10.0)) + 1.0;
}

gs::Weight djkCalc::findWeight(const Point& p, const Point& q, Point& forward, double sliceDegree)
{
	auto len = (p - q).length1();
	return nhMultiplier(p, q, forward, sliceDegree) * len;
	return len;
}