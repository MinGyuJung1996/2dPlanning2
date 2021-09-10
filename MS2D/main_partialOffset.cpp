#include <iostream>
#include <fstream>
#include "voronoi.hpp"
#include "AStarOnVorDiag.h"
#include "collision detection.hpp"
#include "support.hpp"
#include "dijkstra.hpp"
#include "distance field.hpp"
using namespace std;
extern double vbRad;

#define COUT(x) #x ": " << x 
#define ONLY_USE_VORONOI_IN_GRAPH true

constexpr double _h_offVal = 0.01; // offset value of offCrv
constexpr double _h_offMin = 0.001; // min offset needed for offCrv to be calculated 
constexpr double _h_offTan = 0.05; // vor with clr of this val will shoot tangents

/*******************************************************************************************************************************
**																															  **
**														EMPTY BOX															  **
**																															  **
*******************************************************************************************************************************/

#pragma region NS
namespace ms
{
	int main(int argc, char *argv[]);
	int main2(int argc, char* argv[]);
	void renderMinkVoronoi(
		int argc,
		char* argv[],
		std::vector<decltype(ms::Model_Result)>& MRs,
		std::vector<decltype(ms::ModelInfo_Boundary)>& MIBs,
		std::vector<std::vector<planning::output_to_file::v_edge>>& v_edges,
		decltype(planning::voronoiBoundary)& voronoiBoundary
		);
	void renderPath(
		int argc,
		char* argv[],
		std::vector<double>& path
		);
	void renderRefinementCollisionTest(
		int argc,
		char* argv[],
		std::vector<decltype(ms::Model_Result)>& MRs,
		std::vector<decltype(ms::ModelInfo_Boundary)>& MIBs,
		std::vector<std::vector<planning::output_to_file::v_edge>>& v_edges,
		decltype(planning::voronoiBoundary)& voronoiBoundary,
		std::vector<planning::VR_IN> VRINs
		);
	namespace renderPathGlobal
	{
		extern vector<Vertex> vecVertices;
		extern map<Vertex, int, VertexLessFn> mapLookup;
		extern Graph theGr;
	}
}
namespace rendering3D
{
	void renderCSObject(
		int argc,
		char* argv[],
		std::vector<decltype(ms::Model_Result)>& MRs,
		std::vector<decltype(ms::ModelInfo_Boundary)>& MIBs
	);
	namespace main3
	{
		int main3(int argc, char* argv[]);
	}
}
namespace graphSearch
{

	int main(int argc, char** argv);
	void searchTest();
	int main2(int argc, char* argv[]);
}
namespace sceneEditor
{
	int main(int argc, char* argv[]);
}
namespace modelEditor
{
	int main(int argc, char* argv[]);
}
#pragma endregion

namespace graphSearch_partialOffset
{
	using v_edge = planning::output_to_file::v_edge;
	std::vector<cd::lineSegmentCollisionTester> testers;
	std::vector<cd::pointCollisionTester> testers2;

	using lseg = v_edge;
	//std::vector<std::vector<std::vector<lseg>>> offCrv; // offCrv[offIdx][sliceNo][lsegIdx];

	std::vector<decltype(ms::Model_Result)>					MRs(ms::numofframe); // data collected for checking
	std::vector<decltype(ms::ModelInfo_Boundary)>			MIBs(ms::numofframe); // data collected for checking
	std::vector<decltype(ms::InteriorDisks_Convolution)>	IDC(ms::numofframe); // save conv/disk

	std::vector<planning::VR_IN>							VRINs(ms::numofframe);
	std::vector<std::vector<v_edge>>						VROUTs(ms::numofframe);

	std::vector<std::vector<v_edge>>						OffCrv(ms::numofframe); // store offset curves
	std::vector<std::vector<v_edge>>						VorTan(ms::numofframe); // store voronoi-offset tangent := offCon
	std::vector<std::vector<v_edge>>						ComTan(ms::numofframe); // store tangent between models
	std::vector<std::vector<v_edge>>						Connector(ms::numofframe); // connection between geometries

	std::vector<std::vector<SupportComplex>>				SCs(ms::numofframe);
	
	//std::vector<decltype(ms::Model_Result)>				offsetMRs(ms::numofframe);	// data collected for checking
	//std::vector<decltype(ms::ModelInfo_Boundary)>			offsetMIBs(ms::numofframe); // data collected for checking
	//std::vector<decltype(ms::InteriorDisks_Convolution)>	offsetIDC(ms::numofframe);  // save conv/disk

	djkCalc dijk;

	// nav rob 
	// nav obs

	int robotType = 1;
	/*
	0 : rect
	1 : 8
	2 : sq
	3 : 8-2
	*/
	int obsType = 15; // nav obs
	/* options: (obs)
	3: curved S corridor
	1,4 1,5
	1,6~11 (made with modelEditor)
	*12: pg2021
	13: maze
	*15: magnified 3
	*16: maze refined
	*17: ICRA
	*/

	/*
	Def:
		Calc Geometry for Motion Planning
		Then, call renderPath
	Summary:
		1. build Mink
		2. build Voronoi
			2-1. build clr of each pt.
		3. build CollisionTester
		4. Find tangents
			4-1. find emitting points
		5. Build offset curves.
		6. Make searchGraph?
		7. call renderPath?
	TODO:
		1. sup-class for point           |
		2. comTanCalc for point/model	 |-> handled with arc with very small rad, maybe 0.0;
		3. offset for some case / criteria clr =/= offset dist ???
		4. connect graph
		5. do gs
		6. render path
	*/
	int main(int argc, char* argv[])
	{
		/*
		
		*/

		// 0. preprocessing stuff
		{
			planning::_h_fmdsp_g1 = 1e-6;
			initializeRobotObstacles(robotType, obsType);	// read data from files.
			ms::ModelInfo_CurrentModel = std::make_pair(1, 7);
			ms::postProcess(ms::ModelInfo_CurrentModel.first, ms::ModelInfo_CurrentModel.second); // process arcs to satisfy conditions.
			planning::output_to_file::flag = true; // flag to enable : gathering data inside global var, during v-diagram construction.
			planning::output_to_file::v_edges.resize(0); //empty
			planning::output_to_file::v_edges.resize(ms::numofframe); //pre-allocate
			planning::output_to_file::bifur_points.resize(0); // dummy
			planning::output_to_file::bifur_points.resize(ms::numofframe); //dummy
			planning::drawVoronoiSingleBranch = false; //disable drawing for now.
		}

		/*******************************************************************************************************************************
		**																															  **
		**												1 ~ 3 . mink / vor / collision												  **
		**																															  **
		*******************************************************************************************************************************/

		/***
		
		*/


		// 1 + 2. build mink
		//auto hrc0 = std::chrono::high_resolution_clock::now();
		auto mvStartTime = clock();

		for (size_t i = 0, length = ms::numofframe /* = 360*/; i < length; i++)
		{
			// 1-1. minks
			ms::t2 = i;
			ms::minkowskisum(i, 7);

			// save data for checking
			MRs [i] = ms::Model_Result;
			MIBs[i] = ms::ModelInfo_Boundary;
			IDC [i] = ms::InteriorDisks_Convolution;

			// 1-2. Vor
			planning::VR_IN& vrin = VRINs[i];
			planning::_Convert_MsOut_To_VrIn(ms::Model_Result, ms::ModelInfo_Boundary, vrin);
			voronoiCalculator vc;
			vc.initialize();
			vc.setInput(vrin.arcs, vrin.left, vrin.color);

			// call diff with vType
			constexpr int vType = 1;
			if (vType == 1)
			{
				vector<deque<VoronoiEdge>> v_res;
				vc.setOutput(v_res);
				vc.calculate();

				//auto& result = planning::output_to_file::v_edges[i];
				auto& result = VROUTs[i];
				for (auto& a : v_res)
					for (auto& b : a)
						result.push_back(b);
			}
			else if (vType == 2)
			{
				voronoiCalculatorResultG vcg;
				vc.setOutputG(vcg);
				vc.calculateG();
				
				//auto& result = planning::output_to_file::v_edges[i];
				auto& result = VROUTs[i];
				for (auto& a : vcg.E())
					for (auto& b : a)
						result.push_back(b);
			}
		}

		planning::output_to_file::flag = false;
		auto mvEndTime = clock();
		cout << "mink + voronoi time : "
			<< double(mvEndTime - mvStartTime) / 1000 << "s" << endl;

		// 2-1. set clearance.
		for (int i = 0; i < ms::numofframe; i++)
			planning::setClearance(VROUTs[i], VRINs[i].arcs);

		// 3. build collisionTesters 
		auto cdInitStartTime = clock();

		// 3-1. line tester
		testers.resize(ms::numofframe);
		for (size_t i = 0; i < ms::numofframe; i++)
		{
			testers[i].initialize(VRINs[i]);
		}

		// 3-2. point tester
		testers2.resize(ms::numofframe);
		for (size_t i = 0; i < ms::numofframe; i++)
		{
			vector<bool> convexityList(VRINs[i].arcs.size());
			for (int j = 0; j < convexityList.size(); j++)
			{
				convexityList[j] = !VRINs[i].arcs[j].ccw;
			}

			testers2[i].initialize(VRINs[i], convexityList, IDC[i]);
		}
		auto cdInitEndTime = clock();

		cout << "collision tester init time : "
			//<< ms::numofframe << " : " 
			<< double(cdInitEndTime - cdInitStartTime) / 1000 << "s" << endl;


		//// renderminkvoronoi : mink/ vor only
		//if (0)
		//{
		//	std::vector<std::vector<v_edge>>&
		//	//	v_edges = planning::output_to_file::v_edges;
		//		v_edges = VROUTs;
		//
		//	_h_fmdsp_g1 = 1e-3;
		//	ms::renderMinkVoronoi(argc, argv, MRs, MIBs, v_edges, planning::voronoiBoundary);
		//	//ms::renderMinkVoronoi(argc, argv, offsetMRs, offsetMIBs, v_edges, planning::voronoiBoundary);
		//}


		/*******************************************************************************************************************************
		**																															  **
		**													4. build offset curve													  **
		**																															  **
		*******************************************************************************************************************************/

		/*
		Assume clr built at this point
		*/

		// 4-1. build offset
		for (int i = 0; i < ms::numofframe; i++)
		{
			OffCrv[i] = planning::findOffsetCurve(VROUTs[i], VRINs[i].arcs, _h_offVal, _h_offMin);
		}

		/*******************************************************************************************************************************
		**																															  **
		**													5. build tangent/connection												  **
		**																															  **
		*******************************************************************************************************************************/

		// 5-0. build support complex

		// for each slice
		for (int i = 0; i < ms::numofframe; i++)
		{
			auto& mr = MRs[i];
			auto& mib = MIBs[i];
			CommonTangentCalculator ctc;

			// for each loop
			for (int j = 0; j < mr.size(); j++)
			{
				// if (inner loop) continue;
				if (mib[j] == false)
					continue;

				SCs[i].emplace_back();
				SupportComplex& sc = SCs[i].back();
				{
					vector<CircularArc> temp;
					for (auto& as : mr[j])
						for (auto& arc : as.Arcs)
							temp.push_back(arc);

					sc = ctc.findEnvelope(temp);
				}
			}
		}

		// 5-1. find tangent candidates
		auto offsetStartTime = clock();
		vector<vector<Point>> tangentCandidates(ms::numofframe);
		{
			auto& vrOuts = planning::output_to_file::v_edges;

			// for (each slice)
			for (int i = 0; i < ms::numofframe; i++)
			{
				auto& tc = tangentCandidates[i];
				auto& vrOut = vrOuts[i];

				// for (each vorEdge)
				for (int j = 0; j < vrOut.size(); j++)
				{
					// check if clrX < _h_offset1 < clrY holds
					auto& ve = vrOut[j];

					bool flag0 = ve.clr0() < _h_offTan;
					bool flag1 = ve.clr1() > _h_offTan;

					// if (not between clrs)
					if (flag0 ^ flag1)
					{
						// do nothing
					}
					else
					{
						// insert one with larger
						if (ve.clr0() > ve.clr1())
						{
							tc.push_back(ve.v0);
							// TODO : conn? instead of nearest
						}
						else
						{
							tc.push_back(ve.v1);
							// TODO : conn?
						}
					}
				}
			}
		}

		// 5-2. for each candidate, find tangents
		vector<vector<lseg>> tangentLines;
		{
			for (int i = 0; i < ms::numofframe; i++)
			{
				auto& tc = tangentCandidates[i];
				for (int j = 0; j < tc.size(); j++)
				{
					auto& p = tc[j];

					//Support sup(p); // TODO
					//SupportComplex supc(sup);
				}
			}
		}
		/*******************************************************************************************************************************
		**																															  **
		**													  6. build Dijkstra Graph												  **
		**																															  **
		*******************************************************************************************************************************/
		auto gt0 = clock();
		//dijk.buildGraph(ms::numofframe, vorBackup, offCso, offCon, comTan, comTanCon, robotForward);
		auto gt1 = clock();
		cout << COUT(robotForward) << endl;
		cout << "Dijk init time : "
			<< double(gt1 - gt0) / 1000 << "s" << endl;
		{
			auto gt0 = clock();
			auto dijk2 = dijk.getGraph();
			auto gt1 = clock();
			cout << "Dijk copy time : "
				<< double(gt1 - gt0) / 1000 << "s" << endl;
		}


		///* test if result is same : print all v-edges to file and compare it*/
		//std::ofstream fout("ve_out.txt");
		//for (size_t i = 0, length = ms::numofframe; i < length; i++)
		//{
		//	fout << "voronoi " << i << std::endl;
		//	fout << planning::output_to_file::v_edges[i].size() << std::endl;
		//	for (auto e : planning::output_to_file::v_edges[i])
		//		fout << e.v0 << "\t" << e.v1 << std::endl;

		//}
		//std::cout << "END OF TEST" << std::endl;

		// 2. Alias for long names

		//using v_edge = planning::output_to_file::v_edge;
		std::vector<std::vector<v_edge>>&
			v_edges = planning::output_to_file::v_edges;

		if (0)
		{
			_h_fmdsp_g1 = 1e-3;
			ms::renderMinkVoronoi(argc, argv, MRs, MIBs, v_edges, planning::voronoiBoundary);
			//ms::renderMinkVoronoi(argc, argv, offsetMRs, offsetMIBs, v_edges, planning::voronoiBoundary);
		}
		auto gsInitStartTime = clock();

		std::vector<double> renderedPath;
		if (false)
		{
			//Graph theGr = create_VorGraph(v_edges);
			//std::vector<v_edge> path = invoke_AStar(theGr);
			std::vector<std::vector<v_edge>> v_edges_sparced;
			std::vector<std::vector<v_edge>>::iterator iterEdgesLayer = v_edges.begin();
			for (size_t i = 0; i < v_edges.size(); ++i, ++iterEdgesLayer)
				//if (0 == (i % 10)) // filter out intermediate layers
				v_edges_sparced.push_back(*iterEdgesLayer);

			vector<Vertex> vecVertices;
			map<Vertex, int, VertexLessFn> mapLookup;
			::Graph theGr = create_VorGraph(v_edges_sparced, vecVertices, mapLookup);
			//Vertex ptnSrc(-0.70296182284311681, -0.30610712038352472, 0.0);
			//Vertex ptnDst(0.76932775901415118, 0.36524457774288216, 320.0);

			// Set starting pt, ending pt
			Vertex ptnSrc(-0.9, -0.1, 0.0);
			Vertex ptnDst(+0.5, +0.4, 160.0);
			if (obsType == 2)
			{

				// 4
				ptnSrc = Vertex(-1.8, -1.5, 0.0);
				ptnDst = Vertex(1.0, -1.1, 180.0);

				// 3
				//ptnSrc = Vertex(-1.8, -1.5, 0.0);
				//ptnDst = Vertex(2.2, 1.3, 90.0);

				// 2
				//ptnSrc = Vertex(-1.8, -1.5,  0.0);
				//ptnDst = Vertex( 2.5, 2, 90.0); 
				//
				//// 1
				//ptnSrc = Vertex(-1.8, -1.5, 0.0);
				//ptnDst = Vertex(1.0, -1.1, 90.0);
			}
			if (/*robotType == 0 && */obsType == 1)
			{
				ptnSrc = Vertex(0.858, -1.07, 0.0);
				ptnDst = Vertex(-0.502, 0.509, 0.0);
			}
			// Code to find the closest point in G(V,E)
			{
				Vertex src(0, 0, 0), dst(0, 0, 0);
				double d0 = 1e100, d1 = 1e100;
				for (auto& a : vecVertices)
				{
					auto& v = a;
					if (v.z == ptnSrc.z)
					{
						if (v.dist(ptnSrc) < d0)
						{
							d0 = v.dist(ptnSrc);
							src = v;
						}
					}
					if (v.z == ptnDst.z)
					{
						if (v.dist(ptnDst) < d1)
						{
							d1 = v.dist(ptnDst);
							dst = v;
						}
					}
				}
				ptnSrc = src;
				ptnDst = dst;
			}

			auto gsInitEndTime = clock();
			cout << "slice connection time : "
				<< double(gsInitEndTime - gsInitStartTime) / 1000 << "s" << endl;

			auto iaStart = clock();
			std::vector<Vertex> path = invoke_AStar(theGr, vecVertices, mapLookup, ptnSrc, ptnDst);
			auto iaEnd = clock();
			cout << "A-star time : "
				<< double(iaEnd - iaStart) / 1000 << "s" << endl;

			//std::cout << path.size() << std::endl;

			// 3. call functions from namespace graphSearch (AStarOnVorDiag.cpp)
			 // triplet of (path[3n+0], path[3n+1], path[3n+2]) represents a vertice in path. 
			// vertices = ...


			// 4~. do sth with the path....

			// 4-1. Just to check whether mink/vor was constructed properly.
			// uncomment below to begin renderLoop for mink/voronoi calculated above.
			//ms::renderMinkVoronoi(argc, argv, MRs, MIBs, v_edges, planning::voronoiBoundary);
			//ms::renderRefinementCollisionTest(argc, argv, MRs, MIBs, v_edges, planning::voronoiBoundary, VRINs);
			//rendering3D::renderCSObject(argc, argv, MRs, MIBs);

			// 4-2. render robot's path
			// Path found on step 3 should be used instead of dummy_path
			{
				// just a fake path to check program pipeline.
				//for (int i = 0; i < 100; i++)
				for (auto v : path)
				{
					//double x = -0.7 + 1.4 / 100.0 * i; // x-coord of Robot center
					//double y = -0.5 + (i / 100.0) * (i / 100.0); // y-coord of Robot center
					//double z = log10(i + 1) * 180;	// should be rotation in degrees

					renderedPath.push_back(v.x);
					renderedPath.push_back(v.y);
					renderedPath.push_back(v.z);
				}
			}
			bool savePathAsFile = true;
			if (savePathAsFile)
			{
				ofstream fout("path.txt");
				fout << renderedPath.size() / 3 << endl;
				for (int i = 0; i < renderedPath.size(); i += 3)
				{
					fout
						<< renderedPath[i + 0] << " "
						<< renderedPath[i + 1] << " "
						<< renderedPath[i + 2] << endl;
				}
			}

			// ready for renderPath
			ms::renderPathGlobal::vecVertices = vecVertices;
			ms::renderPathGlobal::mapLookup = mapLookup;
			ms::renderPathGlobal::theGr = theGr;

			//{
			//	for (int i = 0; i< v_edges.size(); i++)
			//	{
			//		auto& ves = v_edges[i];
			//		for (auto& ve : ves)
			//		{
			//			auto& p = ve.v0;
			//			if (ve.idx[0] < 0)
			//				continue;
			//			auto& c = VRINs[i].arcs[ve.idx[0]];
			//			auto d = (p - c.c.c).length2();
			//			d = sqrt(d);
			//
			//			ve.clr0() = fabs(d - c.cr());
			//		}
			//	}
			//}
		}

		ms::renderPath(argc, argv, renderedPath);

		return 0;
	}

}
