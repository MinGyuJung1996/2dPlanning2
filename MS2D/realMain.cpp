#include <iostream>
#include <fstream>
#include "voronoi.hpp"
#include "AStarOnVorDiag.h"
#include "collision detection.hpp"

using namespace std;
extern double vbRad;

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

int main(int argc, char *argv[]) {

	//Test stuff:
	if (false)
	{
		using namespace ms;
		vector<CircularArc> cycle;

		CircularArc c0, c1, c2;
		/*
		c0.ccw = false;
		c1.ccw = false;
		c2.ccw = false;
		c0.c.c = Point( 2 + 0, 0 + 0);
		c0.c.r = 1;
		c1.c.c = Point(-2 + 0, 0 + 0);
		c1.c.r = 1;
		c2.c.c = Point( 0 + 0, 2 + 0);
		c2.c.r = 1;*/
		
		/*c0.c.c = Point(0, 4);
		c0.c.r = 1;
		c1.c.c = Point(7, 7);
		c1.c.r = 3 * sqrt(2) - 3;
		c2.c.c = Point(4, -1);
		c2.c.r = 2;
		cycle.push_back(c0);
		cycle.push_back(c1);
		cycle.push_back(c2);*/
		/*
		c0.ccw = false;
		c1.ccw = false;
		c2.ccw = true;
		c0.c.c = Point(0.564474, -0.103322);
		c0.c.r = 0.246477;
		c1.c.c = Point(0.549305, -0.529539);
		c1.c.r = 0.0423848;
		c2.c.c = Point(0, 0);
		c2.c.r = 2;*/

		c0.ccw = false;
		c1.ccw = false;
		c2.ccw = true;
		c0.c.c = Point(0.5, -0.1);
		c0.c.r = 0.2;
		c1.c.c = Point(0.5, -0.5);
		c1.c.r = 0.1;
		c2.c.c = Point(0, 0);
		c2.c.r = 2;


		cycle.push_back(c0);
		cycle.push_back(c1);
		cycle.push_back(c2);


		Point bifur(0, 0);
		double R;
		{
			// + if 
			// # := +- 
			// $ := -+
			// 3 unkown x, y, R (point, R)
			// a0^2 - 2a0x + x^2 + b0^2 - 2b0y + y^2 = r0^2 + 2 * s0 * r0R + R^2 ... (1)
			// a1^2 - 2a1x + x^2 + b1^2 - 2b1y + y^2 = r1^2 + 2 * s1 * r1R + R^2 ... (2)
			// a2^2 - 2a2x + x^2 + b2^2 - 2b2y + y^2 = r2^2 + 2 * s2 * r2R + R^2 ... (3)
			//
			// (a0^2 - a1^2) -2(a0 - a1)x + (b0^2 - b1^2) -2(b0 - b1)y = (r0^2 - r1^2) + 2(s0r0-s1r1)R ...(4)
			// (a1^2 - a2^2) -2(a1 - a2)x + (b1^2 - b2^2) -2(b1 - b2)y = (r1^2 - r2^2) + 2(s1r1-s2r2)R ...(5)
			//
			// [ 2(a0 - a1) 2(b0 - b1) ] [x]  =  [ -2(s0r0-s1r1)R + {(a0^2 - a1^2) + (b0^2 - b1^2) - (r0^2 - r1^2)}
			// [ 2(a1 - a2) 2(b1 - b2) ] [y]     [ -2(s1r1-s2r2)R +	{(a1^2 - a2^2) + (b1^2 - b2^2) - (r1^2 - r2^2)} ... (6)
			// Az = cR + d; [x, y] = [f(r), g(r)]
			// if 3 circle centers are in a straight line... (do something else)
			//
			// plugging x, y in (6) to (1)~(3) will give two solutions of R... it will give one pos & one negative i guess...? -> since solution for (s0,s1,s2) case r0,r1 has duality with solution for (-s0,-s1,-s2)'s -r0,-r1 


			// 1-1. first check if three circle centers are colinear.
			bool isColinear;
			{
				double det = (cycle[0].c.c - cycle[1].c.c).normalize() ^ (cycle[1].c.c - cycle[2].c.c).normalize();
				isColinear = fabs(det) < 1e-10;
			}

			if (true) // TODO? when colinear
			{
				// 1-2. make x, y = cr + d

				// 1-2-1. find sign ($ above)
				double sign0 = cycle[0].ccw ? -1 : 1;
				double sign1 = cycle[1].ccw ? -1 : 1;
				double sign2 = cycle[2].ccw ? -1 : 1;

				// 1-2-2. apply inv mat
				double
					&a0 = cycle[0].c.c.P[0], &b0 = cycle[0].c.c.P[1], &r0 = cycle[0].c.r,
					&a1 = cycle[1].c.c.P[0], &b1 = cycle[1].c.c.P[1], &r1 = cycle[1].c.r,
					&a2 = cycle[2].c.c.P[0], &b2 = cycle[2].c.c.P[1], &r2 = cycle[2].c.r;

				auto det = (a0 - a1) * (b1 - b2) - (b0 - b1) * (a1 - a2);
				auto c0 = -(b1 - b2) * (sign0 * r0 - sign1 * r1) + (b0 - b1) * (sign1 * r1 - sign2 * r2); c0 /= det;
				auto c1 = +(a1 - a2) * (sign0 * r0 - sign1 * r1) - (a0 - a1) * (sign1 * r1 - sign2 * r2); c1 /= det;
				auto e0 = +(a0 - a1) * (a0 + a1) + (b0 - b1) * (b0 + b1) - (r0 - r1) * (r0 + r1);
				auto e1 = +(a1 - a2) * (a1 + a2) + (b1 - b2) * (b1 + b2) - (r1 - r2) * (r1 + r2);
				auto d0 = +(b1 - b2) * e0 - (b0 - b1) * e1; d0 /= 2 * det;
				auto d1 = -(a1 - a2) * e0 + (a0 - a1) * e1; d1 /= 2 * det;

				// 1-3. solve (1) (= aR^2+bR+c) with x = c0*R + d0, y = c1*R + d1;
				// TODO : can save some inst by using b/2 = b'
				d0 -= a0;
				d1 -= b0;

				double
					a = c0*c0 + c1*c1 - 1 * 1,
					b = 2 * c0*d0 + 2 * c1*d1 - 2 * sign0*r0,
					c = d0*d0 + d1*d1 - r0*r0;

				if (fabs(a) < 1e-100) R = -c / b;
				else R = min(((-b + sqrt(b*b - 4 * a*c)) / a) * 0.5, ((-b - sqrt(b*b - 4 * a*c)) / a) * 0.5);

				// 1-4. compute point of bifurcation
				d0 += a0;
				d1 += b0;

				bifur = Point(c0*R + d0, c1*R + d1);

				cout << a << " " << b << " " << c << endl;

				
			}
			cout << "dist0 " << sqrt((c0.c.c - bifur).length()) << endl;
			cout << "dist1 " << sqrt((c1.c.c - bifur).length()) << endl;
			cout << "dist2 " << sqrt((c2.c.c - bifur).length()) << endl;

		}

		cout << "bifur " << bifur << endl;
		cout << "R " << R << endl;
	}

	//cout << "fake func" << endl;

	//sceneEditor::main(argc, argv);	// interactive editing
	//rendering3D::main3::main3(argc, argv);	// draw 3d c-space
	graphSearch::main2(argc, argv);		// draw moving robot
	//graphSearch::main(0, NULL);
	//ms::main2(argc, argv); //RSV
	//ms::main(argc, argv); // mink test

	system("pause");
}

namespace graphSearch
{
	using v_edge = planning::output_to_file::v_edge;
	std::vector<cd::lineSegmentCollisionTester> testers;
	std::vector<cd::pointCollisionTester> testers2;


	/*
	Def:
		Computes Voronoi Edges for each slice, and store them in a vector<vector<v_edge>> v_edges.
		Each Voronoi Edge can be referenced by v_edges[Slice_Number][Edge_Number_In_That_Slice];
	*/
	int main2(int argc, char* argv[])
	{
		int robotType = 0;
		int obsType   = 0;

		//minor stuffs with model
		{
			//if (obsType == 0)
			//	vbRad = 0.7;
			//if (obsType == 2)
			//	ms::zoom = 3.0;
		}

		// 1. build v_edges
		planning::_h_fmdsp_g1 = 1e-8;
		initializeRobotObstacles(robotType, obsType);	// read data from files.
		ms::ModelInfo_CurrentModel = std::make_pair(1, 7);
		ms::postProcess(ms::ModelInfo_CurrentModel.first, ms::ModelInfo_CurrentModel.second); // process arcs to satisfy conditions.
		planning::output_to_file::flag = true; // flag to enable : gathering data inside global var, during v-diagram construction.
		planning::output_to_file::v_edges.resize(0); //empty
		planning::output_to_file::v_edges.resize(ms::numofframe); //pre-allocate
		planning::output_to_file::bifur_points.resize(0); // dummy
		planning::output_to_file::bifur_points.resize(ms::numofframe); //dummy
		planning::drawVoronoiSingleBranch = false; //disable drawing for now.
		std::vector<decltype(ms::Model_Result)>					MRs(ms::numofframe); // data collected for checking
		std::vector<decltype(ms::ModelInfo_Boundary)>			MIBs(ms::numofframe); // data collected for checking
		std::vector<decltype(ms::InteriorDisks_Convolution)>	IDC(ms::numofframe); // save conv/disk
		std::vector<planning::VR_IN>							VRINs(ms::numofframe);
		std::vector<decltype(ms::Model_Result)>					offsetMRs(ms::numofframe);	// data collected for checking
		std::vector<decltype(ms::ModelInfo_Boundary)>			offsetMIBs(ms::numofframe); // data collected for checking
		std::vector<decltype(ms::InteriorDisks_Convolution)>	offsetIDC(ms::numofframe); // save conv/disk

		auto mvStartTime = clock();
		for (size_t i = 0, length = ms::numofframe /* = 360*/; i < length; i++)
		{
			// 1-1. minks
			ms::t2 = i;
			ms::minkowskisum(i, 7);
			// save data for checking
			MRs[i] = ms::Model_Result;
			MIBs[i] = ms::ModelInfo_Boundary;
			IDC[i] = ms::InteriorDisks_Convolution;

			// 1-2. Vor
			planning::VR_IN& vrin = VRINs[i];
			planning::_Convert_MsOut_To_VrIn(ms::Model_Result, ms::ModelInfo_Boundary, vrin);
			//planning::_Medial_Axis_Transformation(vrin);
			voronoiCalculator vc;
			vector<deque<VoronoiEdge>> v_res;
			vc.initialize();
			vc.setInput(vrin.arcs, vrin.left, vrin.color);
			vc.setOutput(v_res);
			vc.calculate();

			auto& result = planning::output_to_file::v_edges[i];
			for (auto& a : v_res)
				for (auto& b : a)
					result.push_back(b);

			// 1-3. Offset Curve;

			//// dbg_out
			//std::cout << "voronoi " << i << " "
			//	<< planning::output_to_file::v_edges[i].size() << std::endl;
		}
		planning::output_to_file::flag = false;
		auto mvEndTime = clock();
		cout << "mink + voronoi time : "
			<< double(mvEndTime - mvStartTime) / 1000 << "s" << endl;

		
		//// 1-2. build collisionTesters 
		auto cdInitStartTime = clock();
		testers.resize(ms::numofframe);
		for (size_t i = 0; i < ms::numofframe; i++)
		{
			testers[i].initialize(VRINs[i]);
		}

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
			<< double(cdInitEndTime - cdInitStartTime) / 1000 <<"s"<< endl;

		// 1-3. build offset curves.
		auto offsetStartTime = clock();
		
		//// mink sum approach
		//{
		//	// use 2 for offsetCircle idx;
		//	int cirIdx = 2;
		//	double offsetRadius = 0.1;
		//	{
		//		// disk
		//		vector<Circle> interiorDisk;
		//		Circle c = Circle(Point(0, 0), offsetRadius * 0.9);
		//		interiorDisk.push_back(c);
		//
		//		// vca
		//		vector<CircularArc> vca;
		//		double eps = 1e-6;
		//		CircularArc c0 = cd::constructArc(Point(0, 0), offsetRadius, 0  + eps, PI  -eps);
		//		CircularArc c1 = cd::constructArc(Point(0, 0), offsetRadius, PI + eps, PI2 -eps);
		//		vca.push_back(c0);
		//		vca.push_back(c1);
		//		
		//		// vas
		//		vector<ArcSpline> vas;
		//		planning::_Convert_VectorCircularArc_To_MsInput(vca, vas, 0.0);
		//
		//		// Assign
		//		ms::Models_Approx			[cirIdx] = vas;
		//		ms::InteriorDisks_Imported	[cirIdx] = interiorDisk;
		//		ms::Model_vca				[cirIdx] = vca;
		//		ms::Model_from_arc			[cirIdx] = true;
		//	}
		//
		//	// set model;
		//	for (int i = 0; i< ms::numofframe; i++)
		//	{
		//		ms::InteriorDisks_Rotated[i] = IDC[i];
		//
		//		ms::Models_Rotated_Approx[i].clear();
		//		vector<CircularArc> vca = VRINs[i].arcs;
		//		auto size = vca.size() - voronoiBoundary.size();
		//		vca.resize(size);
		//		_Convert_VectorCircularArc_To_MsInput(vca, ms::Models_Rotated_Approx[i], 0.0);
		//
		//		//for (auto& loop : MRs[i])
		//		//for (auto& as: loop)
		//		//{
		//		//	ms::Models_Rotated_Approx[i].push_back(as);
		//		//}
		//	}
		//
		//	// do mink and save
		//	for (int i = 0; i < ms::numofframe; i++)
		//	{
		//		//use 3 and 4 here;
		//		minkowskisum(i, cirIdx);
		//
		//
		//		offsetMRs	[i] = ms::Model_Result;
		//		offsetMIBs	[i] = ms::ModelInfo_Boundary;
		//		offsetIDC	[i] = ms::InteriorDisks_Convolution;
		//	}
		//}

		// voronoi approach : use clearance info.
		/* for each voronoi, 6 edges are added
			VEdge : v0-v1 is given.
			Construct 4 vertices: 
				v00 : offset v0 toward arc0, with (its_clearance - minimumClearance)
				v01 : offset v0 toward arc1, with (its_clearance - minimumClearance)
				v10 : offset v1 toward arc0, with (its_clearance - minimumClearance)
				v11 : offset v1 toward arc1, with (its_clearance - minimumClearance)
			construct 6 edges:
			v00 - v0 - v01
			 |          |
			v10 - v1 - v11
			
			WARNING!!!:
				(v00, v10) or (v01, v11) might have intersection with csobs if minClr is small.
				(+) line seg connecting start and closest-vor-pt might also have intersection.
		*/
		using lseg = v_edge;
		vector<vector<lseg>> offCso(ms::numofframe); // Def: offset config-space-obs; offCso[sliceNo][segIdx]
		vector<vector<lseg>> offCon(ms::numofframe); // def: offset connector(to voronoi)
		{
			auto& vess = planning::output_to_file::v_edges;
			double minClr = 0.01; // minimum clearance.
			// for (each slice)
			for (int sn = 0; sn < ms::numofframe; sn++) // sn := sliceNo
			{
				auto& ves = vess[sn];
				set<Point> inserted;
				set<int> boundaryCase;

				// for (each ve)
				for  (int i = 0; i < ves.size(); i++)
				{
					auto& ve = ves[i];

					// alias
					auto& v0 = ve.v0;
					auto& v1 = ve.v1;
					auto& arc0 = VRINs[sn].arcs[ve.idx[0]];
					auto& arc1 = VRINs[sn].arcs[ve.idx[1]];
					
					// get direction and clearance
					auto dir00 = arc0.cc() - v0;
					auto dir01 = arc1.cc() - v0;
					auto dir10 = arc0.cc() - v1;
					auto dir11 = arc1.cc() - v1;

					auto dist00 = dir00.length1();
					auto dist01 = dir01.length1();
					auto dist10 = dir10.length1();
					auto dist11 = dir11.length1();

					// if (arc = concave) flip direction
					if (dist00 < arc0.cr()) 
					{
						dir00 = -dir00;
						dir10 = -dir10;
					}
					if (dist01 < arc1.cr())
					{
						dir01 = -dir01;
						dir11 = -dir11;
					}

					// check No. of valid endpoints
					int nValid = 2;
					auto& clr0 = ve.clr0();
					clr0 = fabs(dist00 - arc0.cr());
					if (clr0 < minClr) nValid--;
					auto& clr1 = ve.clr1();
					clr1 = fabs(dist11 - arc1.cr());
					if (clr1 < minClr) nValid--;

					// if both endpoints are valid => we are good to go.
					// else take care
					if (nValid == 0)
						continue;
					else if (nValid == 1)
					{
						// take care of it later
						boundaryCase.insert(i);
						continue;
					}

					// normalize
					dir00 = dir00 / dist00;
					dir01 = dir01 / dist01;
					dir10 = dir10 / dist10;
					dir11 = dir11 / dist11;

					// get vXX
					auto v00 = v0 + (clr0 - minClr) * dir00;
					auto v01 = v0 + (clr0 - minClr) * dir01;
					auto v10 = v1 + (clr1 - minClr) * dir10;
					auto v11 = v1 + (clr1 - minClr) * dir11;

					//if (clr0 - minClr < 0 || clr1 - minClr < 0)
					//	cerr << "!!!!!!!!!!!!!!!!ERR : Direction changed" << endl;

					// offset
					lseg off0;
					off0.v0 = v00;
					off0.v1 = v10;

					lseg off1;
					off1.v0 = v01;
					off1.v1 = v11;

					offCso[sn].push_back(off0);
					offCso[sn].push_back(off1);

					// conn
					if (inserted.find(v0) != inserted.end())
					{
						lseg con00;
						con00.v0 = v0;
						con00.v1 = v00;

						lseg con01;
						con01.v0 = v0;
						con01.v1 = v01;

						offCon[sn].push_back(con00);
						offCon[sn].push_back(con01);

						
					}
					if (inserted.find(v1) != inserted.end())
					{
						lseg con10;
						con10.v0 = v1;
						con10.v1 = v10;

						lseg con11;
						con11.v0 = v1;
						con11.v1 = v11;

						offCon[sn].push_back(con10);
						offCon[sn].push_back(con11);
					}

					// insert those who are done
					inserted.insert(v0);
					inserted.insert(v1);

					//dbg_out
					if (false && sn == 0)
					{
						Point err(-1.04205, -0.392571);
						auto a = (v00 - err).length1();
						auto b = (v01 - err).length1();
						auto c = (v10 - err).length1();
						auto d = (v11 - err).length1();
						if (a < 1e-2 || b < 1e-2 || c < 1e-2 || d < 1e-2)
						{
							cout 
								<< "len :" << endl
								<< a << endl
								<< b << endl
								<< c << endl
								<< d << endl;
							cout
								<< "ver : " << endl
								<< ve.v0 << endl
								<< ve.v1 << endl;
							cout 
								<< "arc : " << endl
								<< arc0 << endl
								<< arc1 << endl;
							cout
								<< "RES: " << endl
								<< v00 << endl
								<< v01 << endl
								<< v10 << endl
								<< v11 << endl
								<< "DIR: " << endl
								<< dir00 << endl
								<< dir01 << endl
								<< dir10 << endl
								<< dir11 << endl
								<< "clr: " << endl
								<< clr0 << endl
								<< clr1 << endl;
						}
					}
				}

				// for (those only that have a single valid endpoint)
				for (auto i : boundaryCase)
				{
					auto& ve = ves[i];

					// alias
					auto& v0 = ve.v0;
					auto& v1 = ve.v1;
					auto& arc0 = VRINs[sn].arcs[ve.idx[0]];
					auto& arc1 = VRINs[sn].arcs[ve.idx[1]];

					// get direction and clearance
					auto dir00 = arc0.cc() - v0;
					auto dir01 = arc1.cc() - v0;
					auto dir10 = arc0.cc() - v1;
					auto dir11 = arc1.cc() - v1;

					auto dist00 = dir00.length1();
					auto dist01 = dir01.length1();
					auto dist10 = dir10.length1();
					auto dist11 = dir11.length1();

					// if (arc = concave) flip direction

					// This part is different from above, as clearance = 0 can be erroneous for this case
					if (dist00 < arc0.cr())
						dir00 = -dir00;
					if (dist10 < arc0.cr())
						dir10 = -dir10;

					if (dist01 < arc1.cr())
						dir01 = -dir01;
					if (dist11 < arc1.cr())
						dir11 = -dir11;

					// find clearance
					auto& clr0 = ve.clr0();
					clr0 = fabs(dist00 - arc0.cr());
					auto& clr1 = ve.clr1();
					clr1 = fabs(dist11 - arc1.cr());

					// we know that one of the clearnace is below min and one is over min.
					double t = (minClr - clr0) / (clr1 - clr0);
					double t1 = 1 - t;

					// Assume: chage of clr is linear in small line seg
					Point vNew = t1 * v0 + t * v1;

					// normalize
					dir00 = dir00 / dist00;
					dir01 = dir01 / dist01;
					dir10 = dir10 / dist10;
					dir11 = dir11 / dist11;

					// get vXX

					auto v00 = v0 + (clr0 - minClr) * dir00;
					auto v01 = v0 + (clr0 - minClr) * dir01;
					auto v10 = v1 + (clr1 - minClr) * dir10;
					auto v11 = v1 + (clr1 - minClr) * dir11;

					// offset : changes with which vert has larger clr
					bool v0larger = clr0 > clr1;

					lseg off0;
					lseg off1;
					if (v0larger)
					{
						off0.v0 = v00;
						off0.v1 = vNew;
						off1.v0 = v01;
						off1.v1 = vNew;
					}
					else
					{
						off0.v0 = vNew;
						off0.v1 = v10;
						off1.v0 = vNew;
						off1.v1 = v11;
					}
					
					offCso[sn].push_back(off0);
					offCso[sn].push_back(off1);

					// conn
					if (v0larger && (inserted.find(v0) != inserted.end()))
					{
						lseg con00;
						con00.v0 = v0;
						con00.v1 = v00;

						lseg con01;
						con01.v0 = v0;
						con01.v1 = v01;

						offCon[sn].push_back(con00);
						offCon[sn].push_back(con01);


					}
					if (!v0larger && (inserted.find(v1) != inserted.end()))
					{
						lseg con10;
						con10.v0 = v1;
						con10.v1 = v10;

						lseg con11;
						con11.v0 = v1;
						con11.v1 = v11;

						offCon[sn].push_back(con10);
						offCon[sn].push_back(con11);
					}

					// insert those who are done
					inserted.insert(v0);
					inserted.insert(v1);

					//dbg_out : problem was aobut
					if (false && sn == 0)
					{
						Point err(-1.04205, -0.392571);
						auto a = (v00 - err).length1();
						auto b = (v01 - err).length1();
						auto c = (v10 - err).length1();
						auto d = (v11 - err).length1();
						if (a < 1e-2 || b < 1e-2 || c < 1e-2 || d < 1e-2)
						{
							cout
								<< "len :" << endl
								<< a << endl
								<< b << endl
								<< c << endl
								<< d << endl;
							cout
								<< "ver : " << endl
								<< ve.v0 << endl
								<< ve.v1 << endl;
							cout
								<< "arc : " << endl
								<< arc0 << endl << arc0.ccw << endl
								<< arc1 << endl << arc1.ccw << endl;
							cout
								<< "RES: " << endl
								<< v00 << endl
								<< v01 << endl
								<< v10 << endl
								<< v11 << endl
								<< "DIR: " << endl
								<< dir00 << endl
								<< dir01 << endl
								<< dir10 << endl
								<< dir11 << endl
								<< "clr: " << endl
								<< clr0 << endl
								<< clr1 << endl;
							cout 
								<< "clr: " << endl
								<< dist00 - arc0.cr() << endl
								<< dist01 - arc1.cr() << endl
								<< dist10 - arc0.cr() << endl
								<< dist11 - arc1.cr() << endl;
							cout
								<< "v1 : " << v1 << endl
								<< "dir11 : " << dir11 << endl
								<< "cl1 : " << clr1 << endl
								<< "minClr: " << minClr << endl
								<< "dir11* a:" << (clr1 - minClr) * dir11 << endl
								<< "v11 : " << v1 + (clr1 - minClr) * dir11 << endl;

							auto norm0 = arc0.x1() - arc0.cc();
							auto norm1 = arc1.x0() - arc1.cc();
							auto t0 = std::atan2(norm0.y(), norm0.x());
							auto t1 = std::atan2(norm1.y(), norm1.x());


							cout << "t0 : " << t0 << endl
								<< "t1 : " << t1 << endl;

							cout << "sn : " << sn << endl
								<< "i : " << i << endl;
							
							//offCso[sn].pop_back();
							//offCso[sn].pop_back();
						}
					}

				}
			}
		}

		// WARNING:
		// simply add offset to vorEdges
		{
			auto& ves = planning::output_to_file::v_edges;
			for (int i = 0; i < ms::numofframe; i++)
			{
				//ves[i].clear();
				ves[i].insert(ves[i].end(), offCon[i].begin(), offCon[i].end());
				ves[i].insert(ves[i].end(), offCso[i].begin(), offCso[i].end());
			}
		}

		auto offsetEndTime = clock();
		cout << "Offset init time : "
			<< double(offsetEndTime - offsetStartTime) / 1000 << "s" << endl;

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

		//Graph theGr = create_VorGraph(v_edges);
		//std::vector<v_edge> path = invoke_AStar(theGr);
		std::vector<std::vector<v_edge>> v_edges_sparced;
		std::vector<std::vector<v_edge>>::iterator iterEdgesLayer = v_edges.begin();
		for ( size_t i = 0; i < v_edges.size(); ++i, ++iterEdgesLayer)
			//if (0 == (i % 10)) // filter out intermediate layers
				v_edges_sparced.push_back(*iterEdgesLayer);

		vector<Vertex> vecVertices;
		map<Vertex, int, VertexLessFn> mapLookup;
		Graph theGr = create_VorGraph(v_edges_sparced, vecVertices, mapLookup);
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
		std::vector<double> renderedPath;
		{
			// just a fake path to check program pipeline.
			//for (int i = 0; i < 100; i++)
			for(auto v : path)
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

		ms::renderPath(argc, argv, renderedPath);

		return 0;
	}

	

	void searchTest()
	{
		double coords2d[9][2] = { {0.,0.}, {1., 0.}, {10., 0.},
								 {0.,1.}, {1., 1.}, {10., 1.},
								 {0.,10.}, {1., 10.}, {10., 10.} };
		int layerEdgeIndeces[12][2] = { {0,1}, {1,2}, {3,4}, {4,5}, {6,7}, {7,8},
										{0,3}, {1,4}, {2,5}, {3,6}, {4,7}, {5,8} };
		vector<v_edge> layer0;
		vector<v_edge> layer1;
		for (int i = 0; i < 12; ++i)
		{
			int idxP = layerEdgeIndeces[i][0];
			int idxQ = layerEdgeIndeces[i][1];
			ms::Point p(coords2d[idxP][0], coords2d[idxP][1]);
			ms::Point q(coords2d[idxQ][0], coords2d[idxQ][1]);
			v_edge currEdge;
			currEdge.v0 = p;
			currEdge.v1 = q;
			currEdge.idx[0] = -1;
			currEdge.idx[1] = -1;
			layer0.push_back(currEdge);
			layer1.push_back(currEdge);
		}
		vector<vector<v_edge>> v_edges;
		v_edges.push_back(layer0);
		v_edges.push_back(layer1);

		vector<Vertex> vecVertices;
		map<Vertex, int, VertexLessFn> mapLookup;
		Graph theGr = create_VorGraph(v_edges, vecVertices, mapLookup);
		Vertex ptnSrc(0, 0, 0);
		Vertex ptnDst(10, 10, 180.);
		std::vector<Vertex> path = invoke_AStar(theGr, vecVertices, mapLookup, ptnSrc, ptnDst);
		if (path.empty())
		{
			cout << "Path not found" << endl;
		}
		else
		{
			for (auto v : path)
			{
				cout << v << endl;
			}
		}

	}
}