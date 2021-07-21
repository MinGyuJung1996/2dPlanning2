#include "MS2D.h"
#include <iostream>
#include <iomanip>
#include "voronoi.hpp"
#include "collision detection.hpp"
#include "AStarOnVorDiag.h"
#include "support.hpp"
#include "dijkstra.hpp"

#define derr cerr
#define COUT(x) #x " : " << x

planning::coneVoronoi coneVor;

extern Point robotForward;

namespace graphSearch
{
	extern std::vector<cd::lineSegmentCollisionTester> testers;
	extern std::vector<cd::pointCollisionTester> testers2; 
	
	using lseg = v_edge;
	extern std::vector<std::vector<std::vector<lseg>>> offCrv; // offCrv[offIdx][sliceNo][lsegIdx];

	extern std::vector<decltype(ms::Model_Result)>				MRs			; // data collected for checking
	extern std::vector<decltype(ms::ModelInfo_Boundary)>		MIBs		; // data collected for checking
	extern std::vector<decltype(ms::InteriorDisks_Convolution)>	IDC			; // save conv/disk
	extern std::vector<planning::VR_IN>							VRINs		;
	extern std::vector<decltype(ms::Model_Result)>				offsetMRs	;	// data collected for checking
	extern std::vector<decltype(ms::ModelInfo_Boundary)>		offsetMIBs	; // data collected for checking
	extern std::vector<decltype(ms::InteriorDisks_Convolution)>	offsetIDC	;  // save conv/disk

	extern djkCalc dijk;
}

void
	refinePathBiarc(vector<Vertex>& in, int nBiarc, vector<CircularArc>& outArc, vector<double>& outRot);
void 
	tessPathBiarc(vector<CircularArc>& in_pathBiarc, vector<double>& in_pathBiarcRot, double in_tessLen, vector<Vertex>& out_pathTess);
void 
	tessPathAlignTheta(vector<Vertex>& pathTess, std::vector<cd::pointCollisionTester>& pTesters, Point forwardDir);
void 
	tessPathClear(vector<Vertex>& in_path, double in_tessLen, vector<Vertex>& out_path);

/*
Def:

Warnings:

	collision test maybe wrong for extremely complicated scenario
	use glo-var nofframe
*/
void optimizePath(vector<xyt>& _in_pos, vector<double>& _in_clr, vector<gs::EdgeType>& _in_ets, vector<xyt>& _out_opt_path)
{
	auto& pos = _in_pos;
	auto& clr = _in_clr;
	auto& ets = _in_ets;
	auto& out = _out_opt_path;

	// 1. find intervals that need optimization
	using Itv = std::pair<int, int>;

	vector<Itv> intervals; // edgeNo, half-open-interval [)
	vector<bool> doOpt;

	if(ets.size() > 0)
	{
		// 1-1. lambda : is this edge optimizable?
		auto isOpt = [](gs::EdgeType et) ->bool
		{
			switch (et)
			{
			// case: do not optimize
			case gs::etEndVer:
			case gs::etStrVer:
			case gs::etStrEnd:
				return false;

			// case: do optimize
			case gs::etOff:
			case gs::etVor:
			case gs::etOcn:
			case gs::etTan:
			case gs::etTcn:
			case gs::etIsc:
			case gs::etStrHor:
			case gs::etEndHor:
			case gs::etUnknow:
				return true;
			}
		};

		// 1-2. find intervals
		int 
			state = isOpt(ets[0]), 
			beg = 0;
		for (int i = 1; i < ets.size(); i++)
		{
			int new_state = isOpt(ets[i]);
			if (new_state != state)
			{
				intervals.push_back(make_pair(beg, i));
				doOpt.push_back(state);

				state = new_state;
				beg = i;
			}
		}
		// take care of last interval;
		intervals.push_back(make_pair(beg, ets.size()));
		doOpt.push_back(state);
	}
	else return;

	// 2. do optimization
	vector<xyt> temp;
	{
		for (int i = 0; i < intervals.size(); i++)
		{
			// 2-1. alias

			// considering edges = [ei0, ei1)
			auto ei0 = intervals[i].first;
			auto ei1 = intervals[i].second; // might not be valid (when ei1 = ets.size())

			// considering verts = [vi0, vi1)
			auto vi0 = ei0;
			auto vi1 = ei1;	// always valid

			// 2-2.
			if (doOpt[i])
			{
				// triv
				if (vi1 - vi0 == 1)
				{
					temp.push_back(pos[i]);
					continue;
				}

				// for optimizing case, first build reduced path and add it to temp
				vector<xyt> reducedPath;

				int start = vi0;
				for (auto i = vi0 + 2; i <= vi1; i++)
				{
					// 2-2-1. check whether reducing vert [start, i) is okay
					bool reduceable = true;
					{
						// samples from straight line
						vector<xyt> lineSamples;
						for (int j = start; j < i; j++)
						{
							// warning : col test maybe unsafe?
							double t = double(j - start) / (i - start);
							xyt lineP = pos[start] * (1-t) + pos[i] * t;
							lineP.t() = pos[j].t();
							lineSamples.push_back(lineP);
						}

						// test collision
						for (int j = 0; j < lineSamples.size(); j++)
						{
							int sliceNo = lineSamples[j].t() / 360.0 * ms::numofframe; // warning!!! numofframe
							Point tested(lineSamples[j].x(), lineSamples[j].y());
							Point closest;
							bool ctRes = gs::testers2[sliceNo].testPrecise(tested, closest);
							if (ctRes || (tested-closest).length2() < 1e-3 )
							{
								reduceable = false;
								break;
							}

						}

					}

					if (reduceable)
					{
						// do nothing. go to next i 
					}
					else
					{
						// we know that [start, i-1) is reduceable => reduce
						vector<xyt> lineSamples;
						for (int j = start; j < (i - 1); j++)
						{
							// warning : col test maybe unsafe?
							double t = double(j - start) / ((i - 1) - start);
							xyt lineP = pos[start] * (1 - t) + pos[i-1] * t;
							lineP.t() = pos[j].t();
							lineSamples.push_back(lineP);
						}
						reducedPath.insert(reducedPath.end(), lineSamples.begin(), lineSamples.end());

						start = i - 1; // notice that next i-value will be i+1. => the edge length is always > =2
					}
				}

				// take care of last edge piece
				if (start == vi1)
				{
					//do nothing
				}
				else if (start == vi1 - 1)
				{
					reducedPath.push_back(pos[start]);
				}
				else // start <= vi1 - 2
				{
					vector<xyt> lineSamples;
					for (int j = start; j < vi1; j++)
					{
						// warning : col test maybe unsafe?
						double t = double(j - start) / (vi1 - start);
						xyt lineP = pos[start] * (1 - t) + pos[vi1] * t;
						lineP.t() = pos[j].t();
						lineSamples.push_back(lineP);
					}
					reducedPath.insert(reducedPath.end(), lineSamples.begin(), lineSamples.end());
				}

				// add to temp;
				temp.insert(temp.end(), reducedPath.begin(), reducedPath.end());
			}
			else // (not  an optimizable edge set)
			{
				for (auto i = vi0; i < vi1; i++)
				{
					temp.push_back(pos[i]);
				}
			}
		}
		// take care of last vert
		temp.push_back(pos.back());
	}

	// 9. swap
	out = temp;
}

namespace ms {

	/* OpenGL Parameters */
	GLsizei wd = 900, ht = 600;
	double zoom = 1.4;
	double tx, ty;
	using namespace std;

	/* Comments - MS2D.h 참고 */

#pragma region global variables related to time.
	int ModelInfo_CurrentFrame = 64; // 355;//355; //148; // 44;		// 7,7,261 case
	int ModelInfo_CurrentFrameOld; // due to the program structure... when we process frame i, modelInfo_CurrentFrame = i + 1...(due to pre increment in some func...) => save val.
	//int ModelInfo_CurrentFrame = 93;			// [0, 360)
	pair<int, int> ModelInfo_CurrentModel;	// [0, 8) x [0, 8)

	int
		&t0 = ModelInfo_CurrentModel.first,
		&t1 = ModelInfo_CurrentModel.second,
		&t2 = ModelInfo_CurrentFrame;	

	int numoftest = 0;						// actually a bool, if 1 -> end test.
#pragma endregion

#pragma region global variables containing data (~~ input for minkowskiSum)
// Models are 1. imported. 2. segmented. 3. rotated.
	vector<BezierCrv>	Models_Imported[8];					// def: result of read_file.							// set@ void initialize()
	vector<Circle>		InteriorDisks_Imported[8];			// def: result of read_file.							// set@ void initialize()

	vector<BezierCrv>	Models[8];							// def: segmented version of Models_Imported.			// set@ initialize	// used@ minkowskiSum
	vector<ArcSpline>	Models_Approx[8];					// def: arc approximated verison of above				// set@ initialize	// used@ minkowskiSum
	vector<BezierCrv>	Models_Rotated[numofframe];			// def: this[i] = i degree rotated First model segments // set@ postProcess // used@ minkowskiSum
	vector<ArcSpline>	Models_Rotated_Approx[numofframe];	// def: arc approximated verison of above				// set@ postProcess // used@ minkowskiSum

	vector<Circle>		InteriorDisks_Rotated[numofframe];	// def: this[i] = i degree rotated disks of First model // set@ postProcess // used@ minkowskiSum
#pragma endregion

#pragma region global variables used during minkowskiSum (~~ intermediate)
	vector<Circle> InteriorDisks_Convolution;				// def: ???												// set&used@ minkowskiSum (and minkowskisum_id)
	Grid Grid_Trimming;										// def: ???												// set@ minkowski // used@ trimmingTest (minkowski->overlapTest->convolution_ArcSplin->)
	CacheCircles Cache_Trimming;							// def: ???												// set@ minkowski // used@ trimmingTest
#pragma endregion

#pragma region global variables storing result (~~ ouput of minkowskiSum)
	vector<deque<ArcSpline>> Model_Result;					// def: resultant minkowskiSum
	vector<bool> ModelInfo_Boundary;						// def: whether each element in model_result is boundary??? (_Result.size == _Boundary)
															//	-> true when outer boundary, not inner loop (i guess...)
#pragma endregion

#pragma region global varibles (misc)
	FILE *f;					// "time.txt"
	bool curves = false;		// not used var.
	bool ModelInfo_Identical;	// flag set true when mod1=mod2 && rotation = 0. -> calls minkowski_id instead of original.
	vector<BezierCrv> Crvr;		// not used.
	clock_t now, last;			// only used in animate_func
#pragma endregion

#define dbg
	dbg int dbgcnt = 0;

	void animate_func()
	{
		now = clock();
		clock_t seconds = now - last;

		// if last_drawn's rotation = 0, show it for 1second.
		if (ModelInfo_CurrentFrame == 1) {
			while (clock() - last < 1000);
		}

		//calculate minkowski sum.
		if ((ModelInfo_CurrentFrame == 0) && (ModelInfo_Identical))
		{
			minkowskisum_id(ModelInfo_CurrentFrame, ModelInfo_CurrentModel.second);
			cout << "minkowski_id\n";
		}
		else
			minkowskisum(ModelInfo_CurrentFrame, ModelInfo_CurrentModel.second);
		
		if(planning::forwardTime)
			ModelInfo_CurrentFrame++;
		dbgcnt = 0;
		//debug 

		//~debug
		

		
		last = now;

		glutPostRedisplay();
	}


	void reshape_callback(GLint w, GLint h)
	{
		if (w * 2 < h * 3) {
			wd = w;
			ht = w * 2 / 3;
			glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
		}
		else {
			wd = h * 3 / 2;
			ht = h;
			glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
		}
	}


	void setup_viewvolume()
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(-zoom + tx, zoom + tx, -zoom + ty, zoom + ty);
	}

	void setup_transform()
	{
	}

	void display_callback(void)
	{
		// if (rotated 360 degree) change model not rotating (=model2)
		if (ModelInfo_CurrentFrame == numofframe) {
			ModelInfo_CurrentFrame = 0;
			ModelInfo_CurrentModel.second++;

			// if (traversed all model2) change model rotating (=model1)
			if (ModelInfo_CurrentModel.second == 8) {
				ModelInfo_CurrentModel.second = 0;
				ModelInfo_CurrentModel.first++;
			}

			// if (traversed all model1) pause screen.
			if (ModelInfo_CurrentModel.first == 8) {
				numoftest++;
				if (numoftest == 1) {
					save();
					system("pause");
					exit(0);
				}
				ModelInfo_CurrentModel.first = ModelInfo_CurrentModel.second = 0;
			}
			postProcess(ModelInfo_CurrentModel.first, ModelInfo_CurrentModel.second);
		}

		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		setup_viewvolume();
		setup_transform();


		// Rendering code....
		glLineWidth(2.0f);
		glPointSize(2.8f);
		glColor3f(0.0f, 0.0f, 0.0f);


		glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
		//glViewport(wd * 1 / 3, 0, 128, 128);

		/* Boundary 위에 Model을 올려놓은 그림을 그리는 코드 */
		//if (curves) {
		//	for (int u = 0; u < (int)Crv[ModelInfo_CurrentModel.second].size(); u++)
		//		Crv[ModelInfo_CurrentModel.second][u].draw();
		//}

		//if (curves) {
		//	glColor3f(1.0f, 0.0f, 0.0f);
		//	for (int i = 0; i < (int)Model_Result.size(); i++)
		//		for (int j = 0; j < (int)Model_Result[i].size(); j++)
		//			for (int u = 0; u < (int)Models_Rotated[ModelInfo_CurrentFrame].size(); u++)
		//				Models_Rotated[ModelInfo_CurrentFrame][u].draw(Model_Result[i][j].Models_Rotated_Approx.front().x[0]);

		//}

		glColor3f(0.0f, 0.0f, 0.0f);

		////////////////////////////////////////////////////////////////////////
		// plannning::output_to_file
		if(OUTPUT_TO_FILE && t2 == 1)
			planning::output_to_file::start();

		vector<vector<CircularArc>>& pushed = planning::output_to_file::ms_obj[t2];
		if(planning::output_to_file::flag)
		for (size_t i = 0, length = Model_Result.size(); i < length; i++)
		{
			if (ModelInfo_Boundary[i])
			{
				vector<CircularArc> v;
				for (auto& a : Model_Result[i])
				{
					for (auto&b : a.Arcs)
						v.push_back(b);
				}
				pushed.push_back(v);
			}
			else continue;
		}

		////////////////////////////////////////////////////////////////////////


		//debug : see order
		static int cnt = 11110;
		cnt++;
		int c = 0;

		// debug : to see if results are connected in order : conclusion : somewhat connected & not?
		if (0 && ModelInfo_CurrentFrame == 62)
		{
			if (Model_Result.size() > 0)
			{
				for (size_t i = 0; i < Model_Result.size(); i++)
				{
					for (size_t j = 0; j < Model_Result[i].size(); j++)
					{
						for (size_t k = 0; k < Model_Result[i][j].Arcs.size(); k++)
						{
							auto& a = Model_Result[i][j].Arcs[k];
							cout << i << ", " << j << ", " << k << " : " << Model_Result[i][j].Arcs[k].x[0].P[0] << ", " << Model_Result[i][j].Arcs[k].x[1].P[0] << ",, " << a.globalccw << " " << a.ccw;
							if (a.lhs)
							{
								cout
									<< ", lhs: r " << a.lhs->c.r << " quad " << a.lhs->n[0] << " ccw " << a.lhs->ccw
									<< ", rhs: r " << a.rhs->c.r << " quad " << a.rhs->n[0] << " ccw " << a.rhs->ccw;
							}
							cout << endl;

						}
						//cout << i << ", " << j << " : " << Model_Result[i][j].Arcs[0].x[0].P[0] << ", " << Model_Result[i][j].Arcs[Model_Result[i][j].Arcs.size()-1].x[1].P[0] << endl;
					}
				}
			}
		}
		//~debug

		// dbg_out
		static int minkInterest = -1;
		if (planning::keyboardflag['2'] && !planning::keyboardflag_last['2']) minkInterest = minkInterest + 1;
		if (planning::keyboardflag['3'] && !planning::keyboardflag_last['3']) minkInterest = minkInterest + 10;
		if (planning::keyboardflag['4'] && !planning::keyboardflag_last['4']) minkInterest = minkInterest + 100;
		if (planning::keyboardflag['5'] && !planning::keyboardflag_last['5']) minkInterest = minkInterest - 1;
		if (planning::keyboardflag['6'] && !planning::keyboardflag_last['6']) minkInterest = minkInterest - 10;
		if (planning::keyboardflag['7'] && !planning::keyboardflag_last['7']) minkInterest = minkInterest - 100;
		cout << "mink interest arc no. :" << minkInterest << endl;
		// ~debug
		
		// DRAW MINK
		glLineWidth(5.0);
		if(planning::drawMinkowski)
		for (int i = 0; i < (int)Model_Result.size(); i++) {

			int arcNo = 0;

			for (int j = 0; j < (int)Model_Result[i].size(); j++)
			{
				for (int k = 0; k < (int)Model_Result[i][j].Arcs.size(); k++)
				{
					dbg if (c > cnt) break;
					dbg c++;

					if (arcNo == minkInterest)
						glColor3d(0.7, 1, 0);
					else
					{
						if (ModelInfo_Boundary[i])
							glColor3d(0.0, 0.0, 0.0);
						else
							glColor3d(0.0, 0.0, 1.0);
					}
					Model_Result[i][j].Arcs[k].draw();
					arcNo++;
				}
				dbg if (c > cnt) break;
			}

			////debug
			//vector<CircularArc> ordered;
			//planning::convertMsOutput_Clockwise(Model_Result[i], ordered);
			//for (size_t i = 0, length = ordered.size(); i < length; i++)
			//{
			//	if (c > cnt) break;
			//	c++;

			//	/*if (i == 97)
			//		glColor3f(0.5f, 0.5f, 1.0f);
			//	if (i == 99)
			//		glColor3f(0, 0, 0);*/
			//	ordered[i].draw();
			//}
		}

		// TEST
		/*glBegin(GL_TRIANGLES);
		glVertex3f(1, 1, -1);
		glVertex3f(0, 1, 0);
		glVertex3f(1, 0, -1);
		glEnd();
*/
		//TEST : conc internal pixel format of opengl can be thought of as unsigned byte(1byte per component, 4byte(rgba) per pixel).
		//	data isn't change if glColor/glReadPixels are done with 1byte-components.
		//	unsigned short changes (255,256,257, 280 are all changed to 257... but why 257?)
		/*{
			static int temp = -3;
			temp+=3;

			glColor3ub(temp, temp+1, temp+2);
			glBegin(GL_TRIANGLES);
			glVertex3d(-100, 100, 0.1);
			glVertex3d(+100, 100, 0.1);
			glVertex3d(0, -300, 0.1);
			glEnd();

			unsigned char b[24];
			glReadPixels(300, 300, 2, 3, GL_RGBA, GL_UNSIGNED_BYTE, b);
			for (auto i : b)
				cerr << int(i) << " ";
			cerr << endl;
			if (b[0] != temp) cerr << "eeeeeeee: " << b[0] << " " << temp << endl;
			if (b[1] != temp+1) cerr << "eeeeeeee: " << b[1] << " " << temp+ 1 << endl;
			if (b[2] != temp + 2) cerr << "eeeeeeee: " << b[2] << " " << temp + 2 << endl;
		}*/

		/////////////////////////////////////////////////////////////////////////
		// DRAW BOUND

		if(planning::drawBoundary)
		{
			auto temp = RES;
			RES = 100;
			double r = 2;
			//CircularArc a(Point(0, 0), r, Point(+1, +0), Point(+0, +1)); a.draw();
			//CircularArc b(Point(0, 0), r, Point(+0, +1), Point(-1, +0)); b.draw();
			//CircularArc c(Point(0, 0), r, Point(-1, +0), Point(+0, -1)); c.draw();
			//CircularArc d(Point(0, 0), r, Point(+0, -1), Point(+1, +0)); d.draw();
			for (auto& arc : planning::voronoiBoundary)
			{
				arc.draw();
			}
			RES = temp;
		}

		glColor3f(0.0f, 0.7f, 0.0f);
		static int vr_cnt = 0;
		vr_cnt++;
		
		/*if(vr_cnt > 2)
		{*/
			planning::VR_IN vrin;
			planning::_Convert_MsOut_To_VrIn(Model_Result, ModelInfo_Boundary, vrin);
			planning::_Medial_Axis_Transformation(vrin);
			cout << "---------------------------------------------number of circular arcs ("<< ModelInfo_CurrentModel.first << ", " <<ModelInfo_CurrentModel.second << ", " << ModelInfo_CurrentFrame << ") : " << vrin.arcs.size() << endl;
			cout << "---------------------------------------------number of voronoi line segments : " << planning::lineSegCnt << endl;

			/////////////////////////////////////////////////////////////////////////
			//output_to_file
			if (planning::output_to_file::flag)
				planning::output_to_file::vrIn[t2] = vrin;

			/////////////////////////////////////////////////////////////////////////
			//draw coneVoronoi
			coneVor.colorType = 1;
			if (planning::keyboardflag['q'])
				coneVor.drawVoronoi(vrin);
			/////////////////////////////////////////////////////////////////////////

			// dbg
			if (1)
			{
				// tag0544
				for (auto a : vrin.arcs)
				{
					Point p(0.595017, -0.75676);
					double eps = 1e-8;
					if ((a.x[0] - p).length() < eps || (a.x[1] - p).length() < eps)
					{
						cerr << "normals " << a.n[0] << "       " << a.n[1] << endl;
						a.draw();
					}
				}
			}
		/*}*/
		glColor3f(0, 0, 0);
		/////////////////////////////////////////////////////////////////////////

		//DRAW circlesToDraw;
		if(planning::keyboardflag['e'])
		for (auto a : planning::circlesToDraw)
			a.draw();

		/////////////////////////////////////////////////////////////////////////

		// debug : to see if results are connected in order : conclusion : somewhat connected & not?
		if (0 && ModelInfo_CurrentFrame == 1)
		{
			if (Model_Result.size() > 0)
			{
				for (size_t i = 0; i < Model_Result.size(); i++)
				{
					for (size_t j = 0; j < Model_Result[i].size(); j++)
					{
						for (size_t k = 0; k < Model_Result[i][j].Arcs.size(); k++)
						{
							auto& a = Model_Result[i][j].Arcs[k];
							cout << i << ", " << j << ", " << k << " : " << Model_Result[i][j].Arcs[k].x[0].P[0] << ", " << Model_Result[i][j].Arcs[k].x[1].P[0] << ",, " << a.globalccw << " " << a.ccw;
							if (a.lhs)
							{
								cout
									<< ", lhs: r " << a.lhs->c.r << " quad " << a.lhs->n[0] << " ccw " << a.lhs->ccw
									<< ", rhs: r " << a.rhs->c.r << " quad " << a.rhs->n[0] << " ccw " << a.rhs->ccw;
							}
							cout << endl;
						}
						//cout << i << ", " << j << " : " << Model_Result[i][j].Arcs[0].x[0].P[0] << ", " << Model_Result[i][j].Arcs[Model_Result[i][j].Arcs.size()-1].x[1].P[0] << endl;
					}
				}
			}
		}
		//~debug

		////////////////////////////////////////////////////////////////////////
		// plannning::output_to_file
		if (planning::output_to_file::flag && t2 == 0)
			planning::output_to_file::end();

		////////////////////////////////////////////////////////////////////////
		//minksum approx
		if (planning::keyboardflag['c'])
		{
			glClearColor(1.0, 1.0, 1.0, 1.0);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			setup_viewvolume();
			setup_transform();

			// Rendering code....
			glLineWidth(2.0f);
			glPointSize(2.8f);


			glColor3f(1.0f, 0.0f, 0.0f);
			for (auto& as : Models_Approx[t1])
			{
				auto translation = as.init();
				for(auto& as : Models_Rotated_Approx[t2])
					for (auto& arc : as.Arcs)
					{
						cd::translateArc(arc, translation).draw();
					}
			}

			glColor3f(0.0f, 0.0f, 1.0f);
			for (auto& as : Models_Approx[t1])
				as.draw();
		}

		////////////////////////////////////////////////////////////////////////

		//interest arcspline
		//debug
		static int interest = 0, interest2 = 0, cpressed = -1, vpressed = -1;
		{
			if (planning::keyboardflag['c'] && !cpressed)
			{
				interest++;
				cpressed = 1;
			}
			if (!planning::keyboardflag['c']) cpressed = 0;
			if (planning::keyboardflag['v'] && !vpressed)
			{
				interest2++;
				vpressed = 1;
			}
			if (!planning::keyboardflag['v']) vpressed = 0;
			//interest = interest % Models_Rotated_Approx[ModelInfo_CurrentFrame].size();
			//interest2 = interest2 % Models_Approx[ModelInfo_CurrentModel.second].size();
		}

		//dbg_out 
		cout << "interest ArcSpline index 1 & 2 : " << interest << "   " << interest2 << endl;

		// debug time(1,1,257), interest (5, 4)

		////////////////////////////////////////////////////////////////////////
		/* Left View ports */
		////////////////////////////////////////////////////////////////////////
		glColor3f(0.0f, 0.0f, 0.0f);

		if (false) // true: draw original, false : draw Approx version
		{
			glViewport(0, 0, wd * 1 / 3, ht / 2);
			for (int i = 0; i < (int)Models_Rotated[ModelInfo_CurrentFrame].size() && i < cnt; i++)
				Models_Rotated[ModelInfo_CurrentFrame][i].draw();

			glViewport(0, ht / 2, wd * 1 / 3, ht / 2);
			for (int i = 0; i < (int)Models_Imported[ModelInfo_CurrentModel.second].size() && i < cnt; i++) {
				Models_Imported[ModelInfo_CurrentModel.second][i].draw();
			}
		}
		else
		{
			glViewport(0, 0, wd * 1 / 3, ht / 2);
			if (planning::keyboardflag['x'])
			{
				for (auto& as : Models_Approx[ModelInfo_CurrentModel.second])
					for (auto& arc : as.Arcs)
					{
						auto a2 = arc;
						a2.ccw = (arc.n[0] ^ arc.n[1]) > 0; // tag maccw
						a2.draw2();
					}
			}
			else
			{
				for (int i = 0; i < (int)Models_Rotated_Approx[ModelInfo_CurrentFrame].size() && i < cnt; i++)
				{
					if (i == interest)
						glColor3f(0, 1, 1);
					else
						glColor3f(0, 0, 0);
					Models_Rotated_Approx[ModelInfo_CurrentFrame][i].draw();

				}
				if (planning::keyboardflag['z'])
					for (auto& idr : InteriorDisks_Rotated[t2])
						idr.draw();
			}
			glViewport(0, ht / 2, wd * 1 / 3, ht / 2);
			if (planning::keyboardflag['x']) // if (x) draw sweep
			{
				vector<CircularArc> temp;
				for (auto& as : Models_Approx[ModelInfo_CurrentModel.second])
					for (auto& arc : as.Arcs)
					{
						// arc.
						auto a2 = arc;
						a2.n0() = (a2.x0() - a2.c.c).normalize();
						a2.n1() = (a2.x1() - a2.c.c).normalize();
						a2.ccw = (a2.n[0] ^ a2.n[1]) > 0;
						a2.convex = a2.ccw;

						temp.push_back(a2);


						////dbg_out (tag maccw)
						//bool ccwEqual = ((arc.n[0] ^ arc.n[1]) > 0 == arc.ccw);
						//if (!ccwEqual) cout << "ccw not equal, cross, ccw :" << asin(arc.n[0] ^ arc.n[1])*180/3.1415 << " " << arc.ccw << endl;
					}
				
				////test
				////result : arc[n].x1 = arc[n+1].x0 in model_approx
				//double accumulate = 0;
				//for (size_t i = 0, length = temp.size(); i < length; i++)
				//{
				//	auto arc0 = temp[i];
				//	auto arc1 = temp[i + 1];
				//	if (i != length - 1)
				//	{
				//		auto len = (arc0.x1() - arc1.x0()).length2();
				//		//dbg_out 
				//		cout << arc0.x1() << " " << arc1.x0() << endl;
				//		accumulate += len;
				//	}
				//	temp[i].convex = !temp[i].ccw;
				//}
				////dbg_out
				//cout << "accumulate : " << accumulate << endl;

				// dbg
				/*temp.resize(0);
				double r = 0.5;*/
				//CircularArc a(Point(0, 1), r, Point(+1, +0), Point(+0, +1)); temp.push_back(a);
				////CircularArc b(Point(-0.5, 1.5), -r, Point(+0, +1), Point(-1, +0)); b.ccw = true; temp.push_back(b);
				//CircularArc b(Point(0, 1), r, Point(+0, +1), Point(-1, +0)); temp.push_back(b);
				//CircularArc c(Point(0, 1), r, Point(-1, +0), Point(-0, -1)); temp.push_back(c);
				//CircularArc d(Point(0, 1), r, Point(+0, -1), Point(+1, +0)); temp.push_back(d);

				//// convex case
				//CircularArc a(Point(0, 0.4) + Point( 0.5,  0.5), -r, Point(+1, +0), Point(+0, +1)); a.ccw = true; a.convex = false; temp.push_back(a);
				//CircularArc b(Point(0, 0.4) + Point(-0.5,  0.5), -r, Point(+0, +1), Point(-1, +0)); b.ccw = true; b.convex = false; temp.push_back(b);
				//CircularArc c(Point(0, 0.4) + Point(-0.5, -0.5), -r, Point(-1, +0), Point(-0, -1)); c.ccw = true; c.convex = false; temp.push_back(c);
				//CircularArc d(Point(0, 0.4) + Point( 0.5, -0.5), -r, Point(+0, -1), Point(+1, +0)); d.ccw = true; d.convex = false; temp.push_back(d);

				//// dbg
				//temp.resize(0);
				//double r = 0.5;
				//CircularArc a(Point(0, 0.6) + Point(+0.4, 0), r, Point(-0.8, +0.6), Point(-0.8, -0.6)); a.ccw = true; a.convex = true; temp.push_back(a);
				//CircularArc b(Point(0, 0.6) + Point(-0.4, 0), r, Point(+0.8, -0.6), Point(+0.8, +0.6)); b.ccw = true; b.convex = true; temp.push_back(b);

				auto sweep = cd::getRotationBoundary(temp, double(ms::t2));
				//dbg

				// dbg_out
				cout << sweep.size() << " " << temp.size() << endl;
				for (auto i : sweep)
					i.draw2();
			}
			else //else draw original
			{
				for (int i = 0; i < (int)Models_Approx[ModelInfo_CurrentModel.second].size() && i < cnt; i++) 
				{
					if (i == interest2)
						glColor3f(0, 1, 1);
					else
						glColor3f(0, 0, 0);
				
					Models_Approx[ModelInfo_CurrentModel.second][i].draw();
				}
			if (planning::keyboardflag['z']) //if(z) draw internal disk
				for (auto& id : InteriorDisks_Imported[ModelInfo_CurrentModel.second])
					id.draw();
			}
		}
		
		// keyboardflag
		for (size_t i = 0; i < 256; i++)
		{
			planning::keyboardflag_last[i] = planning::keyboardflag[i];

		}
		glutSwapBuffers();
	}

	void hit_index(int x, int y)
	{
	}

	Point clickedPoint;
	vector<Point> clickedPoints;
	void mouse_callback(int button, int action, int x, int y)
	{

		if (button == 3)
			zoom *= 1.05;

		if (button == 4)
			zoom *= 0.95;

		if (button == GLUT_LEFT_BUTTON && action == GLUT_UP)
		{
			auto t0 = clock();
			auto fx = float(x - 2 * wd / 3) / (wd/3);
			auto fy = float(-y + ht / 2) / (ht / 2);
			fx = zoom * fx + tx;
			fy = zoom * fy + ty;
			cerr
				<< "***************************************************************" << endl
				<< "clicked point : " << fx << ", " << fy << "   ( zoom = " << zoom << " ) " << endl
				<< "***************************************************************" << endl;
			clickedPoint.x() = fx;
			clickedPoint.y() = fy;

			clickedPoints.push_back(clickedPoint);
			while ((clock() - t0) < CLOCKS_PER_SEC)
			{

			}
		}

		static int grb = 0;
		static Point c;
		if (button == GLUT_RIGHT_BUTTON)
		{
			auto fx = float(x - 2 * wd / 3) / (wd / 3);
			auto fy = float(-y + ht / 2) / (ht / 2);
			fx = zoom * fx + tx;
			fy = zoom * fy + ty;

			if (grb % 2 == 0)
			{
				c = Point(fx, fy);
			}
			else
			{
				planning::circlesToDraw.push_back(Circle(c, Point(fx, fy)));
			}

			grb++;
		}

		glutPostRedisplay();
	}

	void keyboard_callback(unsigned char a, int b, int c) {
		/*
		currently used :
		w
		a s d f g
		z x c v b n m
		*/

		using namespace planning;
		if (a == 'f')
			forwardTime = (forwardTime + 1) % 2;
		if (a == 'v')
			drawVoronoiSingleBranch = (drawVoronoiSingleBranch + 1) % 2;
		if (a == 'b')
			drawBifurCircle = (drawBifurCircle + 1) % 2;
		if (a == 'm')
			drawMinkowski = (drawMinkowski + 1) % 2;
		if (a == 'n')
			drawBoundary = (drawBoundary + 1) % 2;
		if (a == 'g')
			drawTransition = (drawTransition + 1) % 2;

		keyboardflag[a] = !keyboardflag[a];

		if (a == 'r')
		{
			rfbTerminationEps *= 10.0;
			cerr << "rfbTerminationEps = " << rfbTerminationEps << endl;
		}
		if (a == 't')
		{
			rfbTerminationEps /= 10;
			cerr << "rfbTerminationEps = " << rfbTerminationEps << endl;
		}


		if (a == 'a')
			tx += 0.05 * zoom;
		if (a == 'd')
			tx -= 0.05 * zoom;
		if (a == 's')
			ty += 0.05 * zoom;
		if (a == 'w')
			ty -= 0.05 * zoom;
		//if (a == 'm')
		//	curves ^= true;
		glutPostRedisplay();
	}

	int main(int argc, char *argv[])
	{
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
		glutInitWindowSize(wd, ht);
		glutCreateWindow("Minkowski Sum");
		initialize();
		glutReshapeFunc(reshape_callback);
		glutDisplayFunc(display_callback);
		glutMouseFunc(mouse_callback);
		glutKeyboardFunc(keyboard_callback);
		glutIdleFunc(animate_func);

		glEnable(GL_DEPTH_TEST);

		glutMainLoop();
		return 0;
	}


	/*
	main func for testing RSV
	*/
	int main2(int argc, char* argv[])
	{
		wd = 1200, ht = 600;

		auto reshapeFunc = [](GLint w, GLint h)->void
		{
			if (w  < h * 2) {
				wd = w;
				ht = w / 2;
				//glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
			}
			else {
				wd = h * 2;
				ht = h;
				//glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
			}
		};
		auto displayFunc = [](void) -> void
		{
			
			glClearColor(1.0, 1.0, 1.0, 1.0);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			setup_viewvolume();
			setup_transform();

			// Rendering code....
			glLineWidth(2.0f);
			glPointSize(2.8f);
			glColor3f(0.0f, 0.0f, 0.0f);

			// build temp
			vector<CircularArc> temp;
			{
				for (auto& as : Models_Approx[ModelInfo_CurrentModel.second])
					for (auto& arc : as.Arcs)
					{
						// arc.
						auto a2 = arc;
						a2.n0() = (a2.x0() - a2.c.c).normalize();
						a2.n1() = (a2.x1() - a2.c.c).normalize();
						a2.ccw = (a2.n[0] ^ a2.n[1]) > 0;
						a2.convex = a2.ccw;

						temp.push_back(a2);


						////dbg_out (tag maccw)
						//bool ccwEqual = ((arc.n[0] ^ arc.n[1]) > 0 == arc.ccw);
						//if (!ccwEqual) cout << "ccw not equal, cross, ccw :" << asin(arc.n[0] ^ arc.n[1])*180/3.1415 << " " << arc.ccw << endl;

						
					}

				////test
				////result : arc[n].x1 = arc[n+1].x0 in model_approx
				//double accumulate = 0;
				//for (size_t i = 0, length = temp.size(); i < length; i++)
				//{
				//	auto arc0 = temp[i];
				//	auto arc1 = temp[i + 1];
				//	if (i != length - 1)
				//	{
				//		auto len = (arc0.x1() - arc1.x0()).length2();
				//		//dbg_out 
				//		cout << arc0.x1() << " " << arc1.x0() << endl;
				//		accumulate += len;
				//	}
				//	temp[i].convex = !temp[i].ccw;
				//}
				////dbg_out
				//cout << "accumulate : " << accumulate << endl;

				// dbg
				/*temp.resize(0);
				double r = 0.5;*/
				//CircularArc a(Point(0, 1), r, Point(+1, +0), Point(+0, +1)); temp.push_back(a);
				////CircularArc b(Point(-0.5, 1.5), -r, Point(+0, +1), Point(-1, +0)); b.ccw = true; temp.push_back(b);
				//CircularArc b(Point(0, 1), r, Point(+0, +1), Point(-1, +0)); temp.push_back(b);
				//CircularArc c(Point(0, 1), r, Point(-1, +0), Point(-0, -1)); temp.push_back(c);
				//CircularArc d(Point(0, 1), r, Point(+0, -1), Point(+1, +0)); temp.push_back(d);

				//// convex case
				//CircularArc a(Point(0, 0.4) + Point( 0.5,  0.5), -r, Point(+1, +0), Point(+0, +1)); a.ccw = true; a.convex = false; temp.push_back(a);
				//CircularArc b(Point(0, 0.4) + Point(-0.5,  0.5), -r, Point(+0, +1), Point(-1, +0)); b.ccw = true; b.convex = false; temp.push_back(b);
				//CircularArc c(Point(0, 0.4) + Point(-0.5, -0.5), -r, Point(-1, +0), Point(-0, -1)); c.ccw = true; c.convex = false; temp.push_back(c);
				//CircularArc d(Point(0, 0.4) + Point( 0.5, -0.5), -r, Point(+0, -1), Point(+1, +0)); d.ccw = true; d.convex = false; temp.push_back(d);

				//// dbg
				//temp.resize(0);
				//double r = 0.5;
				//CircularArc a(Point(0, 0.6) + Point(+0.4, 0), r, Point(-0.8, +0.6), Point(-0.8, -0.6)); a.ccw = true; a.convex = true; temp.push_back(a);
				//CircularArc b(Point(0, 0.6) + Point(-0.4, 0), r, Point(+0.8, -0.6), Point(+0.8, +0.6)); b.ccw = true; b.convex = true; temp.push_back(b);
			}

			// build internal circles
			double r = 0.2 - 1e-8;
			Point
				center0(0, 0),
				center1(0, 0.3),
				center2(0, -0.3);
			Circle
				circle0(center0, r),
				circle1(center1, r),
				circle2(center2, r);

			//dbg_out : n of conv/conc
			int conv = 0, conc = 0;
			for (auto& a : temp)
				if (a.convex)
					conv++;
				else
					conc++;
			cout << " temp: No of conv/conv : " << conv << "   " << conc << endl;
			// ~dbg_out

			// dbg : test info
			//ms::t2 = 1;
			cout << "****** Model : " << ModelInfo_CurrentModel.second << "   /   Rotation(degree) :  " << ms::t2 << endl;
			// ~dbg

			// dbg : write rectangle to file
			static bool write = false;
			ofstream aR ;
			ofstream aRR;
			ofstream cRR;

			if (write)
			{
				aR.open ("arcModel.txt");
				aRR.open("arcModelRSV.txt");
				cRR.open("circModelRSV.txt");


				auto writeArc = [&](CircularArc& c)	//write arc to file
				{
					aR << scientific << setprecision(20);
					aR << c.c.c.P[0] << " " << c.c.c.P[1] << " " << c.c.r << " " << atan2(c.n[0].P[1], c.n[0].P[0]) << " " << atan2(c.n[1].P[1], c.n[1].P[0]) << " " << c.ccw << endl;
				};

				aR << temp.size() << endl;
				for (auto arc : temp)
				{
					writeArc(arc);
				}
			}
			if (write)
			{
				auto writeCircle = [&](Circle& c)	//write arc to file
				{
					cRR << fixed << setprecision(20);
					cRR << c.c.P[0] << " " << c.c.P[1] << " " << c.r << endl;
				};

				cRR << 3 << endl;
				writeCircle(circle0);
				writeCircle(circle1);
				writeCircle(circle2);
			}
			// ~dbg

			// VIEWPORT 0
			glViewport(0, 0, wd / 2, ht);
			{
				auto sweep = cd::getRotationBoundary(temp, double(ms::t2), 1);
				for (size_t i = 0, length = sweep.size(); i < length; i++)
				{
					//debug
					if (sweep[i].convex)
						glColor3f(0, 0, 0);
					else
						glColor3f(0, 1, 0);

					if (i == 0)
					{
						glColor3f(0, 1, 0);
						cout << "interest arc endpoint : " << sweep[i].x0() << "    " << sweep[i].x1() << endl;
					}
					else
						glColor3f(0, 0, 0);
					sweep[i].draw2(0.5);

				}

				////dbg_out : draw circ
				//circle0.draw();
				//circle1.draw();
				//circle2.draw();

				// dbg_out
				cout << "before trimming size : " << sweep.size() << endl;
			}

			// VIEWPORT 1
			glViewport(wd / 2, 0, wd / 2, ht);
			{
				auto sweep2 = cd::getRotationBoundary(temp, double(ms::t2), 0);
				// dbg_out
				cout << "divArc size : " << sweep2.size() << endl;
				glColor3f(0, 0, 0);
				glLineWidth(2.0f);
				//debug
				static int interest = 0, cpressed = 0, vpressed = 0;
				{
					if (planning::keyboardflag['c'] != cpressed)
					{
						interest++;
						cd::trsiInterestIdx++;
					}
					cpressed = planning::keyboardflag['c'];
					if (planning::keyboardflag['v'] != vpressed)
					{
						interest--;
						cd::trsiInterestIdx--;
					}
					vpressed = planning::keyboardflag['v'];

					if(sweep2.size() > 0)
					interest = interest % sweep2.size();
				}
				for (size_t i = 0, length = sweep2.size(); i < length; i++)
				{
					////debug
					//if (sweep2[i].convex)
					//	glColor3f(0, 0, 0);
					//else
					//	glColor3f(0, 1, 0);

					// debug
					{

						if (i == interest)
						{
							glColor3f(0, 1, 0);
							cout << "interest " << interest << " : interest arc endpoint : " << sweep2[i].x0() << "    " << sweep2[i].x1() << endl;
						}
						else
							glColor3f(0, 0, 0);
					}
					sweep2[i].draw2(0.5);
				}

				if (planning::keyboardflag['x'])
				{
					glColor3f(0, 0, 1);
					for (auto as : Models_Approx[ModelInfo_CurrentModel.second])
						for (auto a : as.Arcs)
							a.draw();
				}

				// dbg : write rectangle
				if (write)
				{
					auto writeArc = [&](CircularArc& c)	//write arc to file
					{
						aRR << fixed << setprecision(20);
						aRR << c.c.c.P[0] << " " << c.c.c.P[1] << " " << c.c.r << " " << atan2(c.n[0].P[1], c.n[0].P[0]) << " " << atan2(c.n[1].P[1], c.n[1].P[0]) << " " << c.ccw << endl;
					};

					// build list of indices of arcs to use.
					set<int> untrimmed;
				
					untrimmed.insert(1);
					untrimmed.insert(3);
					untrimmed.insert(6);
					untrimmed.insert(8);
					untrimmed.insert(11);
					untrimmed.insert(13);
					untrimmed.insert(16);
					untrimmed.insert(18);

					for (int i = 20; i < 32; i++)
						untrimmed.insert(i);

					// fout
					aRR << untrimmed.size() << endl;
					for (int i = 0; i < 32; i++)
						if (untrimmed.find(i) != untrimmed.end())
							writeArc(sweep2[i]);
				}
				// ~dbg
			}

			// dbg
			if (write)
			{
				write = false; // this is last write;
				aR.close();
				aRR.close();
				cRR.close();
			}

			// ~dbg

			glutSwapBuffers();
		};
		auto    idleFunc = []() -> void
		{
			clock_t now = clock();

			// if last_drawn's rotation = 0, show it for 1second.
			if (ModelInfo_CurrentFrame == numofframe - 1) {
				while (clock() - now < 500);
			}
			
			if (planning::forwardTime)
			ModelInfo_CurrentFrame++;
			if (ModelInfo_CurrentFrame == numofframe)
			{
				ModelInfo_CurrentModel.second = (ModelInfo_CurrentModel.second + 1) % 8;
				ModelInfo_CurrentFrame = 0;
			}
			
			// debug
			CircularArc test;
			test.convex = false;
			CircularArc B(test);
			CircularArc C;
			C = test;
			cout << "convexity test " << test.convex << B.convex << C.convex << endl;
			// ~debug
			
			glutPostRedisplay();
		};


		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
		glutInitWindowSize(wd, ht);
		glutCreateWindow("RSV");
		initialize();
		t2 = t1 = t0 = 0;
		t1 = 2;
		t0 = 90;
		//override model0
		{
			double
				r0 = 10,
				r1 = 0.1;

			Point // from 1st quadrant, in ccw dir
				p0(0.2, 0.5),
				p1(-0.2, 0.5),
				p2(-0.2, -0.5),
				p3(0.2, -0.5);
			Point
				c0(0, -r0),
				c1(+r0, 0),
				c2(0, +r0),
				c3(-r0, 0);

			vector<CircularArc>
				vec;

			// semi-straight arcs
			vec.push_back(cd::constructArc(c0, p0 - Point(r1, 0), p1 + Point(r1, 0), true));
			vec.push_back(cd::constructArc(c1, p1 - Point(0, r1), p2 + Point(0, r1), true));
			vec.push_back(cd::constructArc(c2, p2 + Point(r1, 0), p3 - Point(r1, 0), true));
			vec.push_back(cd::constructArc(c3, p3 + Point(0, r1), p0 - Point(0, r1), true));

			// corner-arcs
			vec.push_back(cd::constructArc(p0 + Point(-r1, -r1), r1, PI * 0.0, PI * 0.5));
			vec.push_back(cd::constructArc(p1 + Point(+r1, -r1), r1, PI * 0.5, PI * 1.0));
			vec.push_back(cd::constructArc(p2 + Point(+r1, +r1), r1, PI * 1.0, PI * 1.5));
			vec.push_back(cd::constructArc(p3 + Point(-r1, +r1), r1, PI * 1.5, PI * 2.0));

			for (auto& a : vec)
				a.convex = true;

			ArcSpline as;
			as.Arcs = vec;
			vector<ArcSpline> asvec;
			asvec.push_back(as);

			Models_Approx[0] = asvec;

			//debug
			auto
				t0 = p2 - c2,
				t2 = p3 - c2,
				t1 = Point(0, 0) - c2;
			auto
				th0 = atan2(t0.y(), t0.x()),
				th1 = atan2(t1.y(), t1.x()),
				th2 = atan2(t2.y(), t2.x());

			cout << "Thetas in square " << th0 << " " << th1 << " " << th2 << endl;
		}

		glutReshapeFunc(reshapeFunc);
		glutDisplayFunc(displayFunc);
		glutMouseFunc(mouse_callback);
		glutKeyboardFunc(keyboard_callback);
		glutIdleFunc(idleFunc);

		glEnable(GL_DEPTH_TEST);

		glutMainLoop();
		return 0;
	}


	namespace renderMinkVoronoiGlobal
	{
		// as GLUT needs global variables...
		// not the best implementaion, but these part will only be used for testing/error-checking

		std::vector<decltype(ms::Model_Result)> MRs;
		std::vector<decltype(ms::ModelInfo_Boundary)> MIBs;
		std::vector<std::vector<planning::output_to_file::v_edge>> v_edges;
		decltype(planning::voronoiBoundary) voronoiBoundary;
		int& sliceIdx = ms::t2;

	}
	/*
	similar struct to main funcs

	Def: when called, render minkowskisum and voronoi diagram
	Just to check calculation's correctness.
	*/
	void renderMinkVoronoi(
		int argc,
		char* argv[],
		std::vector<decltype(ms::Model_Result)>& MRs,
		std::vector<decltype(ms::ModelInfo_Boundary)>& MIBs,
		std::vector<std::vector<planning::output_to_file::v_edge>>& v_edges,
		decltype(planning::voronoiBoundary)& voronoiBoundary
		)
	{
		wd = 1200, ht = 800;

		//inefficiency... copy data
		renderMinkVoronoiGlobal::MRs = MRs;
		renderMinkVoronoiGlobal::MIBs = MIBs;
		renderMinkVoronoiGlobal::v_edges = v_edges;
		renderMinkVoronoiGlobal::voronoiBoundary = voronoiBoundary;


		auto reshapeFunc = [](GLint w, GLint h) ->void
		{
			double ratioWH = double(wd)/ht;
			if (w < h * ratioWH) {
				wd = w;
				ht = w / ratioWH;
				//glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
			}
			else {
				wd = h * ratioWH;
				ht = h;
				//glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
			}
		};
		auto    idleFunc = []() -> void
		{
			//manage speed
			static clock_t lastTime = clock();
			double timeInterval = 1000.0 / 36; // 10.0 * 1000 / ms::numofframe; //milliseconds
			while (clock() - lastTime < timeInterval)
			{

			}
			lastTime = clock();

			// managing input
			if(planning::keyboardflag['f'])	// toggle with f
				renderMinkVoronoiGlobal::sliceIdx++;

			static bool lastG = false, lastH = false;
			if (planning::keyboardflag['g'] != lastG)
				renderMinkVoronoiGlobal::sliceIdx++;
			if (planning::keyboardflag['h'] != lastH)
				renderMinkVoronoiGlobal::sliceIdx--;
			lastG = planning::keyboardflag['g'];
			lastH = planning::keyboardflag['h'];

			if (renderMinkVoronoiGlobal::sliceIdx >= ms::numofframe) renderMinkVoronoiGlobal::sliceIdx = 0;
			if (renderMinkVoronoiGlobal::sliceIdx < 0 ) renderMinkVoronoiGlobal::sliceIdx = ms::numofframe - 1;

			cout << "Current slice Number : " << renderMinkVoronoiGlobal::sliceIdx << endl;
			glutPostRedisplay();
		};
		auto displayFunc = [](void) -> void
		{
			using namespace renderMinkVoronoiGlobal;
			glClearColor(1.0, 1.0, 1.0, 1.0);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			setup_viewvolume();
			setup_transform();

			// Rendering code....
			glLineWidth(2.0f);
			glPointSize(2.8f);
			glColor3f(0.0f, 0.0f, 0.0f);

			// build temp

			auto sceneIdx = 7;
			auto robotIdx = 1;
			
			// VIEWPORT 0 (Upper-Left) : draw scene
			glViewport(0, 0, wd * 1 / 3, ht / 2);
			auto& scene = Models_Approx[sceneIdx];
			glColor3f(0, 0, 0);
			for (auto& as : scene)
				for (auto& arc : as.Arcs)
					arc.draw();
			
			// VIEWPORT 1 (Lower-Left) : darw robot
			glViewport(0, ht / 2, wd * 1 / 3, ht / 2);
			auto& robot = Models_Approx[robotIdx];
			glColor3f(0, 0, 0);
			for(auto& as: robot)
				for (auto& arc : as.Arcs)
				{
					cd::rotateArc(arc, sliceIdx * 360.0 / ms::numofframe).draw();
				}
			
			// VIEWPORT 2 (RIGHT) : draw mink/voronoi
			glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
			auto& MR		= renderMinkVoronoiGlobal::MRs[sliceIdx];
			auto& MIB		= renderMinkVoronoiGlobal::MIBs[sliceIdx];
			auto& voronoi	= renderMinkVoronoiGlobal::v_edges[sliceIdx];
			auto& boundary	= renderMinkVoronoiGlobal::voronoiBoundary;

			static vector<CircularArc> ch;
			static vector<CircularArc> chOff;

			// 1. draw mink
			if(!planning::keyboardflag['1'])
			{
				glColor3f(0, 0, 0);
				size_t length = MR.size();
				for (size_t i = 0; i < length; i++)
				{
					//if (MIB[i])
					//	glColor3f(0, 0, 0);
					//else
					//	glColor3f(0, 0, 1);

					for (auto& as : MR[i])
						as.draw();
				}
			}

			// 2. draw voronoi
			if (!planning::keyboardflag['2'])
			{
				glColor3f(0, 0, 1);
				glBegin(GL_LINES);
				for (auto& ve : voronoi)
				{
					glVertex2dv(ve.v0.P);
					glVertex2dv(ve.v1.P);
				}
				glEnd();
			}

			// 3. draw Boundary for voronoi edges
			glColor3f(0, 0, 0);
			if (!planning::keyboardflag['3'])
			{
				auto& arcs = graphSearch::VRINs[sliceIdx].arcs;
				auto& color = graphSearch::VRINs[sliceIdx].color;
				for (int i = arcs.size() - 1; i >= 0; i--)
				{
					if (color[i] == 1.0)
						arcs[i].draw();
					else
						break;
				}
				//for (auto& arc : boundary)
				//	arc.draw();
			}

			// 4. draw clicked point
			static bool clicked = false;
			if (!clicked)
			{
				if (ms::clickedPoint.x() != 0 || ms::clickedPoint.y() != 0)
					clicked = true;
			}
			else if(!planning::keyboardflag['0'])
			{
				
				Point cp;
				bool t0, t1;
				t0 = graphSearch::testers2[sliceIdx].test(ms::clickedPoint);
				t1 = graphSearch::testers2[sliceIdx].testPrecise(ms::clickedPoint, cp);

				// result of t0
				glPointSize(5.5);
				if (t0)
					glColor3f(1, 0, 0);
				else
					glColor3f(0, 0, 1);
				glBegin(GL_POINTS);
				glVertex2dv(ms::clickedPoint.P);
				glEnd();
				glPointSize(1.0);


				// result of t1
				if(t1)
					glColor3f(1, 0, 0);
				else
					glColor3f(0, 0, 1);
				glBegin(GL_LINES);
				glVertex2dv(ms::clickedPoint.P);
				glVertex2dv(cp.P);
				glEnd();

				double rad = sqrt((ms::clickedPoint - cp).length2());
				CircularArc arc = cd::constructArc(ms::clickedPoint, rad, 0, PI2 - 1e-4);
				arc.draw2();

			}
			else
			{
				//draw maxtouch
				if (clickedPoints.size() % 3 == 2)
				{
					cout << "mt calc query" << endl;

					auto it = clickedPoints.end();
					it--;
					auto p1 = *it;
					it--;
					auto p0 = *it;
					auto n0 = (p1 - p0).normalize();
					double par, bestPar;
					CircularArc bestArc;

					double rmin = 1e8;

					
					for (auto& arc : graphSearch::VRINs[sliceIdx].arcs)
					{
						auto rad = arc.maxTouchRad(p0, n0, par);
						if (rad > 1e-300)
						{
							if (rad < rmin)
							{
								rmin = rad;
								bestPar = par;
								bestArc = arc;
							}
						}
					}

					if (rmin != 1e8)
					{
						auto c = p0 + rmin * n0;
						CircularArc ca;
						ca = cd::constructArc(c, rmin, 1e-10, PI2 - 1e-10);
						ca.draw2();

						glPointSize(5.0);
						glBegin(GL_POINT);
						glVertex2dv(p0.P);
						glEnd();
						glPointSize(1.0);


						glBegin(GL_LINES);
						glVertex2dv(p0.P);
						glVertex2dv(ca.cc().P);
						glEnd();

						cout << "minRad : " << rmin << " Par : " << bestPar << endl;
						cout << bestArc << endl;
					}

				}

			}

			// 5. draw disk
			if (planning::keyboardflag['4'])
			{
				glColor3f(0, 0, 0);
				auto& disks = graphSearch::testers2[sliceIdx].getConvDisk();
				for (auto& disk : disks)
					disk.draw();
			}

			// 6. draw offset
			if (planning::keyboardflag['5'])
			{
				static int offSetMax = 0;
				char temp;
				temp = '[';
				if (planning::keyboardflag[temp] != planning::keyboardflag_last[temp])
					offSetMax--;
				temp = ']';
				if (planning::keyboardflag[temp] != planning::keyboardflag_last[temp])
					offSetMax++;

				auto& offCrv = graphSearch::offCrv;
				glLineWidth(1.0f);
				glBegin(GL_LINES);
				for (int i = 0; i < offSetMax && i < offCrv.size(); i++)
				{
					// color
					double t = double(i) / offCrv.size();
					float
						r = t * 1.0,
						g = (1 - t) * 1.0,
						b = 0.0;
					glColor3f(r, g, b);

					// draw
					auto& offs = offCrv[i][sliceIdx];
					for (auto& l : offs)
					{
						glVertex2d(l.v0.x(), l.v0.y());
						glVertex2d(l.v1.x(), l.v1.y());
					}
				}
				glEnd();
			}

			// 7. test support func
			if (planning::keyboardflag['6'])
			{
				//dbg
				//sliceIdx = 296;

				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				glColor3f(0, 0, 0);

				static int hlt = -1;
				char temp;
				temp = '[';
				if (planning::keyboardflag[temp] != planning::keyboardflag_last[temp])
					hlt--;
				temp = ']';
				if (planning::keyboardflag[temp] != planning::keyboardflag_last[temp])
					hlt++;

				static int mc = 1000;
				temp = '-';
				if (planning::keyboardflag[temp] != planning::keyboardflag_last[temp])
					mc--;
				temp = '=';
				if (planning::keyboardflag[temp] != planning::keyboardflag_last[temp])
					mc++;


				vector<SupportComplex> vsc;
				int c = 0;
				for (auto& loop : MR)
				{
					for (auto& as : loop)
					{
						for (auto& arc : as.Arcs)
						{
							Support s(arc);
							vsc.emplace_back(s);

							c++;
							if (c == mc)
								break;
						}
						if (c == mc)
							break;
					}
					if (c == mc)
						break;
				}

				if (Key('p'))
				{
					for (auto& s : vsc)
						s.flip();
				}

				////debug
				//{
				//	CircularArc a0, a1;
				//	a0 = cd::constructArc(Point(1, 0), Point(0, 1), Point(0, 1));
				//	a1 = cd::constructArc(Point(0, -1), Point(-1, 0), Point(0, 1));
				//
				//	vsc.resize(0);
				//	vsc.emplace_back(Support(a0));
				//	vsc.emplace_back(Support(a1));
				//
				//	Support(a0).print();
				//	Support(a1).print();
				//
				//	vsc[0].print();
				//	vsc[1].print();
				//}

				SupportComplex env = vsc.front();
				for (int i = 1; i < vsc.size(); i++)
				{

					////dbg_out derr
					//{
					//	if (sliceIdx == 296 && i == 7)
					//	{
					//		derr << "!!!!!!!!!env " << endl;
					//		env.print();
					//		derr << "!!!!!!!!!vsc[i] " << endl;
					//		vsc[i].print();
					//		derr << "!!!!!!!!!new env" << endl;
					//		debugBlock.push_back(1.0);
					//		env.upperEnvelope(vsc[i]).print();
					//		debugBlock.resize(0);
					//	}
					//}

					env = env.upperEnvelope(vsc[i]);
				}


				// do drawing
				for (int i = 0; i< vsc.size(); i++)
				{
					if (hlt == i)
						glLineWidth(5.0);
					else
						glLineWidth(1.0);
					vsc[i].draw();
				}

				if (hlt == -1)
					glLineWidth(10.0);
				else
					glLineWidth(1.0);
				glColor3f(0, 1, 0);
				env.draw();

				cout << "highlight idx : " << hlt << endl;


				goto end;
			}

			// 8. test COnvex Hull
			if (planning::keyboardflag['7'])
			{
				vector<CircularArc> temp;

				size_t length = MR.size();
				for (size_t i = 0; i < length; i++)
				{
					if(MIB[i])
					for (auto& as : MR[i])
						for (auto& a : as.Arcs)
							temp.push_back(a);
				}

				CvxHullCalc cc;
				cc.setInput(temp);

				cc.setOutput(ch);

				cc.calcAlg();

				for (auto& a : ch)
				{
					a.draw();
				}

				// drawOffset too
				if (Key('8'))
				{
					
					static double offset = 0.3;
					if (Key3('-'))
						offset -= 0.01;
					if (Key3('='))
						offset += 0.01;
					if (offset < 0.01)
						offset = 0.01;
					CvxHullCalc::offset(ch, chOff, offset);
					cout << "cvxH size : " << ch.size() << endl;
					cout << "CvxHull offset : " << offset << "    , size : " << chOff.size() << endl;
					for (auto& a : chOff)
					{
						a.draw();
					}

					if (Key3('p'))
					{
						cout << "!!!!!!!!!!!!! CVX HULL" << endl;
						int i = 0;
						for (auto a : ch)
						{
							cout << i << endl;
							cout << a << endl<<endl;
							i++;
						}

						cout << "!!!!!!!!!!!!!!!!!!! OFF CVX HULL " << endl;
						for (auto a : chOff)
							cout << a << endl<<endl;

						cout << "!!!! NON G1 " << endl;
						auto& in = ch;
						for (int i = 0; i < in.size(); i++)
						{
							// 1. find next arc
							int iNext = i + 1;
							if (iNext == in.size())
								iNext == 0;

							// 3. check g1
							const double errG1 = 1e-8;
							bool g1;
							{
								auto dot = in[i].n1() * in[iNext].n0();
								if (fabs(dot) > 1 - errG1) // note : theoeretically, fabs is not needed.
									g1 = true;
								else
									g1 = false;
							}

							// 4. if(not g1, add intermediate arc)
							if (!g1)
							{
								cout << "non g1 between : " << i << " , " << iNext << endl;
							}

						}
					}
				}

				
				//
			}

			// 9. draw common tangent
			if (Key('9'))
			{
				for(int i = 0; i < MR.size(); i++)
					for (int j = i + 1; j < MR.size(); j++)
					{
						if (MIB[i] && MIB[j])
						{
							auto& l0 = MR[i];
							auto& l1 = MR[j];

							vector<CircularArc>
								loop0, loop1;

							for (auto& as : l0)
								for (auto& a : as.Arcs)
									loop0.push_back(a);

							for (auto& as : l1)
								for (auto& a : as.Arcs)
									loop1.push_back(a);

							ComTanCalc cal;
							vector<Point> pt;
							cal.setInput0(loop0);
							cal.setInput1(loop1);
							cal.setOutput(pt);
							cal.calc();

							glBegin(GL_LINES);
							for (auto& p : pt)
								glVertex2dv(p.P);
							glEnd();

						}
					}
			}

			end:
			// 99. input
			for (int i = 0; i < 256; i++)
				planning::keyboardflag_last[i] = keyboardflag[i];
			glutSwapBuffers();
		};

		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
		glutInitWindowSize(wd, ht);
		glutCreateWindow("Data Check");
		glutReshapeFunc(reshapeFunc);
		glutDisplayFunc(displayFunc);
		glutMouseFunc(mouse_callback);
		glutKeyboardFunc(keyboard_callback);
		glutIdleFunc(idleFunc);

		glEnable(GL_DEPTH_TEST);

		glutMainLoop();
	}

	namespace renderPathGlobal
	{
		// path info
		double* pathPtr;
		int pathSize;
		int pathIdx = 0;
		int pathType = 0;

		vector<Vertex> pathGraphSearch;
		vector<Vertex> mcPath;

		// I/O
		int clickStatus = 0;
		bool clickNow = false;
		double
			sx, sy, sz = 0,
			ex, ey, ez = 0;
		double
			clickedZ = 0.0;

		// Astar
		vector<Vertex> vecVertices;
		map<Vertex, int, VertexLessFn> mapLookup;
		Graph theGr;
		Vertex ver0, ver1;

		// Biarc approx
		vector<CircularArc> arcPath;
		vector<double> arcPathRot;
		vector<Vertex> arcPathTess;

		/*
		Def:
		
		Note:
			Does not care loop feature of theta
		*/
		void refinePath(vector<Vertex>& in, vector<Vertex>& out)
		{
			auto& lTest = graphSearch::testers;
			auto& pTest = graphSearch::testers2;

			// Naive 3d line tester
			// Def : test 3d line
			// return true iff collision
			auto lTest3D = [&pTest, &lTest](Vertex& p0, Vertex& p1)
			{
				int sliceIdx0 = p0.z / 360.0 * ms::numofframe;
				int sliceIdx1 = p1.z / 360.0 * ms::numofframe;

				// simple case
				if (sliceIdx0 == sliceIdx1)
				{
					return lTest[sliceIdx0].test(Point(p0.x, p0.y), Point(p1.x, p1.y));
				}

				double
					x0 = p0.x,
					y0 = p0.y,
					z0 = p0.z,
					dx = p1.x - p0.x,
					dy = p1.y - p0.y,
					dz = p1.z - p0.z;
				dx /= dz;
				dy /= dz;
				dz = 1.0;

				Point v0(p0.x, p0.y);
				Point v1;
				v1.x() = v0.x() + dx;
				v1.y() = v0.y() + dy;

				// loop
				{
					//ready for loop
					int increment = 1;
					if (sliceIdx0 > sliceIdx1)
						increment = -1;

					// just sample points and do tests
					for (int i = sliceIdx0; i != sliceIdx1; i += increment)
					{
						if (lTest[i].test(v0, v1))return true;
						if (lTest[i + increment].test(v0, v1))return true;


						// v0, v1
						v0 = v1;
						v1.x() = v0.x() + dx;
						v1.y() = v0.y() + dy;
					}
					
				}
				return false;
			};

			// 1. big refine first (rverses order)
			vector<Vertex> temp2;
			{
				const int loopLimit = 7;

				int sz = in.size();
				auto temp = in;

				while (temp.size() > 1)
				{
					int eIdx = temp.size() - 1;
					int bIdx = 0;
					int sampleLength = temp.size() - 1;

					int sizeBefore = temp.size();

					for (int i = 0; i < loopLimit; i++)
					{
						bIdx = eIdx - sampleLength;
						
						// err check
						if (sampleLength <= 1)
						{
							temp2.push_back(temp[eIdx]);
							temp.pop_back();
							break;
						}

						if (bIdx < 0)
						{
							sampleLength /= 2;
							continue;
						}

						// test
						int testRes = lTest3D(temp[bIdx], temp[eIdx]);
						if (testRes)
						{
							//when collision
							if (i == loopLimit - 1)
							{
								//// last loop
								//for (int j = eIdx; j > bIdx; j--)
								//{
								//	temp2.push_back(temp[j]);
								//}
								//temp.resize(bIdx + 1);
								//break;

								temp2.push_back(temp[eIdx]);
								temp.pop_back();
								break;
							}
							else
							{
								sampleLength /= 2;
								continue;
							}

						}
						else
						{
							//no collision -> reduce to a line;
							temp2.push_back(temp[eIdx]);
							temp.resize(bIdx + 1);
							break;

						}
					}

					int sizeAfter = temp.size();

					if (sizeAfter >= sizeBefore)
					{
						cerr << "!!!ERR : in refinePath : size of temp did not decrease" << endl;
						break;
					}
				}

				for (int i = temp.size() - 1; i > -1; i--)
				{
					temp2.push_back(temp[i]);
				}
			}


			// 99. tesselate;
			reverse(temp2.begin(), temp2.end());
			//dbg
			{
				out = temp2;
				return;
			}

			auto& beforTess = temp2;
			auto& afterTess = out;
			for (int i = 0; i < temp2.size() - 1; i++)
			{
				const double dl = 0.01;

				auto& v0 = temp2[i];
				auto& v1 = temp2[i+1];
				auto dv = v0 - v1;
				auto l = dv.norm();
				int n = l / dl;
				if (n < 2)
				{
					afterTess.push_back(temp2[i]);
				}
				else
				{
					for (int j = 0; j < n; j++)
					{
						double t = double(j) / n;
						double t1 = 1.0 - t;
						auto v = v0 * t + v1 * t1;
						afterTess.push_back(v);
					}
				}
			}


		}
	}

	void renderPath(
		int argc, 
		char* argv[],
		vector<double>& path
		)
	{
		using namespace renderPathGlobal;

		zoom = 1.0;
		wd = ht = 600;

		pathPtr = &(path[0]);
		pathSize = (path.size() / 3);

		auto reshapeFunc = [](GLint w, GLint h) ->void
		{
			double ratioWH = double(wd) / ht;
			if (w < h * ratioWH) {
				wd = w;
				ht = w / ratioWH;
			}
			else {
				wd = h * ratioWH;
				ht = h;
			}
		};
		auto    idleFunc = []() -> void
		{
			//manage speed
			int millisecondPerFrame = 100;
			static clock_t lastTime = clock();
			while (clock() - lastTime < millisecondPerFrame)
			{

			}
			lastTime = clock();

			// change pathIdx
			if (!planning::keyboardflag['f'])
				pathIdx++;
			static bool lastG = false, lastH = false;
			if (planning::keyboardflag['g'] != lastG)
				pathIdx++;
			if (planning::keyboardflag['h'] != lastH)
				pathIdx--;
			lastG = planning::keyboardflag['g'];
			lastH = planning::keyboardflag['h'];


			if (pathIdx >= pathSize) pathIdx = 0;
			if (pathIdx < 0) pathIdx = pathSize - 1;

			cout << "Current PathIdx : " << pathIdx << endl;
			if (pathIdx < mcPath.size() && pathIdx >= 0)
				cout 
				<< "   " << mcPath[pathIdx].x 
				<< "   " << mcPath[pathIdx].y 
				<< "   " << mcPath[pathIdx].z << endl;

			glutPostRedisplay();
		};
		auto displayFunc = []()
		{
			glClearColor(1.0, 1.0, 1.0, 1.0);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			setup_viewvolume();
			setup_transform();

			// Rendering code....
			glLineWidth(2.0f);
			glPointSize(2.8f);
			glColor3f(0.0f, 0.0f, 0.0f);
			glViewport(0, 0, wd, ht);

			auto sceneIdx = 7;
			auto robotIdx = 1;

			// 1. draw scene
			auto& scene = Models_Approx[sceneIdx];
			for (auto& as : scene)
				as.draw();

			//// 2. draw robot
			//auto& robot = Models_Approx[robotIdx];
			//Point translation(pathPtr[3 * pathIdx + 0], pathPtr[3 * pathIdx + 1]);
			//double rotationDegree = pathPtr[3 * pathIdx + 2];
			//for (auto& as : robot)
			//	for (auto& arc : as.Arcs)
			//	{
			//		cd::translateArc(cd::rotateArc(arc, rotationDegree + 180.0), translation).draw();
			//	}
			//glPointSize(4.0f);
			//glBegin(GL_POINTS);
			//glVertex3d(translation.P[0], translation.P[1], -0.5);
			//glEnd();

			if (planning::keyboardflag['1'])
			{
				//// 3. draw path
				//glBegin(GL_LINE_STRIP);
				//for (size_t i = 0; i < pathSize; i++)
				//{
				//	double ratio = (pathPtr[3 * i + 2]) / 360.0;
				//	glColor3f(1 - ratio, ratio, 0);
				//	glVertex2d(pathPtr[3 * i + 0], pathPtr[3 * i + 1]);
				//}
				//glEnd();
				//
				//// 4.draw more robots;
				//glColor3f(0.8, 0.8, 0.8);
				//int nsample = pathSize / 50;
				//if (nsample > 0)
				//{
				//	int d = pathSize / nsample;
				//	for (int i = 0; i < nsample; i++)
				//	{
				//		int idx = d * i;
				//
				//		auto& robot = Models_Approx[robotIdx];
				//		Point translation(pathPtr[3 * idx + 0], pathPtr[3 * idx + 1]);
				//		double rotationDegree = pathPtr[3 * idx + 2];
				//		for (auto& as : robot)
				//			for (auto& arc : as.Arcs)
				//			{
				//				cd::translateArc(cd::rotateArc(arc, rotationDegree + 180.0), translation).draw();
				//			}
				//		glPointSize(4.0f);
				//		glBegin(GL_POINTS);
				//		glVertex3d(translation.P[0], translation.P[1], -0.5);
				//		glEnd();
				//	}
				//
				//	int idx = pathSize - 1;
				//	{
				//		auto& robot = Models_Approx[robotIdx];
				//		Point translation(pathPtr[3 * idx + 0], pathPtr[3 * idx + 1]);
				//		double rotationDegree = pathPtr[3 * idx + 2];
				//		for (auto& as : robot)
				//			for (auto& arc : as.Arcs)
				//			{
				//				cd::translateArc(cd::rotateArc(arc, rotationDegree + 180.0), translation).draw();
				//			}
				//		glPointSize(4.0f);
				//		glBegin(GL_POINTS);
				//		glVertex3d(translation.P[0], translation.P[1], -0.5);
				//		glEnd();
				//	}
				//}

				// draw CSO
				if(pathIdx < mcPath.size())
				{
					int sliceNo = mcPath[pathIdx].z / 360.0 * ms::numofframe;
					while (sliceNo >= ms::numofframe)
						sliceNo -= ms::numofframe;
					while (sliceNo < 0)
						sliceNo += ms::numofframe;

					//cout << "Sno: " <<  sliceNo << endl;
					glColor3f(0.7, 0.7, 0.7);
					for (auto& arc : graphSearch::VRINs[sliceNo].arcs)
						arc.draw2();
				}
				
				//debugging test
				if (false)
				{
					vector<CircularArc> temp;
					vector<Circle> res;

					int sliceNo = mcPath[pathIdx].z / 360.0 * ms::numofframe;
					while (sliceNo >= ms::numofframe)
						sliceNo -= ms::numofframe;
					while (sliceNo < 0)
						sliceNo += ms::numofframe;

					// build temp
					for (int i = 0; i < graphSearch::VRINs[sliceNo].arcs.size(); i++)
					{
						if (graphSearch::VRINs[sliceNo].color[i] == 1.0)
							break;

						
						temp.push_back(flipArc(graphSearch::VRINs[sliceNo].arcs[i]));
					}

					// build res
					findInteriorDisks(temp, 1, res);

					for (auto& c : res)
						c.draw();
				}
			}



			// 5. draw robot when robot is being pressed;
			if (clickNow)
				if (clickStatus == 0 || clickStatus == 1)
					{
						if (clickStatus == 0)
						{
							glColor3f(1, 0, 0);
						}
						else
						{
							glColor3f(0, 0, 1);
						}

						auto& robot = Models_Approx[robotIdx];
						Point translation = clickedPoint;
						double rotationDegree = clickedZ;
						for (auto& as : robot)
							for (auto& arc : as.Arcs)
							{
								cd::translateArc(cd::rotateArc(arc, rotationDegree + 180.0), translation).draw();
							}
						glPointSize(4.0f);
						glBegin(GL_POINTS);
						glVertex3d(translation.P[0], translation.P[1], -0.5);
						glEnd();
					}

			// 6. draw mcPath
			if (planning::keyboardflag['2'])
			{
				glColor3f(0.3, 0.3, 1.0);
				{
					glBegin(GL_LINE_STRIP);
					for (auto& v : mcPath)
					{
						glVertex2d(v.x, v.y);
					}
					glEnd();
				}
			}

			// 7. draw optimized path  //changed from://draw arcPath
			glColor3f(0, 0, 0);
			if (planning::keyboardflag['3'])
				//{
				//	for (auto& arc : arcPath)
				//		arc.draw2();
				//}
			{
				glColor3f(0, 0.8, 0);
				glBegin(GL_LINE_STRIP);
				for (auto& v : pathGraphSearch)
					glVertex2d(v.x, v.y);
				glEnd();
			}

			// 8. robot following path
			glColor3f(0, 0, 0);
			if(pathType == 1)
			if(pathIdx < arcPathTess.size())
			{
				auto& robot = Models_Approx[robotIdx];
				Point translation(arcPathTess[pathIdx].x, arcPathTess[pathIdx].y);
				double rotationDegree = arcPathTess[pathIdx].z;
				for (auto& as : robot)
					for (auto& arc : as.Arcs)
					{
						cd::translateArc(cd::rotateArc(arc, rotationDegree + 180.0), translation).draw();
					}
				glPointSize(4.0f);
				glBegin(GL_POINTS);
				glVertex3d(translation.P[0], translation.P[1], -0.5);
				glEnd();
			}
			if (pathType == 0 && (pathIdx >= 0 && pathIdx < mcPath.size()))
			{
				auto& robot = Models_Approx[robotIdx];
				Point translation(mcPath[pathIdx].x, mcPath[pathIdx].y);
				double rotationDegree = mcPath[pathIdx].z;
				for (auto& as : robot)
					for (auto& arc : as.Arcs)
					{
						cd::translateArc(cd::rotateArc(arc, rotationDegree + 180.0), translation).draw();
					}
				glPointSize(4.0f);
				glBegin(GL_POINTS);
				glVertex3d(translation.P[0], translation.P[1], -0.5);
				glEnd();

			}

			// 9. draw starting/end point
			if (!planning::keyboardflag['4'] && mcPath.size() > 0)
			{
				// start
				{
					glColor3f(1, 0.8, 0.8);
					auto& robot = Models_Approx[robotIdx];
					Point translation(mcPath.front().x, mcPath.front().y);
					double rotationDegree = mcPath.front().z;
					for (auto& as : robot)
						for (auto& arc : as.Arcs)
						{
							cd::translateArc(cd::rotateArc(arc, rotationDegree + 180.0), translation).draw();
						}
					glPointSize(4.0f);
					glBegin(GL_POINTS);
					glVertex3d(translation.P[0], translation.P[1], -0.5);
					glEnd();
				}

				// end
				{
					glColor3f(0.8, 0.8, 1);
					auto& robot = Models_Approx[robotIdx];
					Point translation(mcPath.back().x, mcPath.back().y);
					double rotationDegree = mcPath.back().z;
					for (auto& as : robot)
						for (auto& arc : as.Arcs)
						{
							cd::translateArc(cd::rotateArc(arc, rotationDegree + 180.0), translation).draw();
						}
					glPointSize(4.0f);
					glBegin(GL_POINTS);
					glVertex3d(translation.P[0], translation.P[1], -0.5);
					glEnd();
				}
			}

			if (Key('5'))
			{
				using namespace gs;
				// draw all edges in graph
				auto& e = dijk.getEdge();
				auto& v = dijk.getVert();

				glBegin(GL_LINES);
				for (auto& e0 : e)
				{
					auto& v0 = v[e0.first];
					auto& v1 = v[e0.second];
					glVertex2d(v0.x(), v0.y());
					glVertex2d(v1.x(), v1.y());
				}
				glEnd();
			}

			// if(6) draw single-slice-edge
			if (Key('6'))
			{
				using namespace gs;
				for (int sn = 0; sn < dijk._vMat.size(); sn++)
				{
					auto& v = dijk._vMat[sn];
					auto& e = dijk._eMat[sn];

					glBegin(GL_LINES);
					for (auto& e0 : e)
					{
						auto& v0 = v[e0.first];
						auto& v1 = v[e0.second];
						glVertex2d(v0.x(), v0.y());
						glVertex2d(v1.x(), v1.y());
					}
					glEnd();
				}
			}

			// if(7) draw interslice edge
			if (Key('7'))
			{
				using namespace gs;
				for (int sn = 0; sn < dijk._vMat.size(); sn++)
				{
					auto snn = sn + 1;
					if (snn == dijk._vMat.size())
						snn = 0;

					auto& v = dijk._vMat[sn];
					auto& v2 = dijk._vMat[snn];
					auto& e = dijk._eMat2[sn];

					//double max = -1.0;

					int i = 0;
					glBegin(GL_LINES);
					for (auto& e0 : e)
					{
						auto& v0 = v[e0.first];
						auto& v1 = v2[e0.second];
						glVertex2d(v0.x(), v0.y());
						glVertex2d(v1.x(), v1.y());

						////dbg_out
						//{
						//	//auto d = (v0 - v1).length1();
						//	auto d = dijk._wMat[sn][i];
						//	if (d > max)
						//		max = d;
						//}

						//if (i == 0)
						//{
						//	cout << COUT(sn) << endl;
						//	cout << COUT(e0.first) << endl;
						//	cout << COUT(e0.second) << endl;
						//	cout << COUT(v0) << endl;
						//	cout << COUT(v1) << endl;
						//}
						//i++;
					}
					glEnd();
					//cout << " longest interslice length: " << max << endl;

				}
			}

			// 98. change mode
			if (planning::keyboardflag['0'] != planning::keyboardflag_last['0'])
			{
				const int nPathType = 2;
				pathType = (pathType + 1) % nPathType;
				pathIdx = 0;

				switch (pathType)
				{
				case 0:
					pathSize = mcPath.size();
					cout << "Biarc Path " << endl;
					break;
				case 1:
					pathSize = arcPathTess.size();
					cout << "MC Path" << endl;
					break;
				}
			}

			// 99. camera
			{
				char temp;
				
				temp = '[';
				if (planning::keyboardflag[temp] != planning::keyboardflag_last[temp])
					zoom *= 0.95;
				temp = ']';
				if (planning::keyboardflag[temp] != planning::keyboardflag_last[temp])
					zoom /= 0.95;

			}

			for (int i = 0; i < 256; i++)
				planning::keyboardflag_last[i] = keyboardflag[i];

			glutSwapBuffers();
			};
		auto mouseFunc = [](int button, int action, int x, int y) ->void
		{

			if (button == 3)
				clickedZ += 15.0;

			if (button == 4)
				clickedZ -= 15.0;

			if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN)
			{
				clickNow = true;
			}

			if (button == GLUT_LEFT_BUTTON && action == GLUT_UP)
			{
				clickNow = false;

				// 1. get clicked point
				auto fx = float(x - wd / 2) / (wd / 2);
				auto fy = float(-y + ht / 2) / (ht / 2);
				fx = zoom * fx + tx;
				fy = zoom * fy + ty;
				clickedPoint.x() = fx;
				clickedPoint.y() = fy;

				// 2.
				switch (clickStatus)
				{
				case 0:
					// set start point
					sx = clickedPoint.x();
					sy = clickedPoint.y();
					sz = clickedZ;
					clickStatus = 1;
					cout << "Start point : " << sx << ", " << sy << ", " << sz << endl;
					break;

				case 1:
					// set end point
					ex = clickedPoint.x();
					ey = clickedPoint.y();
					ez = clickedZ;
					clickStatus = 2;
					cout << "End point : " << ex << ", " << ey << ", " << ez << endl;
					break;

				case 2:
					{
						// Alias
						std::vector<std::vector<v_edge>>&
							v_edges = planning::output_to_file::v_edges;
						int szIdx = sz / 360.0 * ms::numofframe; while(szIdx < 0) szIdx += ms::numofframe; szIdx %= ms::numofframe;
						int ezIdx = ez / 360.0 * ms::numofframe; while(ezIdx < 0) ezIdx += ms::numofframe; ezIdx %= ms::numofframe;
						{
							//just to make 64-bits of z value all same with that of a-star graph
							double dz = 360.0 / v_edges.size();

							sz = 0.0;
							ez = 0.0;
							for (int i = 0; i < szIdx; i++)
								sz += dz;
							for (int i = 0; i < ezIdx; i++)
								ez += dz;
						}

						// i. find closest v_edge

						// i-1. start
						double closestDist = 1E8;
						Point closestPoint;
						Point temp(sx, sy);

						// for (All ve)
						for (int i = 0; i < v_edges[szIdx].size(); i++)
						{
							auto& ve = v_edges[szIdx][i];
							auto& p = ve.v0;
							auto  d = sqrt((temp - p).length2());
							if (d < closestDist && d < ve.clr0())
							{
								closestDist = d;
								closestPoint = p;
							}
						}
						if (closestDist == 1E8)
						{
							cerr << "Invalid Input to Astar" << endl;
							clickStatus = 0;
							break;
						}
						ver0 = Vertex(closestPoint.x(), closestPoint.y(), sz);

						// i-2. end point
						closestDist = 1E8;
						temp = Point(ex, ey);
						// for (All ve)
						for (int i = 0; i < v_edges[ezIdx].size(); i++)
						{
							auto& ve = v_edges[ezIdx][i];
							auto& p = ve.v0;
							auto  d = sqrt((temp - p).length2());
							if (d < closestDist && d < ve.clr0())
							{
								closestDist = d;
								closestPoint = p;
							}
						}
						if (closestDist == 1E8)
						{
							cerr << "Invalid Input to Astar" << endl;
							clickStatus = 0;
							break;
						}
						ver1 = Vertex(closestPoint.x(), closestPoint.y(), ez);
						
						
						// ii. run a-star
						std::vector<Vertex> path;
						vector<xyt> dijkresult;
						std::vector<double> vertClr;
						std::vector<gs::EdgeType> ets;
						bool useDijk = true;
						bool useDijk2 = true;
						if (useDijk)
						{
							if (useDijk2)
							{
								xyt
									s(sx, sy, sz),
									e(ex, ey, ez);

								if (gs::dijk.searchGraph2(s, e, dijkresult, vertClr, ets, robotForward))
								{
									for (auto v : dijkresult)
									{
										Vertex temp;
										temp.x = v.x();
										temp.y = v.y();
										temp.z = v.t();

										path.push_back(temp);
									}
								}
								else
								{
									cout << "Dijkstra: no path" << endl;
									auto t0 = std::chrono::high_resolution_clock::now();
									auto t1 = t0;
									auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
									while (ms < 1000.0)
									{
										t1 = std::chrono::high_resolution_clock::now();
										ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
									}
								}
							}
							else
							{
								vector<int> vertIdxList;
								xyt
									s(sx, sy, sz),
									e(ex, ey, ez);

								// if(path exists)
								if (gs::dijk.searchGraph(s, e, vertIdxList))
								{
									Vertex sv(sx, sy, sz);
									Vertex ev(ex, ey, ez);

									path.push_back(sv);
									for (auto i : vertIdxList)
									{
										auto& v = gs::dijk.getVert()[i];

										Vertex temp;
										temp.x = v.x();
										temp.y = v.y();
										temp.z = v.t();
										path.push_back(temp);
									}
									path.push_back(ev);

									// debug out
									int i = 0;
									for (auto& p : path)
									{
										cout << "dijk path: " << i << ":" << p.x << " , " << p.y << " , " << p.z << endl;
										i++;
									}

								}
								else
								{
									cout << "Dijkstra: no path" << endl;
									auto t0 = std::chrono::high_resolution_clock::now();
									auto t1 = t0;
									auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
									while (ms < 1000.0)
									{
										t1 = std::chrono::high_resolution_clock::now();
										ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
									}
								}
							}
						}
						else
							path = invoke_AStar(theGr, vecVertices, mapLookup, ver0, ver1);
						
						if(false)
						{
							// path : approx biarc
							arcPath.clear();
							arcPathRot.clear();
							refinePathBiarc(path, 10, arcPath, arcPathRot);
							//dbg_out
							{
								for (auto a : arcPath)
									cout << "AP : " << a.x0() << "   " << a.x1() << endl;

								for (auto a : arcPath)
									cout << "AP norm : " << a.n0() << "   " << a.n1() << endl;
							}

							// path : tess
							arcPathTess.clear();
							tessPathBiarc(arcPath, arcPathRot, 0.004, arcPathTess);
							pathSize = arcPathTess.size();
							pathIdx = 0;

							// path : 
							//tessPathAlignTheta(arcPathTess, graphSearch::testers2, Point(0, 1));
						}
						if (true)
						{
							for (int i = 0; i < ets.size(); i++)
							{
								cout 
									<< "edge " << i << " : " << gs::et2str(ets[i]) << endl
									<< "  vert0 : " << path[i  ].x << ", " << path[i  ].y << ", " << path[i  ].z << ", " << vertClr[i  ] << endl
									<< "  vert1 : " << path[i+1].x << ", " << path[i+1].y << ", " << path[i+1].z << ", " << vertClr[i+1] << endl;
							}
						}

						// save mcPath
						pathGraphSearch = path;

						vector<xyt> optPath;
						optimizePath(dijkresult, vertClr, ets, optPath);

						vector<Vertex> path2;
						for (auto v : optPath)
						{
							Vertex temp;
							temp.x = v.x();
							temp.y = v.y();
							temp.z = v.t();

							path2.push_back(temp);
						}

						mcPath.clear(); //drawing path
						tessPathClear(path2, 0.004, mcPath);
						pathSize = mcPath.size();
						pathType = 0;
						pathIdx = 0;

						{
						//// iii. do refinement.
						//vector<Vertex> ref;
						//refinePath(path, ref);
						//
						//bool refiningDone = true;
						//if(!refiningDone)
						//{
						//	//temp
						//	if (mcPath.size() == 0)
						//	{
						//
						//	}
						//	else
						//	{
						//		delete[] pathPtr;
						//	}
						//
						//	pathSize = path.size() + 2;
						//	pathPtr = new double[(path.size()+2) * 3];
						//	pathIdx = 0;
						//
						//	for (int i = 0; i < path.size(); i++)
						//	{
						//		pathPtr[3 * i + 3] = path[i].x;
						//		pathPtr[3 * i + 4] = path[i].y;
						//		pathPtr[3 * i + 5] = path[i].z;
						//	}
						//
						//	pathPtr[0] = sx;
						//	pathPtr[1] = sy;
						//	pathPtr[2] = sz;
						//	pathPtr[(path.size() + 1) * 3 + 0] = ex;
						//	pathPtr[(path.size() + 1) * 3 + 1] = ey;
						//	pathPtr[(path.size() + 1) * 3 + 2] = ez;
						//
						//
						//}
						//else
						//{
						//	//temp
						//	if (mcPath.size() == 0)
						//	{
						//
						//	}
						//	else
						//	{
						//		delete[] pathPtr;
						//	}
						//
						//	pathSize = ref.size();
						//	pathPtr = new double[(pathSize) * 3];
						//	pathIdx = 0;
						//
						//	for (int i = 0; i < pathSize; i++)
						//	{
						//		pathPtr[3 * i + 0] = ref[i].x;
						//		pathPtr[3 * i + 1] = ref[i].y;
						//		pathPtr[3 * i + 2] = ref[i].z;
						//	}
						//
						//	mcPath = path;
						//}
						}

						clickStatus = 0;
						break;
					}
				default:
					break;
				
				}

			}

			if (button == GLUT_RIGHT_BUTTON)
			{
				clickStatus = 0;
			}

			glutPostRedisplay();
		};
		auto motionFunc = [](int x, int y)
		{
			auto fx = float(x - wd / 2) / (wd / 2);
			auto fy = float(-y + ht / 2) / (ht / 2);
			fx = zoom * fx + tx;
			fy = zoom * fy + ty;
			clickedPoint.x() = fx;
			clickedPoint.y() = fy;
		};

		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
		glutInitWindowSize(wd, ht);
		glutCreateWindow("2D planning");
		glutReshapeFunc(reshapeFunc);
		glutDisplayFunc(displayFunc);
		glutMouseFunc(mouseFunc);
		glutKeyboardFunc(keyboard_callback);
		glutMotionFunc(motionFunc);
		glutIdleFunc(idleFunc);

		glEnable(GL_DEPTH_TEST);

		glutMainLoop();
	}


	namespace renderRefinementCollisionTestGlobal
	{
		// as GLUT needs global variables...
		// not the best implementaion, but these part will only be used for testing/error-checking

		std::vector<decltype(ms::Model_Result)> MRs;
		std::vector<decltype(ms::ModelInfo_Boundary)> MIBs;
		std::vector<std::vector<planning::output_to_file::v_edge>> v_edges;
		decltype(planning::voronoiBoundary) voronoiBoundary;
		std::vector<planning::VR_IN> VRINs;
		int& sliceIdx = ms::t2;

		std::vector<cd::lineSegmentCollisionTester> testers;
	}
	/*
	almost a copy of renderMinkVoronoi,
	just that it renders overlayed collision test results.

	*/
	void renderRefinementCollisionTest(
		int argc,
		char* argv[],
		std::vector<decltype(ms::Model_Result)>& MRs,
		std::vector<decltype(ms::ModelInfo_Boundary)>& MIBs,
		std::vector<std::vector<planning::output_to_file::v_edge>>& v_edges,
		decltype(planning::voronoiBoundary)& voronoiBoundary,
		std::vector<planning::VR_IN> VRINs
		)
	{
		wd = 1200, ht = 800;

		//inefficiency... copy data
		renderRefinementCollisionTestGlobal::MRs = MRs;
		renderRefinementCollisionTestGlobal::MIBs = MIBs;
		renderRefinementCollisionTestGlobal::v_edges = v_edges;
		renderRefinementCollisionTestGlobal::voronoiBoundary = voronoiBoundary;
		renderRefinementCollisionTestGlobal::VRINs = VRINs;

		// build testers;
		renderRefinementCollisionTestGlobal::testers.resize(ms::numofframe);
		auto t0 = clock();
		for (int i = 0; i < ms::numofframe; i++)
		{
			cout << "********build BVH : " << i << endl;
			//dbgBlock[0] = i;
			renderRefinementCollisionTestGlobal::testers[i].initialize(VRINs[i]);
		}
		auto t1 = clock();
		cout << "BVH 360 constrcution : " << t1 - t0 << endl;

		ms::t2 = 185;

		auto reshapeFunc = [](GLint w, GLint h) ->void
		{
			double ratioWH = double(wd) / ht;
			if (w < h * ratioWH) {
				wd = w;
				ht = w / ratioWH;
				//glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
			}
			else {
				wd = h * ratioWH;
				ht = h;
				//glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
			}
		};
		auto    idleFunc = []() -> void
		{
			//manage speed
			static clock_t lastTime = clock();
			double timeInterval = 1000.0 / 36; // 10.0 * 1000 / ms::numofframe; //milliseconds
			while (clock() - lastTime < timeInterval)
			{

			}
			lastTime = clock();

			// managing input
			if (planning::keyboardflag['f'])	// toggle with f
				renderRefinementCollisionTestGlobal::sliceIdx++;

			static bool lastG = false, lastH = false;
			if (planning::keyboardflag['g'] != lastG)
				renderRefinementCollisionTestGlobal::sliceIdx++;
			if (planning::keyboardflag['h'] != lastH)
				renderRefinementCollisionTestGlobal::sliceIdx--;
			lastG = planning::keyboardflag['g'];
			lastH = planning::keyboardflag['h'];

			if (renderRefinementCollisionTestGlobal::sliceIdx >= ms::numofframe) renderRefinementCollisionTestGlobal::sliceIdx = 0;
			if (renderRefinementCollisionTestGlobal::sliceIdx < 0) renderRefinementCollisionTestGlobal::sliceIdx = ms::numofframe - 1;

			cout << "Current slice Number : " << renderRefinementCollisionTestGlobal::sliceIdx << endl;
			glutPostRedisplay();
		};
		auto displayFunc = [](void) -> void
		{
			using namespace renderRefinementCollisionTestGlobal;
			glClearColor(1.0, 1.0, 1.0, 1.0);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			setup_viewvolume();
			setup_transform();

			// Rendering code....
			glLineWidth(2.0f);
			glPointSize(2.8f);
			glColor3f(0.0f, 0.0f, 0.0f);

			// build temp

			auto sceneIdx = 7;
			auto robotIdx = 1;

			// VIEWPORT 0 (Upper-Left) : draw scene
			glViewport(0, 0, wd * 1 / 3, ht / 2);
			auto& scene = Models_Approx[sceneIdx];
			glColor3f(0, 0, 0);
			for (auto& as : scene)
				for (auto& arc : as.Arcs)
					arc.draw();

			// VIEWPORT 1 (Lower-Left) : darw robot
			glViewport(0, ht / 2, wd * 1 / 3, ht / 2);
			auto& robot = Models_Approx[robotIdx];
			glColor3f(0, 0, 0);
			for (auto& as : robot)
				for (auto& arc : as.Arcs)
				{
					cd::rotateArc(arc, sliceIdx).draw();
				}

			// VIEWPORT 2 (RIGHT) : draw mink/voronoi
			glViewport(wd * 1 / 3, 0, wd * 2 / 3, ht);
			auto& MR = renderRefinementCollisionTestGlobal::MRs[sliceIdx];
			auto& MIB = renderRefinementCollisionTestGlobal::MIBs[sliceIdx];
			auto& voronoi = renderRefinementCollisionTestGlobal::v_edges[sliceIdx];
			auto& boundary = renderRefinementCollisionTestGlobal::voronoiBoundary;
			auto& vrin = renderRefinementCollisionTestGlobal::VRINs[sliceIdx];

			static int interest = -1;
			static bool lc, lv;
			if (lc != planning::keyboardflag['c'])
				interest++;
			if (lv != planning::keyboardflag['v'])
				interest--;
			lc = planning::keyboardflag['c'];
			lv = planning::keyboardflag['v'];
			cout << "interest : " << interest << endl;

			// 1.5 draw interest;
			if (interest > -1 && interest < vrin.arcs.size())
			{
				glLineWidth(5.0f);
				glColor3f(0, 1, 0);
				auto& arc = vrin.arcs[interest];
				arc.draw();
				glLineWidth(2.0f);
			}

			// 1. draw mink
			{
				size_t length = MR.size();
				for (size_t i = 0; i < length; i++)
				{

					for (size_t j = 0; j < MR[i].size(); j++)
					{
						if (j == interest)
							glColor3f(0, 1, 0);
						else
						{
							if (MIB[i])
								glColor3f(0, 0, 0);
							else
								glColor3f(0, 0, 1);
						}

						MR[i][j].draw();
					}
				}
			}

			// 2. draw voronoi edge
			if(planning::keyboardflag['8'])
			{
				glColor3f(1, 0, 1);
				glBegin(GL_LINES);
				for (auto& ve : voronoi)
				{
					glVertex3d(ve.v0.x(), ve.v0.y(), -0.5);
					glVertex3d(ve.v1.x(), ve.v1.y(), -0.5);
				}
				glEnd();
			}

			// 3. draw Boundary for voronoi edges
			glColor3f(0, 0, 0);
			for (auto& arc : boundary)
				arc.draw();

			// 4. test collision
			if(planning::keyboardflag['5'])
			{
				using namespace renderRefinementCollisionTestGlobal;
				using ve = planning::output_to_file::v_edge;
				static int t = 650;
				static bool last9 = false, last0 = false;
				{
					if (planning::keyboardflag['9'] != last9)
						t += 10;
					if (planning::keyboardflag['0'] != last0)
						t -= 10;

					last9 = planning::keyboardflag['9'];
					last0 = planning::keyboardflag['0'];

					if (t >= renderRefinementCollisionTestGlobal::v_edges[ms::t2].size())
						t = renderRefinementCollisionTestGlobal::v_edges[ms::t2].size() - 1;

					if (t < 0)
						t = 0;
				}
				vector<ve> ves = renderRefinementCollisionTestGlobal::v_edges[ms::t2];
				vector<bool> added(ves.size(), false);

				// build testve
				vector<ve> testve;
				testve.push_back(ves[t]);
				added[t] = true;

				int n = 1000;
				for (int i = 0; i < n; i++)
				{
					Point& p = testve.back().v1;

					int closestidx = -1;
					double closestDist = 1e100;

					for (int j = 0; j < ves.size(); j++)
					{
						if (added[j])
							continue;

						double dist = (ves[j].v0 - p).length2();
						if (dist < closestDist)
						{
							closestidx = j;
							closestDist = dist;
						}

						dist = (ves[j].v1 - p).length2();
						if (dist < closestDist)
						{
							closestidx = j;
							closestDist = dist;
						}
					}

					if (closestidx == -1 || sqrt(closestDist) >0.01)
						break;
					else
					{
						if((ves[closestidx].v0 - p).length2() < (ves[closestidx].v1 - p).length2())
							testve.push_back(ves[closestidx]);
						else
						{
							auto temp = ves[closestidx];
							swap(temp.v0, temp.v1);
							testve.push_back(temp);
						}
						added[closestidx] = true;
					}
				}

				// do collision
				auto t0 = clock();

				//DEBUG : adhoc to test
				cout << ms::t2 % 2 << endl;
				if (ms::t2 % 2 == 0 && testve.size() > 15)
				{
					cout << testve.size() << endl;
					testve.erase(testve.begin(), testve.begin() + 5);
					testve.erase(testve.end() - 5, testve.end());
					cout << testve.size() << endl;

				}
				else if (ms::t2 % 2 == 0 && testve.size() > 5)
				{
					cout << testve.size() << endl;
					testve.erase(testve.begin(), testve.begin() + 1);
					testve.erase(testve.end() - 1, testve.end());
					cout << testve.size() << endl;
				}

				bool lineCol = testers[ms::t2].test(testve.begin(), testve.end());
				auto t1 = clock();

				if (lineCol)
					glColor3f(1, 0, 0);
				else
					glColor3f(0, 0, 1);

				auto p0 = testve.front().v0;
				auto p1 = testve.back().v1;
				
				glBegin(GL_LINES);
				glVertex2dv(p0.P);
				glVertex2dv(p1.P);

				glLineWidth(5.0f);
				glColor3f(0, 1, 0);
				for (auto& v : testve)
				{
					glVertex2dv(v.v0.P);
					glVertex2dv(v.v1.P);
				}
				glLineWidth(2.0f);
				glEnd();

				cout << "collisionTest case : " << t << endl;
				cout << "collisionTest result : " << lineCol << endl;
				cout << "collisionTest time : " << t1 - t0 << endl;
				cout << "testve.size() : " << testve.size() << endl;
				cout << "p0 : " << testve.front().v0 << endl;
			}

			// test2
			static double offset = 0.001;
			static bool lz, lx;
			if (lz != planning::keyboardflag['z'])
				offset += 0.005;
			if (lx != planning::keyboardflag['x'])
				offset -= 0.005;
			lz = planning::keyboardflag['z'];
			lx = planning::keyboardflag['x'];

			if(planning::keyboardflag['7'])
			{
				size_t length = 40;
				double dist = 4.0 / length;
				for (size_t i = 0; i < length; i++)
				{
					for (size_t j = 0; j < length; j++)
					{
						double x = i * dist - 2;
						double y = j * dist - 2;

						Point p = Point(x, y);
						Point q1 = Point(x + dist, y-offset);
						Point q2 = Point(x+offset, y + dist);

						bool
							res1 = testers[ms::t2].test(p, q1),
							res2 = testers[ms::t2].test(p, q2);

						glBegin(GL_LINES);
						
						if (res1)
							glColor3f(1, 0, 0);
						else
							glColor3f(0, 0, 1);
						glVertex2dv(p.P);
						glVertex2dv(q1.P);

						if (res2)
							glColor3f(1, 0, 0);
						else
							glColor3f(0, 0, 1);
						glVertex2dv(p.P);
						glVertex2dv(q2.P);

						glEnd();
					}
				}
			}

			// test3
			if (planning::keyboardflag['6'])
			{
				int j = 16, i = 16;

				size_t length = 40;
				double dist = 4.0 / length;
				{
					double x = i * dist - 2;
					double y = j * dist - 2;

					Point p = Point(x, y);
					Point q1 = Point(x + dist, y);
					Point q2 = Point(x, y + dist);

					bool
						res1 = true, // testers[ms::t2].test(p, q1),
						res2 = testers[ms::t2].test(p, q2);

					glBegin(GL_LINES);

					if (res1)
						glColor3f(1, 0, 0);
					else
						glColor3f(0, 0, 1);
					glVertex2dv(p.P);
					glVertex2dv(q1.P);

					if (res2)
						glColor3f(1, 0, 0);
					else
						glColor3f(0, 0, 1);
					glVertex2dv(p.P);
					glVertex2dv(q2.P);

					glEnd();
				}
			}

			glutSwapBuffers();
		};

		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
		glutInitWindowSize(wd, ht);
		glutCreateWindow("Data Check");
		glutReshapeFunc(reshapeFunc);
		glutDisplayFunc(displayFunc);
		glutMouseFunc(mouse_callback);
		glutKeyboardFunc(keyboard_callback);
		glutIdleFunc(idleFunc);

		glEnable(GL_DEPTH_TEST);

		glutMainLoop();
	}


}