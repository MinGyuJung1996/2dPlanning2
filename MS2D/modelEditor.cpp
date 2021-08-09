#include "MS2D.h"
#include <iostream>
#include <iomanip>
#include "voronoi.hpp"
#include "collision detection.hpp"
#include "xyt-cs.hpp"
#include "refine.hpp"

#include <filesystem>
#include <sstream>
namespace fs = std::experimental::filesystem; // may need to be changed, since it moved out from experimental

#define COUT(x) #x " : " << x


void arcAlignedPathSweepVolume(vector<CircularArc>& _in_path, vector<CircularArc>& _in_robot, vector<CircularArc>& _out);

namespace modelEditor
{
	int wd, ht;
	// mouse
	Point mousePosition;
	bool isMousePositionAtEnd = false;
	bool isMouseAlignedToPrecedingArcs = false;
	vector<Point> mouseAlignedPoints;
	double mxOri, myOri; // current mouse pos in world-coord (original, not modified)
	int mouseState = 0;
		// 0 : currently picking arc start pt // only happens at start of program, as arcs are g0 cont
		// 1 : currently picking arc end   pt
		// 2 : currently picking arc rad
	Point x0, x1;

	// adjustable
	double samplingRate = 0.1;
	double _h_line_rad = 100.0;
	
	// mode : may need to change to enum later
	enum modes { modes_editing = 0, modes_cspace, modes_rsv };
	modes pgMode = modes_editing;

	// result, current edit
	vector<CircularArc> model_arcs;
	vector<Circle>		model_disk;

	// result, saved loops
	vector<vector<CircularArc>> models_arcs;
	vector<vector<Circle>>		models_disk;
	int loopIdx = -1;

	// 
	vector<CircularArc> background_arcs;

	//  cameraMatrix
	float cameraMatrix[16] = {
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1 };
	constexpr int
		cm_zoom = 0,
		cm_x = 12,
		cm_y = 13; // e.g.) camearMatrix[cm_y] := camera's y position

	namespace mode1
	{
		bool robotSet = false;
		bool obsSet = false;
		bool vbSet = false;
		int sliceNo = 0;

		string robotName;
		double robotScale;
	}


	// manual
	string manual =
		R"del(
~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Manual ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Make CCW loop and find interior disks

~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mode 0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~ Draw stuff

Mouse Left: draw
Mouse Right: undo
Mouse Scroll : zoom (TODO)
Arrow : camera position (TODO)
????? : camera reset (TODO)

// Keyboard Left (frequently used)
q: Align:	p2 is picked axis aligned
w: Align:	p2 is in a same line with some point in vec
r: Refine:	refine model g1
a: Axis:	draw axis grid
d: Disk:	calc disk
z: Change line rad
v: Voronoi:	vor test

// Keyboard Right (less frequently used)
p: save loop to loop-sets
\: clear saved loop-sets
-: file: save loop-sets to   file (type no to cancel)
=: file: load loop-sets from file (type no to cancel)
[: loopIdx-- 
]: loopIdx++ 
u,j: moveV  loop
h,k: moveH  loop
y,i: rotate loop
n,m: scale  loop 
<  : copy   loop 
>  : delete loop
;  : recalculate disk for loops
'  : set background display loop

// Numbers (Change visibility)
~: change mode
1: show disk calculated
2: print mouse positions

~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mode 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~ show mink vor


~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mode 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~ show rsv (path = currently drawing arcs)
)del";

	/*
	hyper params: (in this code section)
		guessArc::_h_line_rad

	hyper params: (outside)
		things related to regfine::g1
	*/


	/*
	Def: from glo var x0,x1 and mousPosition, guess arc

	Warning:
		this uses global variables.
	*/
	CircularArc guessArc()
	{
		auto& mp = mousePosition;

		// 1. check whethre line representing arc should be sent
		bool isLine = false;
		{
			// i. when mp is too close
			if ((mp - x0).length2() < 1e-8 || (mp - x1).length2() < 1e-8)
				isLine = true;
			else // if no else -> mp-x0 could be 0 -> normalize -> nan -> propagtes strangely with comparator
			{
				// ii. check whether ordering (x0-mp-x1) holds in arc
				auto tan = (x1 - x0).normalize();
				auto d0 = x0 * tan;
				auto d1 = x1 * tan;
				auto d2 = mp * tan;

				bool f0 = d0 < d2;
				bool f1 = d1 < d2;
				if (!(f0 ^ f1))
					isLine = true;
				else
				{
				// iii. when dot product of (mp-x0), (mp-x1) is too big
				auto tan0 = (mp - x0).normalize();
				auto tan1 = (x1 - mp).normalize();
				if (abs(tan0 * tan1) > 1 - 1e-2)
					isLine = true;
				}
			}

		}

		// 
		if (isLine)
		{
			return cd::constructArc(x0, x1, _h_line_rad, true, true);
		}
		else
		{
			auto l0 = ms::bisector(x0, mp);
			auto l1 = ms::bisector(x1, mp);
			Point cc(l0, l1);
			bool ccw;
			{
				// find ccw
				// travel from t0 in ccw
				//	-> if tt is met first than t1 <==> ccw

				// find param
				auto n0 = x0 - cc;
				auto n1 = x1 - cc;
				auto nt = mp - cc;
				auto t0 = atan2(n0.y(), n0.x());
				auto t1 = atan2(n1.y(), n1.x());
				auto tt = atan2(nt.y(), nt.x());

				// confine it to [t0, t0+2pi)
				if (t1 < t0)
					t1 += PI2;
				if (tt < t0)
					tt += PI2;
				
				// check whether traveling from t0 to t1 confronts tt
				if (tt < t1)
					ccw = true;
				else
					ccw = false;
			}

			return cd::constructArc(cc, x0, x1, ccw );
		}

	}

	/*
	Def: tranlsates loop
	*/
	void 
	translateLoop(
		Point _in_translation,
		vector<CircularArc>& _in_arcs,
		vector<Circle>&		 _in_crcs,
		vector<CircularArc>& _out_arcs,
		vector<Circle>&		 _out_crcs)
	{
		for (auto& a : _in_arcs)
			_out_arcs.push_back(cd::translateArc(a, _in_translation));

		for (auto& c : _in_crcs)
		{
			_out_crcs.push_back(c);
			_out_crcs.back().c = _out_crcs.back().c + _in_translation;
		}
	}

	/*
	Def: rotates loop
	*/
	void
		rotateLoop(
			double _in_rotationDegree,
			vector<CircularArc>& _in_arcs,
			vector<Circle>& _in_crcs,
			vector<CircularArc>& _out_arcs,
			vector<Circle>& _out_crcs)
	{
		if (_in_arcs.size() == 0) return;
		
		//
		auto rotRad = _in_rotationDegree / 180.0 * PI;
		
		// get rotation center, some point that won't change with rotation
		Point p(0, 0);
		for (auto& a : _in_arcs)
			p = p + a.x0();
		p = p / _in_arcs.size();
		Point mp = -p;

		// rotate arc
		for (auto& a : _in_arcs)
			_out_arcs.push_back(
				cd::translateArc(
					cd::rotateArc(
						cd::translateArc(
							a,
							mp), 
						_in_rotationDegree),
					p)
			);

		// rotate crc
		for (auto& c : _in_crcs)
		{
			Point temp = c.c;
			temp = temp - p;
			temp = temp.rotate(rotRad);
			temp = temp + p;

			_out_crcs.push_back(c);
			_out_crcs.back().c = temp;
		}

	}

	/*
	Def: rotates loop
	*/
	void
		scaleLoop(
			double _in_scale,
			vector<CircularArc>& _in_arcs,
			vector<Circle>& _in_crcs,
			vector<CircularArc>& _out_arcs,
			vector<Circle>& _out_crcs)
	{
		if (_in_arcs.size() == 0) return;

		// get scale center(=rotation center), some point that won't change with rotation
		Point p(0, 0);
		for (auto& a : _in_arcs)
			p = p + a.x0();
		p = p / _in_arcs.size();
		Point mp = -p;

		// scale arc
		for (auto& a : _in_arcs)
			_out_arcs.push_back(
				cd::translateArc(
					cd::scaleArc(
						cd::translateArc(
							a,
							mp),
						_in_scale),
					p)
			);

		//
		for (auto& c : _in_crcs)
		{
			Point temp = c.c;
			temp = temp - p;
			temp = temp * _in_scale;
			temp = temp + p;

			_out_crcs.push_back(c);
			_out_crcs.back().c = temp;
			_out_crcs.back().r = c.r * _in_scale;
		}


	}

	/*
	Def: 
		set robot for mode_Cspace
	Ref:
		initializeRobotObstacles
	*/
	void setRobot(string _in_name, double _in_scale)
	{
		using namespace ms;
		// 0. Trivial 
		// Just set some Idx. By default 1: robot, 7: Obstacles
		ModelInfo_CurrentModel.first = 1; // 1
		ModelInfo_CurrentModel.second = 7; // 6

		int rsvIdx = 1, // mink-sum-call idx of robot's RSV 
			oriIdx = 0, // mink-sum-call idx of original robot
			obsIdx = 7; // mink-sum-call idx of obstacles

		// 1-0. Alias
		using namespace std;
		Model_from_arc[rsvIdx] = true;
		Model_from_arc[oriIdx] = true;

		// These 4 variables are the goal of this section
		auto& robOri = Model_vca[oriIdx];		// circArc output (ori)
		auto& robRsv = Model_vca[rsvIdx];		// circArc output (rsv)
		auto& robVAS = Models_Approx[rsvIdx];	// Arcspline output
		auto& robCir = InteriorDisks_Imported[rsvIdx]; // circ output

		// 1-1. init step (set path and scalefactor)
		string path("Robots/");
		double scaleFactor = _in_scale;
		path += _in_name;
		path += "/";

		auto pathRsv(path + "arcRSV.txt");
		auto pathOri(path + "arcOri.txt");
		auto pathCircRsv(path + "circ.txt");
		auto pathCircOri(path + "circ.txt"); // if we intend to use diff circles for arcOri, this part should be changed

		vector<Circle>
			circRead0,
			circRead1,
			garbage;
		vector<CircularArc>
			arcRead0,
			arcRead1,
			temp;

		// 1-2. Do reading
		readArcModel(pathOri.c_str(), pathCircOri.c_str(), arcRead0, circRead0);
		readArcModel(pathRsv.c_str(), pathCircRsv.c_str(), arcRead1, circRead1);

		// 1-3. Apply transformations
		appendArcModel(temp, robCir, arcRead1, circRead1, scaleFactor, 0, Point(0, 0)); // RSV case
		appendArcModel(robOri, garbage, arcRead0, circRead0, scaleFactor, 0, Point(0, 0)); // Original case
		// do RSV

		// 1-4. make temp G0
		vector<CircularArc> g0;
		planning::_Convert_VectorCircularArc_G0(temp, g0, 1);

		// 1-5. make g0 G1
		vector<CircularArc>& g1 = robRsv;
		if (useNewerVersionOfG1conveter)
			planning::_Convert_VectorCircularArc_G1_new(g0, g1);
		else
			planning::_Convert_VectorCircularArc_G1(g0, g1);

		// 1-6. convert to AS
		vector<ArcSpline>& asRSV = robVAS;
		planning::_Convert_VectorCircularArc_To_MsInput(g1, asRSV);

		// 1-7. postprocess
		ms::postProcess(ModelInfo_CurrentModel.first, ModelInfo_CurrentModel.second);
	}

	/*
	Def: 
		set obstacles to what we've made
	Ref:
		initializeRobotObstacles
	*/
	void setObstacles()
	{
		using namespace ms;
		// 0. Trivial 
		// Just set some Idx. By default 1: robot, 7: Obstacles
		ModelInfo_CurrentModel.first = 1; // 1
		ModelInfo_CurrentModel.second = 7; // 6

		int rsvIdx = 1, // mink-sum-call idx of robot's RSV 
			oriIdx = 0, // mink-sum-call idx of original robot
			obsIdx = 7; // mink-sum-call idx of obstacles
		
		// 2. Obstacles
		{
			/*
			What is done:
			1. set flag : Model_from_arc
			2. set model data : Model_vca (will be used in post processing... optional since this is obstacles, not robot...)
			3. convert model and save : Models_approx
			4. also set disks : ID_I

			Notice that there can be up to 2~3 copies of the same model with different representation
				=> Model_vca(vector<CircularArc>, original, read from file)
				=> Models_approx(vector<ArcSpline>, modified(segmented) to be an input of minkowski sum)
				=> Models_approx_rotated(rotate Model_vca and then segment)
			*/

			// 1. alias
			Model_from_arc[obsIdx] = true;
			auto& sceneOriginal = Model_vca[obsIdx];
			auto& sceneProcessed = Models_Approx[obsIdx];
			auto& sceneCircles = InteriorDisks_Imported[obsIdx];

			sceneOriginal.resize(0);
			sceneProcessed.resize(0);
			sceneCircles.resize(0);


			// 2. Load models
			for (auto& loop : models_arcs)
				for (auto& arc : loop)
					sceneOriginal.push_back(arc);
			for (auto& disks : models_disk)
				for (auto& disk : disks)
					sceneCircles.push_back(disk);

			// 3. build vec<AS> sceneProcessed(input to MinkSum)
			planning::_Convert_VectorCircularArc_To_MsInput(sceneOriginal, sceneProcessed);
		}
	}

	void motionFunc(int x, int y)
	{
		// 1. find mousePosition
		mxOri = float(x - wd / 2) / (wd / 2);
		myOri = float(-y + ht / 2) / (ht / 2);
		mousePosition.x() = mxOri;
		mousePosition.y() = myOri;

		isMousePositionAtEnd = false;
		isMouseAlignedToPrecedingArcs = false;
		mouseAlignedPoints.clear();

		// 2. if certain button is pushed, align position
		//cout << COUT(Key('q')) << endl;
		if (Key('q'))
		{
			// if choosing second point of arc
			if (mouseState == 1)
			{
				// i. align x? y?
				bool alignXdir = false;
				auto difVec = mousePosition - x0;
				if (abs(difVec.x()) > abs(difVec.y()))
					alignXdir = true;

				// ii. set
				if (alignXdir)
					mousePosition.y() = x0.y();
				else
					mousePosition.x() = x0.x();
			}
		}

		if (Key('w') && model_arcs.size() > 0 && mouseState == 1)
		{
			//auto& p = model_arcs.front().x0(); 
			//if (abs(mousePosition.x() - p.x()) < 5e-3)
			//{
			//	mousePosition.x() = p.x();
			//	isMouseAlignedToPrecedingArcs = true;
			//}
			//if (abs(mousePosition.y() - p.y()) < 5e-3)
			//{
			//	mousePosition.y() = p.y();
			//	isMouseAlignedToPrecedingArcs = true;
			//}

			for (auto& arc : model_arcs)
			{
				auto& p = arc.x0();
				if (abs(mousePosition.x() - p.x()) < 5e-3)
				{
					mousePosition.x() = p.x();
					isMouseAlignedToPrecedingArcs = true;
					mouseAlignedPoints.push_back(p);
				}
				if (abs(mousePosition.y() - p.y()) < 5e-3)
				{
					mousePosition.y() = p.y();
					isMouseAlignedToPrecedingArcs = true;
					mouseAlignedPoints.push_back(p);
				}
			}

		}


		// if close to start of loop, make it so.
		if (model_arcs.size() > 0)
		{
			if ((mousePosition - model_arcs.front().x0()).length2() < 1e-4)
			{
				mousePosition = model_arcs.front().x0();
				isMousePositionAtEnd = true;
			}
		}

	};
	void reshapeFunc(GLint w, GLint h)
	{
		// later
		if (w > h)
		{
			wd = h;
			ht = h;
		}
		else
		{
			wd = w;
			ht = w;
		}
	};
	void idleFunc()
	{
		// later

		glutPostRedisplay();
	};
	void mouseFunc(int button, int action, int x, int y)
	{
		if(pgMode == modes_editing)
		{
		if(action == GLUT_UP)
		switch (button)
		{
		case GLUT_LEFT_BUTTON:
		{
			switch(mouseState)
			{
			case 0: //initial state : where no arcs are in model_arcs
				// save x0
				mouseState = 1;
				x0 = mousePosition;
				break;
			case 1:
				// save x1
				mouseState = 2;
				x1 = mousePosition;
				break;
			case 2:
				// from x0, x1, and current 
				mouseState = 1;

				//find arc
				model_arcs.push_back(guessArc());

				// next arc's starting point = this arc's end
				x0 = x1; 
				break;
			}
		}
			break;

		case GLUT_RIGHT_BUTTON:
		{
			switch (mouseState)
			{
			case 0:
				if (models_arcs.size() > 0)
				{
					model_arcs = models_arcs.back();
					model_disk.clear();

					models_arcs.pop_back();
					models_disk.pop_back();

					mouseState = 1;
					x0 = model_arcs.back().x1();

				}
				break;
			case 1:
				if (model_arcs.size() > 0)
				{
					x0 = model_arcs.back().x0();
					model_arcs.pop_back();
					mouseState = 1;
				}
				else
				{
					mouseState = 0;
				}
				break;
			case 2:
				mouseState = 1;
				break;
			}
		}
			break;
		}
		}
		else if (pgMode == modes_cspace)
		{

		}
	};
	void keyboardFunc(unsigned char a, int b, int c)
	{
		//cout << "keyboardFunc:" << a << endl;
		planning::keyboardflag[a] = !planning::keyboardflag[a];
		// pressing a once will flip keybaord flag btw 0 and 1
		// pressing it for a long time will make it oscillate 010101
		
	};
	void keyboardUpFunc(unsigned char a, int b, int c)
	{

	}
	void displayFunc()
	{
		// 0.
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glViewport(0, 0, wd, ht);

		glMatrixMode(GL_MODELVIEW);
		glLoadMatrixf(cameraMatrix);

		// 0-1. check change of mode
		if (Key3('~'))
		{
			switch (pgMode)
			{
			default:
			case modes::modes_editing:
				mode1::obsSet = false;
				pgMode = modes::modes_cspace;
				break;
			case modes::modes_cspace:
				pgMode = modes::modes_rsv;
				break;
			case modes::modes_rsv:
				pgMode = modes::modes_editing;
				break;
			}
		}

		if (pgMode == modes_editing)
		{
		if(Key('2'))
			cout << "current mous Pos : " << mousePosition << endl;

		// 1. show added
		{
			glColor3f(0, 0, 0);
			for (auto& a : model_arcs)
				a.draw2();

			int cnt = 0;
			for (auto& l : models_arcs)
			{
				if (cnt == loopIdx)
					glColor3f(0, 0, 1);
				else
					glColor3f(0, 0, 0);

				for (auto& a : l)
					a.draw2();

				cnt++;
			}
		}

		// 2. show current arc
		switch (mouseState)
		{
		case 0:
			break;
		case 1:
			glColor3f(0.8, 0.8, 0.8);
			glBegin(GL_LINES);
			glVertex2dv(mousePosition.P);
			glVertex2dv(x0.P);
			glEnd();
			break;
		case 2:
			glColor3f(0.8, 0.8, 0.8);
			auto arc = guessArc();
			arc.draw2();
			break;
		}

		// 3. show mousPosition (& align pts)
		{
			if (isMousePositionAtEnd)
				glColor3f(0, 1, 0);
			else if (isMouseAlignedToPrecedingArcs)
				glColor3f(0, 1, 1);
			else
				glColor3f(0, 0, 1);
			
			Circle c;
			c.c = mousePosition;
			c.r = 0.01;

			c.draw();

			{
				glBegin(GL_LINES);
				for (auto& p : mouseAlignedPoints)
				{
					glVertex2dv(mousePosition.P);
					glVertex2dv(p.P);
				}
				glEnd();
			}
		}

		// 4. find disks
		if (Key3('d'))
		{
			if (model_disk.size() > 0)
				model_disk.clear();

			if (model_arcs.front().x0() == model_arcs.back().x1())
			{
				findInteriorDisks2(model_arcs, samplingRate, model_disk);
				cout << "disks: disks found " << endl;
			}
			else
			{
				cout << "disks: loop not formed " << endl;
			}
		}

		// 5. draw disks
		if (!Key('1'))
		{
			glColor3f(0.7, .7, .7);
			for (auto& c : model_disk)
				c.draw();

			int cnt = 0;
			for (auto& b : models_disk)
			{
				if (cnt == loopIdx)
					glColor3f(0.7, .7, 1);
				else
					glColor3f(0.7, .7, .7);

				for (auto& c : b)
					c.draw();

				cnt++;
			}

		}

		// 6. refine
		if (Key3('r'))
		{
			if (model_arcs.front().x0() == model_arcs.back().x1())
			{
				vector<CircularArc> temp;
				Refiner rf;
				rf.g1(model_arcs, temp);

				cout << "RefineMent: from " << model_arcs.size() << " to " << temp.size() << endl;

				//dbg_out
				{
				for (auto& a : model_arcs)
					cout << a << endl;
				for (auto& a : temp)
					cout << a << endl;
				}

				model_arcs = temp;
			}
		}

		// 7. axis
		if (Key('a'))
		{
			double z = 1.0;

			glColor3f(.7, .7, .7);
			
			glLineWidth(2.0);
			glBegin(GL_LINES);
			glVertex3d(-1.0, 0.0, z);
			glVertex3d(+1.0, 0.0, z);
			glVertex3d(0.0, -1.0, z);
			glVertex3d(0.0, +1.0, z);
			glEnd();


			glLineWidth(0.5);
			glBegin(GL_LINES);
			for (int i = 1; i < 20; i++)
			{
				if (i == 10) continue;
				
				if (i == 5 || i == 15)
					glColor3f(.8, .8, .8);
				else
					glColor3f(.9, .9, .9);

				double t = double(i - 10) * 0.1;

				glVertex3d(-1.0, t   , z);
				glVertex3d(+1.0, t   , z);
				glVertex3d(t   , -1.0, z);
				glVertex3d(t   , +1.0, z);
			}
			glEnd();

			glLineWidth(1.0);
		}

		// 8. test vor
		if (Key('v') && models_arcs.size() > 0)
		{
			planning::VR_IN vrin;

			// i. vrin models
			int offset = 0;
			for (auto& v : models_arcs)
			{
				// models_arcs is ccw => convert cw
				for (int i = 0; i < v.size(); i++)
				{
					int in = i + 1;
					if (in == v.size())
						in = 0;

					vrin.arcs.push_back(planning::flipArc(v[i]));
					vrin.left.push_back(in + offset);
					vrin.color.push_back(0.0);
				}

				offset += v.size();
			}

			// ii. vrin boundary
			{
				Point p0(1, -1);
				Point p1(1, 1);
				Point p2(-1, 1);
				Point p3(-1, -1);

				auto a0 = cd::constructArc(p0, p1, 10, true, true);
				auto a1 = cd::constructArc(p1, p2, 10, true, true);
				auto a2 = cd::constructArc(p2, p3, 10, true, true);
				auto a3 = cd::constructArc(p3, p0, 10, true, true);

				vrin.arcs.push_back(a0);
				vrin.arcs.push_back(a1);
				vrin.arcs.push_back(a2);
				vrin.arcs.push_back(a3);

				vrin.left.push_back(offset + 3);
				vrin.left.push_back(offset + 0);
				vrin.left.push_back(offset + 1);
				vrin.left.push_back(offset + 2);

				vrin.color.push_back(1.0);
				vrin.color.push_back(1.0);
				vrin.color.push_back(1.0);
				vrin.color.push_back(1.0);
			}

			// iii. do vor
			voronoiCalculator vc;
			vector<deque<VoronoiEdge>> res;

			vc.setInput(vrin.arcs, vrin.left, vrin.color);
			vc.setOutput(res);
			vc.calculate();

			// iv. draw
			glColor3f(1, 0, 1);
			glBegin(GL_LINES);
			for (auto& deq : res)
				for (auto& ve : deq)
				{
					glVertex2dv(ve.v0.P);
					glVertex2dv(ve.v1.P);
				}
			glEnd();


		}
		
		// 9. change line Rad
		if (Key3('z'))
		{
			cout << "Change line rad from " << _h_line_rad << " to: " << endl;
			double d;
			cin >> d;
			_h_line_rad = abs(d);
		}

		// 10. background
		if (Key3('\''))
		{
			cout << "name of background object: (type no to clear) " << endl;
			string name;
			cin >> name;

			if (name == "no")
			{
				background_arcs.clear();
			}
			else
			{
				if (name[0] == '/')
					name = string(name.c_str() + 1); // remove first char if ('/')
				if (name.back() == '/')
					name.pop_back();

				name = "modelEditor/" + name;

				string arcFileName, cirFileName;
				{
					arcFileName = name + "/arc.txt";
					cirFileName = name + "/circ.txt";
				}

				ifstream arcin(arcFileName), crcin(cirFileName);
				if (!arcin || !crcin)
				{
					cout << "Failed to read file" << endl;
				}
				else
				{
					// read data
					vector<Circle> cread;
					readArcModel(arcFileName.c_str(), cirFileName.c_str(), background_arcs, cread);
				}
			}
		}
		{
			glColor3f(0.8, 0.8, 0.8);
			for (auto& a : background_arcs)
			{
				a.draw2();
			}
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//																															  //
		//															loop related													  //
		//																															  //
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		{
			constexpr double
				dt = 0.01, // amount of translation
				dr = 3.00, // amount of rotation
				ds = 0.024; // amount of scaling

			// 20. clear
			if (Key3('\\'))
			{
				cout << "Type 'yes' to clear loops: " << endl;
				string temp;
				cin >> temp;
				if (temp == "yes")
				{
					models_arcs.clear();
					models_disk.clear();
				}
			}

			// 21. change loopIdx
			if (Key3('['))
			{
				loopIdx--;
				if (loopIdx < -1)
					loopIdx = -1; // no loop selected
			}

			if (Key3(']'))
			{
				loopIdx++;
				if (loopIdx == models_arcs.size())
					loopIdx--;
			}

			// 22. move current loop to models_arc
			if (Key3('p'))
			{
				// if(form loop)
				if (model_arcs.size() > 0 && model_arcs.front().x0() == model_arcs.back().x1())
				{

					vector<CircularArc> refined;
					vector<Circle> interiorDisks;

					// i. refine
					Refiner rf;
					rf.g1(model_arcs, refined);

					// ii. find disk
					findInteriorDisks2(refined, samplingRate, interiorDisks);

					// iii. push
					models_arcs.push_back(refined);
					models_disk.push_back(interiorDisks);

					// iv. flush
					model_arcs.clear();
					model_disk.clear();
					mouseState = 0;

					cout << "recording loop with size = " << refined.size() << ", disk = " << interiorDisks.size() << endl;
				}
			}

			// 23. save arc/disk to file
			if (Key3('-'))
			{
				// take in input name
				cout << "Saving to file. Type in new model name: " << endl;
				string name;
				cin >> name;

				if (name != "no")
				{

					if (name[0] == '/')
						name = string(name.c_str() + 1); // remove first char if ('/')
					if (name.back() == '/')
						name.pop_back();

					name = "modelEditor/" + name;

					// make directory
					fs::path relDir(name);
					bool succeed = fs::create_directories(name); // notice that existing folder name is not allowed

					// make file name
					string arcFileName, cirFileName;
					if (!succeed)
					{
						cerr << "Err: Writing to file failed during creation of diretory. Saving to Temp" << endl;

						arcFileName = "modelEditor/arc.txt";
						cirFileName = "modelEditor/circ.txt";
					}
					else
					{
						arcFileName = name + "/arc.txt";
						cirFileName = name + "/circ.txt";
					}

					// get size
					int totalArc = 0, totalCirc = 0;
					vector<int> lidx_arc, lidx_circ; // contiain starting arcIdx of loop
					{
						for (auto& v : models_arcs)
						{
							lidx_arc.push_back(totalArc);
							totalArc += v.size();
						}
						for (auto& v : models_disk)
						{
							lidx_circ.push_back(totalCirc);
							totalCirc += v.size();
						}
					}

					cout << "Writing (arc, disk) = (" << totalArc << ", " << totalCirc << ") to file" << endl;

					// write arc
					{
						ofstream arcout(arcFileName);
						arcout << std::setprecision(20);

						// write arc data
						arcout << totalArc << endl;
						for (auto& v : models_arcs)
							for (auto& a : v)
							{
								auto par = a.param();

								arcout << "d " << a.cx() << " " << a.cy() << " " << a.cr() << " "
									<< par.first * 180.0 / PI << " " << par.second * 180.0 / PI << " " << a.ccw << endl;
							}

						// write meta data
						arcout << endl
							<< "Metadata: " << endl
							<< "LoopIdx: ";
						for (auto i : lidx_arc)
							arcout << i << " ";

					}
					// write disk
					{
						ofstream circout(cirFileName);
						circout << std::setprecision(20);

						//write arc data
						circout << totalCirc << endl;
						for (auto& v : models_disk)
							for (auto& c : v)
							{
								circout << c.c.x() << " " << c.c.y() << " " << c.r << endl;
							}


						// write meta data
						circout << endl
							<< "Metadata: " << endl
							<< "LoopIdx: ";
						for (auto i : lidx_circ)
							circout << i << " ";
					}

				}
			}

			// 24. load file
			if (Key3('='))
			{
				cout << "Loading file. Type in existing model name: " << endl;
				string name;
				cin >> name;

				if (name != "no")
				{
					if (name[0] == '/')
						name = string(name.c_str() + 1); // remove first char if ('/')
					if (name.back() == '/')
						name.pop_back();

					name = "modelEditor/" + name;

					string arcFileName, cirFileName;
					{
						arcFileName = name + "/arc.txt";
						cirFileName = name + "/circ.txt";
					}

					ifstream arcin(arcFileName), crcin(cirFileName);
					if (!arcin || !crcin)
					{
						cout << "Failed to read file" << endl;
					}
					else
					{
						// read data
						vector<CircularArc> aread;
						vector<Circle> cread;
						readArcModel(arcFileName.c_str(), cirFileName.c_str(), aread, cread);

						//read metadata;
						
						// i. 
						vector<int> interval_a, interval_c;
						char buffer[256];
						while (!arcin.eof())
						{
							arcin.getline(buffer, 256);
							stringstream ss;
							ss << buffer;
						
							string firstword;
							ss >> firstword;
						
							if (firstword == "LoopIdx:")
							{
								while (ss)
								{
									int i;
									if (ss >> i) //assume space at end
										interval_a.push_back(i);
								}
							}
						}
						while (!crcin.eof())
						{
							crcin.getline(buffer, 256);
							stringstream ss;
							ss << buffer;
						
							string firstword;
							ss >> firstword;
						
							if (firstword == "LoopIdx:")
							{
								while (ss)
								{
									int i;
									if (ss >> i)
										interval_c.push_back(i);
								}
							}
						}
						
						//if (failed to read loop info)
						if (interval_a.size() == 0 || interval_c.size() == 0 || (interval_a.size() != interval_c.size()))
						{
							cout << "Appending the whole file as a loop: (arcNo, crcNo) = (" << aread.size() << ", " << cread.size() << ")" << endl;
						
							models_arcs.push_back(aread);
							models_disk.push_back(cread);
						}
						else
						{
							cout << "Appending [";
							for (auto i : interval_a)
								cout << i << ", ";
							cout << aread.size();
							cout << ") arcs." << endl;
						
							cout << "Appending [";
							for (auto i : interval_c)
								cout << i << ", ";
							cout << cread.size();
							cout << ") crcs." << endl;
						
							interval_a.push_back(aread.size());
							interval_c.push_back(cread.size());
						
							for (int i = 0; i < interval_a.size() - 1; i++) // size of interval_a,_c is same
							{
								int a0 = interval_a[i + 0];
								int a1 = interval_a[i + 1];
								int c0 = interval_c[i + 0];
								int c1 = interval_c[i + 1];
						
								if (a1 <= a0)
									continue;
								if (c1 <= c0)
									continue;
						
								vector<CircularArc> temp;
								vector<Circle> temp2;
								for (int a = a0; a < a1; a++)
									temp.push_back(aread[a]);
								for (int c = c0; c < c1; c++)
									temp2.push_back(cread[c]);
								models_arcs.push_back(temp);
								models_disk.push_back(temp2);
							}
						
						}
						
						// ii.
					} // if (!arcin || !crcin)
				} // ~if(!= "no")
			} // if (Key3('='))

			// 25. translate/rotate/scale
			{
				if (loopIdx >= 0 && loopIdx < models_arcs.size())
				{
					// 1. find T, R, S

					bool
						doTrn	= false,
						doRot	= false,
						doScl	= false;
					
					Point  T(0, 0);
					double R = 0.0;
					double S = 1.0;

					if (Key3('u') && !Key3('j')) { doTrn = true; T.y() += dt; }
					if (Key3('j') && !Key3('u')) { doTrn = true; T.y() -= dt; }
					if (Key3('h') && !Key3('k')) { doTrn = true; T.x() -= dt; }
					if (Key3('k') && !Key3('h')) { doTrn = true; T.x() += dt; }
					if (Key3('y') && !Key3('i')) { doRot = true; R	 += dr; }
					if (Key3('i') && !Key3('y')) { doRot = true; R	 -= dr; }
					if (Key3('m') && !Key3('n')) { doScl = true; S	 += ds; }
					if (Key3('n') && !Key3('m')) { doScl = true; S	 -= ds; }
						   
					// do T, R, S
					if (doTrn)
					{
						auto& va = models_arcs[loopIdx];
						auto& vc = models_disk[loopIdx];

						vector<CircularArc> tempa;
						vector<Circle> tempc;

						translateLoop(T, va, vc, tempa, tempc);

						va = tempa;
						vc = tempc;
					}
					if (doRot)
					{
						auto& va = models_arcs[loopIdx];
						auto& vc = models_disk[loopIdx];

						vector<CircularArc> tempa;
						vector<Circle> tempc;

						rotateLoop(R, va, vc, tempa, tempc);

						va = tempa;
						vc = tempc;
					}
					if (doScl)
					{
						auto& va = models_arcs[loopIdx];
						auto& vc = models_disk[loopIdx];

						vector<CircularArc> tempa;
						vector<Circle> tempc;

						scaleLoop(S, va, vc, tempa, tempc);

						va = tempa;
						vc = tempc;
					}
				}
			}

			// 26. copy/delete Loop
			if (loopIdx >= 0 && loopIdx < models_arcs.size())
			{
				if (Key3('<'))
				{
					auto& vca = models_arcs[loopIdx];
					auto& vc = models_disk[loopIdx];
					models_arcs.push_back(vca);
					models_disk.push_back(vc);

					Point trans(0.1, -0.1);
					for (auto& a : models_arcs.back())
					{
						a.cc() = a.cc() + trans;
						a.x0() = a.x0() + trans;
						a.x1() = a.x1() + trans;
					}
					for (auto& c : models_disk.back())
					{
						c.c = c.c + trans;
					}

					loopIdx = models_arcs.size() - 1;
				}

				if (Key3('>'))
				{
					auto it = models_arcs.begin() + loopIdx;
					models_arcs.erase(it);

					auto it2 = models_disk.begin() + loopIdx;
					models_disk.erase(it2);
				}
				
			}

			// 27. recalculate Disk for whole loops in loop-set (too support cw loops => holes inside model)
			if (Key3(';'))
			{

				//i. build vec<CA>
				vector<CircularArc> va;
				vector<int> li; // save loop idx
				for (int i = 0; i < models_arcs.size(); i++)
				{
					for (auto& a : models_arcs[i])
					{
						va.push_back(a);
						li.push_back(i);
					}
				}

				// ii. call
				vector<Circle> vc;
				vector<int> cidx;
				findInteriorDisksWithIdx(va, li, samplingRate, vc, cidx);

				// iii. assign
				for (auto& v : models_disk)
					v.clear();

				for (int i = 0; i < vc.size(); i++)
				{
					models_disk[cidx[i]].push_back(vc[i]);
				}
			}
			
		}
		}
		else if (pgMode == modes_cspace)
		{
			using namespace mode1;
			// 0. take care of first call
			if (mode1::robotSet == false)
			{
				robotName = "8-3";
				robotScale = 0.2;
				setRobot(robotName, robotScale);
				mode1::robotSet = true;
			}
			if (mode1::obsSet == false)
			{
				setObstacles();
				mode1::obsSet = true;
			}
			if (vbSet == false)
			{
				Point p0(1, -1);
				Point p1(1, 1);
				Point p2(-1, 1);
				Point p3(-1, -1);

				auto a0 = cd::constructArc(p0, p1, 10, true, true);
				auto a1 = cd::constructArc(p1, p2, 10, true, true);
				auto a2 = cd::constructArc(p2, p3, 10, true, true);
				auto a3 = cd::constructArc(p3, p0, 10, true, true);

				planning::voronoiBoundary.push_back(a0);
				planning::voronoiBoundary.push_back(a1);
				planning::voronoiBoundary.push_back(a2);
				planning::voronoiBoundary.push_back(a3);

				vbSet = true;
			}
			// 1. change sliceNo
			if (Key3('['))
			{
				sliceNo--;
				if (sliceNo < 0)
					sliceNo = 0;
				cout << COUT(sliceNo) << endl;
			}
			if (Key3(']'))
			{
				sliceNo++;
				if (sliceNo >= ms::numofframe)
					sliceNo = ms::numofframe - 1;
				cout << COUT(sliceNo) << endl;
			}

			// 2. draw mink vor
			{
				// 2-1.
				ms::t2 = sliceNo;
				ms::minkowskisum(ms::t2, 7);
				
				glColor3f(0, 0, 0);
				for (auto& a : ms::Model_Result)
					for (auto& as : a)
						for (auto& arc : as.Arcs)
							arc.draw();
				
				// 2-2. Vor
				planning::VR_IN vrin;
				planning::_Convert_MsOut_To_VrIn(ms::Model_Result, ms::ModelInfo_Boundary, vrin);
				
				voronoiCalculator vc;
				vc.initialize();
				vc.setInput(vrin.arcs, vrin.left, vrin.color);

				glColor3f(0, 0, 1);
				int vType = 1;
				if (vType == 1)
				{
					vector<deque<VoronoiEdge>> v_res;
					vc.setOutput(v_res);
					vc.calculate();

					glBegin(GL_LINES);
					for (auto& a : v_res)
						for (auto& b : a)
						{
							glVertex2dv(b.v0.P);
							glVertex2dv(b.v1.P);
						}
					glEnd();
					//cout << COUT(v_res.size()) << endl;
					//cout << COUT(vrin.arcs.size()) << endl;
				}
			}

			// 3. change robot
		}
		else if (pgMode == modes_rsv)
		{
			auto& robot = ms::Model_vca[1];
			//
			robotForward = Point(1, 1).normalize();

			auto path = model_arcs;

			vector<CircularArc> temp;

			arcAlignedPathSweepVolume(path, robot, temp);

			if (!Key('1'))
			{
				glColor3f(0, 0, 0);
				for (auto& a : path)
					a.draw2();
			}
			if (!Key('2'))
			{
				glColor3f(0, 0, 1);
				for (auto& a : temp)
					a.draw2();
			}
			if (!Key('3'))
			{
				glColor3f(0, 1, 0);
				auto& p = path;
				auto& r = robot;

				for (auto& a : p)
				{
					auto temp = r;

					auto t = atan2(robotForward.y(), robotForward.x());
					auto tan0 = a.tan0();
					auto tan1 = a.tan1();
					auto t0 = atan2(tan0.y(), tan0.x()) - t; // rot needed to
					auto t1 = atan2(tan1.y(), tan1.x()) - t; // rot needed to

					if (a.ccw)
					{
						while (t1 < t0)
							t1 += PI2;
					}
					else
					{
						while (t0 < t1)
							t0 += PI2;
					}

					// 1. rotate temp;
					vector<CircularArc> temp2;
					{
						for (auto& b : temp)
						{
							temp2.push_back(cd::rotateArc(b, t0 * 180 / PI));
						}
					}

					// 2. translate temp
					vector<CircularArc> temp3;
					{
						for (auto& b : temp2)
						{
							auto trans = a.x0() - a.cc();
							temp3.push_back(cd::translateArc(b, trans));
						}
					}
					
					// 3. rotate again
					vector<CircularArc> temp4;
					{
						for (auto& b : temp3)
						{
							temp4.push_back(cd::rotateArc(b, (t1-t0) * 180 / PI));
						}
					}

					// 4. draw
					for (auto& b : temp2)
						cd::translateArc(b, a.x0()).draw2();

					for (auto& b : temp4)
						cd::translateArc(b, a.cc()).draw2();

				}

			}

		}

		// 99.
		for (int i = 0; i < 256; i++)
		{
			//Key2(i) = Key(i);
			planning::keyboardflag_last[i] = planning::keyboardflag[i];
		}

		glutSwapBuffers();
	}

	int main(int argc, char* argv[])
	{
		// 0. print manual
		{
			cout << manual << endl;
		}

		// 1. basic stuff
		wd = 800;
		ht = 800;
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
		glutInitWindowSize(wd, ht);
		glutCreateWindow("SceneEditor");
		//glEnable(GL_SCISSOR_TEST);

		// 99. register Func and start loop
		glutReshapeFunc(reshapeFunc);
		glutDisplayFunc(displayFunc);
		glutMouseFunc(mouseFunc);
		glutKeyboardFunc(keyboardFunc);
		glutKeyboardUpFunc(keyboardUpFunc);
		glutIdleFunc(idleFunc);
		glutMotionFunc(motionFunc);
		glutPassiveMotionFunc(motionFunc);

		glutMainLoop();
		
		return 0;
	}
}