#include "distance field.hpp"
#include "stb_image_write.h"

using namespace std;
using namespace ms;

///*****************************************************************************************************************************************************
//**																																					**
//**															  Distance Field Calculator																**
//**																																					**
//*****************************************************************************************************************************************************/
//
///*
//Def: 
//	Per-pixel distance field calc
//Warning:
//	Uses render buffer. (uses glclear & other
//Desc:
//	Save info of Pixel's
//		1. closest obj
//		2. dist
//	by using depth buffer
//Assume:
//	loop is (ccw/cw)
//	loop is g1
//	obj insider [-1, 1] x [-1, 1]
//	No. of arc < 256 * 256 (2byte i guess)
//How to use:
//	1. make scene of circular arcs (globally ccw)
//		They should be inside [-1, 1] x [-1, 1]
//			if not, uniformly resize all of them and scale distance values at final step
//	2. set mat(no perspective) & other render specific stuffs if needed
//	3. setInput(vec<CircularArc>& )
//	4. calc()
//		(does rendering)
//	5. call multiple times: getDist(x, y, ...)
//		if (resizing was done in step 1, same should be done to x,y prior to call)
//		note that only x,y in [-1, 1] x [-1, 1] returns value
//*/
//class distFldCalc
//{
//public:
//	 distFldCalc() = default;
//	~distFldCalc() = default;
//
//	void setInput(vector<CircularArc>& _in);
//	void calc();
//	int  getDist(int	_in_w, int	  _in_h, double& _out_dist, int& _out_arcIdx);
//	int  getDist(double _in_x, double _in_y, double& _out_dist, int& _out_arcIdx);
//	//void setOutput();
//
//private:
//
//	bool _is_inputSet = false;
//	bool _is_calcDone = false;
//
//	int _width, _height; // width, height of render buffer
//	vector<CircularArc>* _in_ptr;
//	
//	vector<unsigned char> _renderBufferColor;
//	vector<float>		  _renderBufferDepth;
//
//	double _h_zBufferScale = 3.0; // zbuffer has dist, but it is confined to [-1, 1]. Therefore it needs some scaling
//		// (value read from z-buffer) * _h_zBufferScale = actual dist
//	double _h_backSign = -1.0; // rather 1.0 or -1.0. For camera configs, where obj behind has higher/lower z-value
//	double _h_sampleDtheta = PI2 / 360.0; // cone is approx with tri-fan. divide theta of circular arcs with this value (during tessellation)
//
//	void _drawCone(CircularArc& _in_arc, unsigned char* _in_colorInside, unsigned char* _in_colorOutside);
//
//};
//using distanceFieldCalculator = distFldCalc;
//using dfCalc = distFldCalc;


/*
Def: simply set intput
*/
void distFldCalc::setInput(vector<CircularArc>& _in, int _in_width, int _in_height)
{
	_in_ptr = &_in;
	_width = _in_width;
	_height = _in_height;
	_is_inputSet = true;
}

/*
Def: 
	find dist info per pixel
Assume:
	camera Matrix is identity. (the alg is dependent on camera mat)
Desc:
	Render cones with
		Color = arcIdx or errorCode/Type
			(0,0,0) : clear color
			(1, n/256, n%256) : n is closest arc idx;
			(255,255,255) : inside color
		
*/
void distFldCalc::calc()
{
	if (!_is_inputSet)
	{
		cerr << "Error: distFldCalc: input not set" << endl;
		return;
	}
	auto& in = *_in_ptr;

	// 1. clear
	glClearColor(rcode::r_notset / 255.0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity(); // err: not identity
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-1, 1, -1, 1);
	glViewport(0, 0, _width, _height);

	// 2. draw cone from arcs
	if(false)
	{
		unsigned char inColor[4] = { 255, 255, 255, 255 };

		for (int i = 0; i < in.size(); i++)
		{
			// 2-1. color
			unsigned char outColor[4];
			outColor[0] = 1;
			outColor[1] = i / 256;
			outColor[2] = i % 256;
			outColor[3] = 255;

			// 2-2. draw
			_drawCone(in[i], inColor, outColor);
		}
	}
	colorVoronoi(true);

	// 3. retrieve data = build _render buffer
	{
		_renderBufferColor.resize(4 * _width * _height);
		_renderBufferDepth.resize(_width * _height);

		
		glReadBuffer(GL_BACK);
		glReadPixels(0, 0, _width, _height, GL_RGBA, GL_UNSIGNED_BYTE, &_renderBufferColor[0]);
		glReadPixels(0, 0, _width, _height, GL_DEPTH_COMPONENT, GL_FLOAT, &_renderBufferDepth[0]);
	}

	if (false)
	{
		stbi_write_png((string("vorIdx.png")).c_str(), _width, _height, 4, &_renderBufferColor[0], 4 * _width);
	}

	// 99. isCalc done
	_is_calcDone = true;
}

/*
Def: from calculated data, find dist per pixel

Return:
	true iff success
*/
int  distFldCalc::getDist(int _in_w, int _in_h, double& _out_dist, int& _out_arcIdx)
{
	// -1. alias
	auto w = _in_w;
	auto h = _in_h;
	auto& d = _out_dist;
	auto& arcIdx = _out_arcIdx;

	// 0. check error
	if (!_is_calcDone)
	{
		// dbg_out
		cerr << "Error: distFldCalc: input not set" << endl;
		calc();
	}

	// 1. check bound
	if (
		w < 0 ||
		h < 0 ||
		w >= _width ||
		h >= _height
		)
	{
		// dbg_out
		cerr << "ERROR: out of bound at w, h = " << _in_w << ", " << _in_h << endl;
		return false;
	}

	// 2. find idx
	int dIdx = w + _width * (_height - h - 1); // WARNING : diff btw pixel coord and pixelRead
	int rIdx = 4 * dIdx + 0;
	int gIdx = 4 * dIdx + 1;
	int bIdx = 4 * dIdx + 2;

	// 3. find data
	//d = abs(-2 * _renderBufferDepth[dIdx] + 1.0) * _h_zBufferScale; // ortho maps [1, -1] to [0, 1], undo this, then, consider z buff scaling
	d = _depthToDist(_renderBufferDepth[dIdx]);
	unsigned char
		r = _renderBufferColor[rIdx],
		g = _renderBufferColor[gIdx],
		b = _renderBufferColor[bIdx];

	if (r == 1)
	{
		// dbg_out
		arcIdx = (256 * int(g)) + int(b);
		cerr << "rgb: " << int(r) << " " << int(g) << " " << int(b) << endl;
		cerr << "dist: " << d << "    arcIdx: " << arcIdx << endl;
		
		return true;
	}
	else
	{
		// dbg_out
		cerr << "rgb: " << int(r) << " " << int(g) << " " << int(b) << endl;
		cerr << "dist: " << d << endl;
		return false;
	}
}

/*
Def: from calculated data, find dist per pixel

Return:
	true iff success
*/
int  distFldCalc::getDist(double _in_x, double _in_y, double& _out_dist, int& _out_arcIdx)
{
	// alias
	auto x = _in_x;
	auto y = _in_y;

	// 0. check eror
	if (!_is_calcDone)
	{
		// dbg_out
		cerr << "Error: distFldCalc: input not set" << endl;
		calc();
	}

	// 1. check bound
	if (abs(_in_x) > 1.0 || abs(_in_y) > 1.0)
	{
		// dbg_out
		cerr << "ERROR: out of bound at x, y = " << _in_x << ", " << _in_y << endl;
		return false;
	}

	// 2. transform [-1, 1] x [-1, 1] -> [0, w] x [0, h]
	int pw, ph;
	{
		pw = (0.5 + 0.5 * x) * _width;
		ph = (0.5 - 0.5 * y) * _height;
	
	}
	
	// 3. call other func
	return getDist(pw, ph, _out_dist, _out_arcIdx);
}

/*
Def: 
	draw cone corresponding to a arc

Assume:
	arc is globally ccw
		used in relationship between color of (in/out) and (convex/concave)
In:
	arc: arc of interest
	cInside: color of insideRegion
	cOutside: color of outsideRegion
*/
void distFldCalc::_drawCone(CircularArc& _in_arc, unsigned char* _in_colorInside, unsigned char* _in_colorOutside)
{
	// alias
	auto& arc = _in_arc;
	bool convex = arc.lccw(); //arc's convexity
	unsigned char
		* colorConvex,
		* colorConcave;
	{
		if (convex)
		{
			colorConvex = _in_colorOutside;
			colorConcave = _in_colorInside;
		}
		else
		{
			colorConvex = _in_colorInside;
			colorConcave = _in_colorOutside;
		}
	}


	// 1. find theta (of sampling);
	vector<double> theta;
	{
		auto par = arc.param();
		auto t0 = par.first;
		auto t1 = par.second;
		
		// change so that t0 < t1 (as we already set variable "convex", we can flip arcs for convenience)
		if (t0 > t1)
		{
			auto temp = t0;
			t0 = t1;
			t1 = temp;
		}

		theta.reserve(abs(t1 - t0) / _h_sampleDtheta + 2.0);

		// 1-1. [t0, t1)
		for (auto t = t0; t < t1; t += _h_sampleDtheta)
		{
			theta.push_back(t);
		}
		
		// 1-2. t1
		theta.push_back(t1);
	}

	// 2. build tri
	Point apex;
	vector<Point> circ0, circ1;
	double apexZ, circ0Z, circ1Z; // z-value corresponding to each point
	// desc: tri-fan of (apex, circ0) forms a cone for concave region
	// desc: tri-strip of (circ0, circ1) forms a trun-cone for convex region
	{
		// 2-1. apex
		apex  = arc.cc();
		apexZ = arc.cr() / _h_zBufferScale * _h_backSign;

		// 2-2. circ0
		for (auto& t : theta)
		{
			circ0.push_back(arc.cc() + arc.cr() * Point(cos(t), sin(t)));
		}
		circ0Z = 0.0;

		// 2-3. circ1
		double rad = (1.0 * _h_zBufferScale) + arc.cr();
		for (auto& t : theta)
		{
			circ1.push_back(arc.cc() + rad * Point(cos(t), sin(t)));
		}
		circ1Z = _h_backSign * 1.0;
	}

	// 3. draw
	// 3-1. concave
	{
		glColor3ubv(colorConcave);
		glBegin(GL_TRIANGLE_FAN);
		glVertex3d(apex.x(), apex.y(), apexZ);
		for (auto& p : circ0)
			glVertex3d(p.x(), p.y(), circ0Z);
		glEnd();
	}
	// 3-2. convex
	{
		glColor3ubv(colorConvex);
		glBegin(GL_TRIANGLE_STRIP);
		for (int i = 0; i < theta.size(); i++)
		{
			glVertex3d(circ0[i].x(), circ0[i].y(), circ0Z);
			glVertex3d(circ1[i].x(), circ1[i].y(), circ1Z);
		}
		glEnd();
	}
}

/*
Def: 
	Choose colors and fill voronoi cells (uses input set data)
	if (useIdxAsColors), arc index data is encoded to colors
	if (!useIdxAsColors), random colors are used for arcs
*/
void distFldCalc::colorVoronoi(bool useIdxAsColors)
{
	if(!useIdxAsColors)
	{
		for (int i = 0; i < _in_ptr->size(); i++) 
		{
			// dbg
			if (dbgVar >= 0 && dbgVar < _in_ptr->size())
			{
				if (i != dbgVar)
					continue;
				else
				{
					cout << (*_in_ptr)[i] << endl;
				}
			}

			unsigned char r = _red(i);
			unsigned char g = _green(i);
			unsigned char b = _blue(i);

			unsigned char c[8] =
			{ r, g, b, 255,
			255 - r, 255 - g, 255 - b, 255 };
			_drawCone((*_in_ptr)[i], c, c + 4);
		}
	}
	else
	{
		for (int i = 0; i < _in_ptr->size(); i++)
		{
			// dbg
			if (dbgVar >= 0 && dbgVar < _in_ptr->size())
			{
				if (i != dbgVar)
					continue;
				else
				{
					cout << (*_in_ptr)[i] << endl;
				}
			}

			unsigned char r = rcode::r_valid;
			unsigned char g = i / 256;
			unsigned char b = i % 256;

			unsigned char c[8] =
			{ rcode::r_inside, g, b, 255,
			r, g, b, 255 };
			_drawCone((*_in_ptr)[i], c, c + 4);
		}
	}
}

/*
Def:
	find pixel-wise distance map use red ~ blue ~ green coloring
*/
void distFldCalc::findDistanceMap(vector<unsigned char>& _out_rgbaBuffer)
{
	constexpr double _h_colorChangeDist = 0.3;
	constexpr double _h_colorChangeDist2 = 1.0;
	float color0[4] = { 1,0,0,1 };
	float color1[4] = { 0,0,1,1 };
	float color2[4] = { 0,1,0,1 };

	// 1. alias & resize
	auto& out = _out_rgbaBuffer;
	out.resize(_renderBufferColor.size());

	// 2. find color
	for (int i = 0; i < _renderBufferColor.size(); i += 4)
	{
		auto& r = out[i + 0];
		auto& g = out[i + 1];
		auto& b = out[i + 2];
		auto& a = out[i + 3];

		// exception
		if (_renderBufferColor[i] != rcode::r_valid)
		{
			r = g = b = a = 255;
			continue;
		}

		// find dist
		auto depth = _renderBufferDepth[i/4];
		auto dist = _depthToDist(depth);

		// if(certain interval) interpolate
		if (dist < _h_colorChangeDist)
		{
			double t = dist / _h_colorChangeDist;
			r = 255 * ((1 - t) * color0[0] + (t)*color1[0]);
			g = 255 * ((1 - t) * color0[1] + (t)*color1[1]);
			b = 255 * ((1 - t) * color0[2] + (t)*color1[2]);
			a = 255 * ((1 - t) * color0[3] + (t)*color1[3]);
		}
		else if (dist < _h_colorChangeDist2)
		{
			double t = (dist - _h_colorChangeDist) / (_h_colorChangeDist2 - _h_colorChangeDist);
			r = 255 * ((1 - t) * color1[0] + (t)*color2[0]);
			g = 255 * ((1 - t) * color1[1] + (t)*color2[1]);
			b = 255 * ((1 - t) * color1[2] + (t)*color2[2]);
			a = 255 * ((1 - t) * color1[3] + (t)*color2[3]);
		}
		else
		{
			r = 255 * color2[0];
			g = 255 * color2[1];
			b = 255 * color2[2];
			a = 255 * color2[3];
		}
	}


}

/*
Def: draw stuffs from findDistanceMap
*/
void distFldCalc::drawDistanceMap(vector<unsigned char>& _in_rgbaBuffer)
{
	glDrawPixels(_width, _height, GL_RGBA, GL_UNSIGNED_BYTE, &(_in_rgbaBuffer[0]));
}

/*
Def:
	find pixel-wise distance map use gray coloring + inner region coloring
	color varies loop-wise
*/
void distFldCalc::findDistanceMap2(vector<unsigned char>& _out_rgbaBuffer, vector<int>& _in_arcIdxToLoopIdx)
{
	constexpr double _h_colorChangeDist  = 0.2;
	constexpr double _h_colorChangeDist2 = 1.0;
	constexpr double mid = 0.8;
	float color0[4] = { 0.0, 0.0, 0.0, 1.0 };
	float color1[4] = { mid, mid, mid, 1.0 };
	float color2[4] = { 1.0, 1.0, 1.0, 1.0 };

	// 1. alias & resize
	auto& out = _out_rgbaBuffer;
	out.resize(_renderBufferColor.size());

	// 2. find color
	for (int i = 0; i < _renderBufferColor.size(); i += 4)
	{
		auto& r = out[i + 0];
		auto& g = out[i + 1];
		auto& b = out[i + 2];
		auto& a = out[i + 3];

		// exception
		if (_renderBufferColor[i] == rcode::r_inside)
		{
			int arcIdx = 256 * int(_renderBufferColor[i+1]) + int(_renderBufferColor[i+2]);
			int loopIdx = _in_arcIdxToLoopIdx[arcIdx];
			r = _red  (loopIdx);
			g = _green(loopIdx);
			b = _blue (loopIdx);
			a = 255;

			////dbg_out
			//if(arcIdx != 0)
			//	cout << "arcIdx, loopIdx: " <<(arcIdx) << " , " << (loopIdx) << endl;

			continue;
		}

		// find dist
		auto depth = _renderBufferDepth[i / 4];
		auto dist = _depthToDist(depth);

		// if(certain interval) interpolate
		if (dist < _h_colorChangeDist)
		{
			double t = dist / _h_colorChangeDist;
			r = 255 * ((1 - t) * color0[0] + (t)*color1[0]);
			g = 255 * ((1 - t) * color0[1] + (t)*color1[1]);
			b = 255 * ((1 - t) * color0[2] + (t)*color1[2]);
			a = 255 * ((1 - t) * color0[3] + (t)*color1[3]);
		}
		else if (dist < _h_colorChangeDist2)
		{
			double t = (dist - _h_colorChangeDist) / (_h_colorChangeDist2 - _h_colorChangeDist);
			r = 255 * ((1 - t) * color1[0] + (t)*color2[0]);
			g = 255 * ((1 - t) * color1[1] + (t)*color2[1]);
			b = 255 * ((1 - t) * color1[2] + (t)*color2[2]);
			a = 255 * ((1 - t) * color1[3] + (t)*color2[3]);
		}
		else
		{
			r = 255 * color2[0];
			g = 255 * color2[1];
			b = 255 * color2[2];
			a = 255 * color2[3];
		}
	}


}

/*
Def:
	find pixel-wise distance map with uniform color regardless of coloring
	color varies loop-wise
*/
void distFldCalc::findDistanceMap3(vector<unsigned char>& _out_rgbaBuffer, vector<int>& _in_arcIdxToLoopIdx)
{
	//constexpr double _h_colorChangeDist = 0.3;
	//constexpr double _h_colorChangeDist2 = 1.0;
	//float color0[4] = { 1,0,0,1 };
	//float color1[4] = { 0,0,1,1 };
	//float color2[4] = { 0,1,0,1 };

	// 1. alias & resize
	auto& out = _out_rgbaBuffer;
	out.resize(_renderBufferColor.size());

	// 2. find color
	for (int i = 0; i < _renderBufferColor.size(); i += 4)
	{
		auto& r = out[i + 0];
		auto& g = out[i + 1];
		auto& b = out[i + 2];
		auto& a = out[i + 3];

		// exception
		if (_renderBufferColor[i] == rcode::r_inside)
		{
			r = g = b = a = 0;
			continue;
		}

		int arcIdx = 256 * int(_renderBufferColor[i + 1]) + int(_renderBufferColor[i + 2]);
		int loopIdx = _in_arcIdxToLoopIdx[arcIdx];
		r = _red  (loopIdx);
		g = _green(loopIdx);
		b = _blue (loopIdx);
		a = 255;
	}
}