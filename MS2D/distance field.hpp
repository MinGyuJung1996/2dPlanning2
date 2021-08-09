#pragma once
#include <vector>
#include "MS2D.h"

using std::vector;
using ms::CircularArc;


/*****************************************************************************************************************************************************
**																																					**
**															  Distance Field Calculator																**
**																																					**
*****************************************************************************************************************************************************/

/*
Def:
	Per-pixel distance field calc
Warning:
	Uses render buffer. (uses glclear & other
Desc:
	Save info of Pixel's
		1. closest obj
		2. dist
	by using depth buffer
Assume:
	loop is (ccw/cw)
	loop is g1
	obj insider [-1, 1] x [-1, 1]
	No. of arc < 256 * 256 (2byte i guess)
How to use:
	1. make scene of circular arcs (globally ccw)
		They should be inside [-1, 1] x [-1, 1]
			if not, uniformly resize all of them and scale distance values at final step
	2. set mat(no perspective) & other render specific stuffs if needed
	3. setInput(vec<CircularArc>& )
	4. calc()
		(does rendering)
	5. call multiple times: getDist(x, y, ...)
		if (resizing was done in step 1, same should be done to x,y prior to call)
		note that only x,y in [-1, 1] x [-1, 1] returns value
*/
class distFldCalc
{
public:
	distFldCalc() = default;
	~distFldCalc() = default;

	void setInput(vector<CircularArc>& _in, int _in_width, int _in_height);
	void calc();

	void colorVoronoi(bool useIdxAsColors = false);
	
	int  getDist(int	_in_w, int	  _in_h, double& _out_dist, int& _out_arcIdx);
	int  getDist(double _in_x, double _in_y, double& _out_dist, int& _out_arcIdx);
	
	void findDistanceMap (vector<unsigned char>& _out_rgbaBuffer);
	void drawDistanceMap (vector<unsigned char>&  _in_rgbaBuffer);
	void findDistanceMap2(vector<unsigned char>& _out_rgbaBuffer, vector<int>& _in_arcIdxToLoopIdx);
	void findDistanceMap3(vector<unsigned char>& _out_rgbaBuffer, vector<int>& _in_arcIdxToLoopIdx);
	//void setOutput();

	int dbgVar = -1;

	// r-color of rendered image contains info of arc or error
	enum rcode
	{
		r_valid = 1,
		r_inside = 255, // color for inner region
		r_notset = 254 // clear color
	};

private:

	bool _is_inputSet = false;
	bool _is_calcDone = false;

	int _width, _height; // width, height of render buffer
	vector<CircularArc>* _in_ptr;

	vector<unsigned char> _renderBufferColor;
	vector<float>		  _renderBufferDepth;

	double _h_zBufferScale = 3.0; // zbuffer has dist, but it is confined to [-1, 1]. Therefore it needs some scaling
		// (value read from z-buffer) * _h_zBufferScale = actual dist
	double _h_backSign = -1.0; // rather 1.0 or -1.0. For camera configs, where obj behind has higher/lower z-value
	double _h_sampleDtheta = PI2 / 360.0; // cone is approx with tri-fan. divide theta of circular arcs with this value (during tessellation)

	void _drawCone(CircularArc& _in_arc, unsigned char* _in_colorInside, unsigned char* _in_colorOutside);
	inline double _depthToDist(double _in_depth) { return abs(-2 * _in_depth + 1.0) * _h_zBufferScale; }

/*debug*/ public:
	inline unsigned char _red  (int i) 
	{ 
		if (i == 0 || i == 3 || i == 4) return 255;
		else if (i < 5) return 0;
		else if (i == 13) return ((61 * i) % 256) * 2;
		else return (61 * i) % 256; 
	}
	inline unsigned char _green(int i) 
	{ 
		if (i == 1 || i == 3) return 255; 
		else if (i < 5) return 0; 
		else if (i == 13) return ((99 * i) % 256) * 2;
		else return (99 * i) % 256;
	}
	inline unsigned char _blue (int i) 
	{ 
		if (i == 2 || i == 4 ) return 255; 
		else if (i < 5) return 0; 
		else return (27 * i) % 256; 
	}
	/*
	1 0 0
	0 1 0
	0 0 1
	1 1 0
	1 0 1
	
	*/
};
using distanceFieldCalculator = distFldCalc;
using dfCalc = distFldCalc;