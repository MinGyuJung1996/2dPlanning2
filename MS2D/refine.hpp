#pragma once
/*****************************************************************************************************************************************************
**																																					**
**																Refinements																			**
**																																					**
**	Support:																																		**
**		Refine models non-g1																														**
**																																					**
*****************************************************************************************************************************************************/

#pragma region Other files
namespace ms
{
	class CircularArc;
	class Point;
}
using ms::CircularArc;
using ms::Point;

#include <vector>
#include <iostream>
#include <iomanip>

using std::vector;


#pragma endregion


/*

*/
class Refiner
{
public:
	static void g1(vector<CircularArc>& _in, vector<CircularArc>& _out);
	static void g0order(vector<CircularArc>& _in, vector<CircularArc>& _out); //TODO
	static void g0cut(vector<CircularArc>& _in, vector<CircularArc>& _out);	  //TODO


	inline static double& g1Error() { return _h_g1_error; }
	inline static bool&   beVerbose() { return _be_verbose; }
private:
	static double _h_g1_error;
	static bool   _be_verbose;
	static double _h_t0_begin;
	static double _h_t1_end;
};