#pragma once

#define BOUNDARY_EXTRACTOR_DO_ERROR_CHECK true

#pragma region include/classProto/alias/using
#include <vector>
#include  <set>;
namespace ms
{
	class CircularArc;
	class Point;
}
using std::vector;
using std::set;
using std::pair;
using ms::CircularArc;
using ms::Point;
#pragma endregion


/*
Def:
	given a set of circular arcs (probably overlaping loops), extract boundaries
Assume:
	all loops are counter-clockwise
Ret:
*/
class bex
{
	using inputType = vector<CircularArc>;
	using outputType = vector<vector<CircularArc>>;

public:
	bex() = default;
	~bex() = default;

	inline void setInput (inputType&  _in)  { this->_in = &_in; }
	inline void setOutput(outputType& _out) { this->_out = &_out; }

	void calculate();
private:
	inputType
		* _in;
	outputType
		* _out;

};
using boundaryExtractor = bex; 
using boundaryExtractionCalculator = bex;