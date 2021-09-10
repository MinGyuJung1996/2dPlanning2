#pragma once

/*
Def:
	A Header that is contained by all files
	Main purpose:
		take care of all hyperParameters in a file.
How to use:
	For files to reference all this var
		Just include
	For a file to actually store var definitions
		#define GLOBAL_HEADERS_DEFINITONS true
		#include "GlobalHeader.hpp"

*/

#ifndef GLOBAL_HEADER_DEFINITIONS
#define GLOBAL_HEADER_DEFINITIONS false
#endif

#if GLOBAL_HEADER_DEFINITIONS == true
#define DEFINE_PARAM(type, name, value) type name = value;
#else
#define DEFINE_PARAM(type, name, value) extern type name;
#endif


/*******************************************************************************************************************************
**																															  **
**														HYPER-PARAMS														  **
**																															  **
*******************************************************************************************************************************/

DEFINE_PARAM(int, _h_tester, 10)

namespace HYPER
{
	DEFINE_PARAM(double, POINT_ARC_RADIUS, 1e-12);



}