#ifndef OGRECOMMON_H
#define OGRECOMMON_H

#define _OgreExport

#include <iostream>
#include <cassert>

#define _USE_MATH_DEFINES
#include <cmath>

#define TWO_PI (M_PI*2)
#define HALF_PI (M_PI/2)
namespace Ogre
{
	class Vector2;
	class Vector3;
	class Vector4;
	class Matrix3;
	class Matrix4;
	class Quaternion;
}

inline float isqrtf(float x)
{
	return 1/sqrtf(x);
}

#endif