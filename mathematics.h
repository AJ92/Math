#ifndef MATH_H
#define MATH_H

#include <cmath>

// Pre-calculated value of PI / 180.
#define kPI180   0.017453

// Pre-calculated value of 180 / PI.
#define k180PI  57.295780

// Converts degrees to radians.
#define degreesToRadians(x) (x * kPI180)

// Converts radians to degrees.
#define radiansToDegrees(x) (x * k180PI)

#include "Vector/vector3.h"
#include "Vector/vector4.h"

#include "Matrix/matrix3x3.h"
#include "Matrix/matrix4x4.h"



#endif // MATH_H
