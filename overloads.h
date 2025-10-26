#ifndef __OVERLOADS__
#define __OVERLOADS__

#include "scene.h"
#include "utils.h"

using namespace std;
using namespace scene;

/* Operator Overloads */
/* VectorFloatTriplet */
VectorFloatTriplet crossProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
float dotProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet normalize(const VectorFloatTriplet& a);
VectorFloatTriplet operator-(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator+(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator*(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet& operator+=(VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet& operator-=(VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet& operator*=(VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet& operator/=(VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator-(const VectorFloatTriplet& a);
VectorFloatTriplet operator*(const VectorFloatTriplet& a, const float& b);
VectorFloatTriplet operator*(const float& a, const VectorFloatTriplet& b);

/* VectorIntTriplet */
VectorIntTriplet crossProduct(const VectorIntTriplet& a, const VectorIntTriplet& b);
int dotProduct(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet normalize(const VectorIntTriplet& a);
VectorIntTriplet operator-(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet operator+(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet operator*(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet& operator+=(VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet& operator-=(VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet& operator*=(VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet& operator/=(VectorIntTriplet& a, const VectorIntTriplet& b);

#endif
