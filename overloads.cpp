#include "overloads.h"
#include <math.h>
#include "scene.h"

using namespace std;
using namespace scene;

/* Operator Overloads */
/* VectorFloatTriplet */
VectorFloatTriplet crossProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; 
}

float dotProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return a.x * b.x + a.y * b.y + a.z * b.z; 
}

VectorFloatTriplet normalize(const VectorFloatTriplet& a) { 
    float length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z); 
    return VectorFloatTriplet{a.x / length, a.y / length, a.z / length}; 
}

VectorFloatTriplet operator-(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x - b.x, a.y - b.y, a.z - b.z}; 
}

VectorFloatTriplet operator+(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x + b.x, a.y + b.y, a.z + b.z}; 
}

VectorFloatTriplet operator*(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x * b.x, a.y * b.y, a.z * b.z}; 
}

VectorFloatTriplet& operator+=(VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}

VectorFloatTriplet& operator-=(VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
    return a;
}

VectorFloatTriplet& operator*=(VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
    return a;
}

VectorFloatTriplet& operator/=(VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
    return a;
}

VectorFloatTriplet operator-(const VectorFloatTriplet& a) { 
    return VectorFloatTriplet{-a.x, -a.y, -a.z}; 
}

VectorFloatTriplet operator*(const VectorFloatTriplet& a, const float& b) { 
    return VectorFloatTriplet{a.x * b, a.y * b, a.z * b}; 
}

VectorFloatTriplet operator*(const float& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a * b.x, a * b.y, a * b.z}; 
}
/* VectorIntTriplet */
VectorIntTriplet crossProduct(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; 
}

int dotProduct(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return a.x * b.x + a.y * b.y + a.z * b.z; 
}

VectorIntTriplet normalize(const VectorIntTriplet& a) { 
    int length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z); 
    return VectorIntTriplet{a.x / length, a.y / length, a.z / length}; 
}

VectorIntTriplet operator-(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x - b.x, a.y - b.y, a.z - b.z}; 
}

VectorIntTriplet operator+(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x + b.x, a.y + b.y, a.z + b.z}; 
}

VectorIntTriplet operator*(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x * b.x, a.y * b.y, a.z * b.z}; 
}

VectorIntTriplet& operator+=(VectorIntTriplet& a, const VectorIntTriplet& b) { 
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}

VectorIntTriplet& operator-=(VectorIntTriplet& a, const VectorIntTriplet& b) { 
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
    return a;
}

VectorIntTriplet& operator*=(VectorIntTriplet& a, const VectorIntTriplet& b) { 
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
    return a;
}

VectorIntTriplet& operator/=(VectorIntTriplet& a, const VectorIntTriplet& b) { 
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
    return a;
}
