#ifndef INTERSECTIONS_H
#define INTERSECTIONS_H

#include "Math/Geometry/aabb.h"
#include "Math/Geometry/sphere.h"

class Intersections
{
public:
    Intersections();

    static double square(double a);
    static bool sphereAABBIntersection(Sphere s, AABB a);
    static bool pointInAABB(Vector3 point , AABB a);
    static bool pointInSphere(Vector3 point , Sphere s);

};

#endif // INTERSECTIONS_H
