#ifndef SPHERE_H
#define SPHERE_H

#include "Math/mathematics.h"

class Sphere
{
public:
    Sphere();
    Sphere(Vector3 pos, double radius);
    Sphere(double pos_x, double pos_y, double pos_z, double radius);

    Vector3 getPos();
    double getRadius();

    void setPos(Vector3 pos);
    void setPos(double pos_x, double pos_y, double pos_z);
    void setRadius(double radius);

private:
    Vector3 pos;
    double radius;
};

#endif // SPHERE_H
