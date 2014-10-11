#include "aabb.h"

AABB::AABB()
{
}

AABB::AABB(Vector3 bmin, Vector3 bmax) :
    bmin(bmin),
    bmax(bmax)
{

}

AABB::AABB( double bmin_x, double bmin_y, double bmin_z,
            double bmax_x, double bmax_y, double bmax_z) :
    bmin(Vector3(bmin_x, bmin_y, bmin_z)),
    bmax(Vector3(bmax_x, bmax_y, bmax_z))
{

}

Vector3 AABB::getBmin(){
    return bmin;
}

Vector3 AABB::getBmax(){
    return bmax;
}

void AABB::setBmin(Vector3 bmin){
    this->bmin = bmin;
}

void AABB::setBmax(Vector3 bmax){
    this->bmax = bmax;
}

void AABB::setBmin(double bmin_x, double bmin_y, double bmin_z){
    this->bmin = Vector3(bmin_x, bmin_y, bmin_z);
}

void AABB::setBmax(double bmax_x, double bmax_y, double bmax_z){
    this->bmax = Vector3(bmax_x, bmax_y, bmax_z);
}
