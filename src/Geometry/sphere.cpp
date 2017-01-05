#include "Geometry/sphere.h"

Sphere::Sphere()
{
    pos = Vector3(-1,-1,-1);
    radius = -1;
}

Sphere::Sphere(Vector3 pos, double radius) :
    pos(pos),
    radius(radius)
{

}

Sphere::Sphere(double pos_x, double pos_y, double pos_z, double radius) :
    pos(Vector3(pos_x,pos_y,pos_z)),
    radius(radius)
{

}


Vector3 Sphere::getPos(){
    return pos;
}

double Sphere::getRadius(){
    return radius;
}


void Sphere::setPos(Vector3 pos){
    this->pos = pos;
}

void Sphere::setPos(double pos_x, double pos_y, double pos_z){
    this->pos = Vector3(pos_x, pos_y, pos_z);
}

void Sphere::setRadius(double radius){
    this->radius = radius;
}
