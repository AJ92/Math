#include "plane.h"

Plane::Plane(){

}

Plane::Plane(Vector3 point1, Vector3 point2, Vector3 point3){
    //construct the plane

    Vector3 v = point2-point1;
    Vector3 u = point3-point1;

    plane_normal = Vector3::crossProduct(v,u);
    plane_normal.normalize();

    d = Vector3::dotProduct(-plane_normal,point1);
    plane_point = point1;
}

Plane::Plane(Vector3 point, Vector3 normal){
    //construct the plane
    plane_point = point;
    plane_normal = normal;
    d = Vector3::dotProduct(-plane_normal,point);
}

double Plane::distance(Vector3 toPoint){
    return Vector3::dotProduct(plane_normal,toPoint) + d;
}

Vector3 Plane::project(Vector3 fromPoint){
    return fromPoint - distance(fromPoint) * plane_normal;
}
