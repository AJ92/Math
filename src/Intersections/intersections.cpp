#include "Intersections/intersections.h"

Intersections::Intersections()
{
}

double Intersections::square(double a){
    return a*a;
}

//static
bool Intersections::sphereAABBIntersection(Sphere s, AABB a){

    Vector3 sphere_pos = s.getPos();
    double sphere_radius = s.getRadius();

    Vector3 aabb_bmin = a.getBmin();
    Vector3 aabb_bmax = a.getBmax();

    double r2 = sphere_radius * sphere_radius;
    double dmin = 0.0;

    //3 dimensions
    for( int i = 0; i < 3; i++ ) {
      if( sphere_pos[i] < aabb_bmin[i] ) {
          dmin += square( sphere_pos[i] - aabb_bmin[i] );
      }
      else if( sphere_pos[i] > aabb_bmax[i] ) {
          dmin += square( sphere_pos[i] - aabb_bmax[i] );
      }
    }
    return dmin <= r2;
}

bool Intersections::pointInAABB(Vector3 point , AABB a){
    Vector3 aabb_bmin = a.getBmin();
    Vector3 aabb_bmax = a.getBmax();

    int checks = 0;
    //3 dimensions
    for( int i = 0; i < 3; i++ ) {
        if( point[i] >= aabb_bmin[i] ) {
            checks += 1;
        }
        if( point[i] <= aabb_bmax[i] ) {
            checks += 1;
        }
    }
    return checks == 6;
}

bool Intersections::pointInSphere(Vector3 point , Sphere s){
    Vector3 sphere_pos = s.getPos();
    double sphere_radius = s.getRadius();

    return (sphere_pos - point).length() <= sphere_radius;
}
