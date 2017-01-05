#include "vector3.h"

Vector3::Vector3()
{
    set_to_null();
}

Vector3::Vector3(double f){
    vec3[0] = f;
    vec3[1] = f;
    vec3[2] = f;
}

Vector3::Vector3(double x, double y, double z){
    vec3[0] = x;
    vec3[1] = y;
    vec3[2] = z;
}

Vector3::Vector3(const double *vec3){
    this->vec3[0] = vec3[0];
    this->vec3[1] = vec3[1];
    this->vec3[2] = vec3[2];
}

Vector3::Vector3(Vector3 const& vec){
    this->vec3[0] = vec.vec3[0];
    this->vec3[1] = vec.vec3[1];
    this->vec3[2] = vec.vec3[2];
}

void Vector3::set_to_null(){
    vec3[0] = 0.0;
    vec3[1] = 0.0;
    vec3[2] = 0.0;
}

bool Vector3::is_null() const{
    return vec3[0] == 0.0 && vec3[1] == 0.0 && vec3[2] == 0.0;
}

void Vector3::set_value(int index, double value){
    if(index < 0 && index > 2){
        //qDebug("void Vector3::set_value(int index, double value) has a wrong index: %i", index);
        return;
    }
    vec3[index] = value;
}

double Vector3::get_value(int index){
    if(index < 0 && index > 2){
        //qDebug("double Vector3::get_value(int index) has a wrong index: %i", index);
        return 0;
    }
    return vec3[index];

}

double Vector3::x() const{
    return vec3[0];
}

double Vector3::y() const{
    return vec3[1];
}

double Vector3::z() const{
    return vec3[2];
}



void Vector3::set_x(double x){
    vec3[0] = x;
}

void Vector3::set_y(double y){
    vec3[1] = y;
}

void Vector3::set_z(double z){
    vec3[2] = z;
}



double Vector3::length() const{
    return sqrt(vec3[0] * vec3[0] + vec3[1] * vec3[1] + vec3[2] * vec3[2]);
}

double Vector3::lengthSquared() const{
    return vec3[0] * vec3[0] + vec3[1] * vec3[1] + vec3[2] * vec3[2];
}

Vector3 Vector3::normalized() const{
    // Need some extra precision if the length is very small.
    double len = double(vec3[0]) * double(vec3[0]) +
                 double(vec3[1]) * double(vec3[1]) +
                 double(vec3[2]) * double(vec3[2]);

    return *this / sqrt(len);

}

void Vector3::normalize(){
    // Need some extra precision if the length is very small.
    double len = double(vec3[0]) * double(vec3[0]) +
                 double(vec3[1]) * double(vec3[1]) +
                 double(vec3[2]) * double(vec3[2]);

     len = sqrt(len);

     vec3[0] /= len;
     vec3[1] /= len;
     vec3[2] /= len;
}


double Vector3::distance(Vector3 vector){
    Vector3 diff(vec3[0] - vector[0],
                vec3[1] - vector[1],
                vec3[2] - vector[2]);
    return diff.length();
}



const double& Vector3::operator[](int index) const{
    if(index < 0 && index > 2){
        //qDebug("const double& Vector3::operator[](int index) const has a wrong index: %i", index);
    }
    return vec3[index];
}

double& Vector3::operator[](int index){
    if(index < 0 && index > 2){
        //qDebug("double& Vector3::operator[](int index) has a wrong index: %i", index);
    }
    //var can be modified after return...
    return vec3[index];
}


const double& Vector3::operator()(int index) const{
    if(index < 0 && index > 2){
        //qDebug("const double& Vector3::operator()(int index) const has a wrong index: %i", index);
    }
    return vec3[index];
}

double& Vector3::operator()(int index){
    if(index < 0 && index > 2){
        //qDebug("double& Vector3::operator()(int index) has a wrong index: %i", index);
    }
    //var can be modified after return...
    return vec3[index];
}


Vector3& Vector3::operator+=(const Vector3 &vector){
    vec3[0] += vector.vec3[0];
    vec3[1] += vector.vec3[1];
    vec3[2] += vector.vec3[2];
    return *this;
}

Vector3& Vector3::operator-=(const Vector3 &vector){
    vec3[0] -= vector.vec3[0];
    vec3[1] -= vector.vec3[1];
    vec3[2] -= vector.vec3[2];
    return *this;
}

Vector3& Vector3::operator*=(double factor){
    vec3[0] *= factor;
    vec3[1] *= factor;
    vec3[2] *= factor;
    return *this;
}

Vector3& Vector3::operator*=(const Vector3 &vector){
    vec3[0] *= vector.vec3[0];
    vec3[1] *= vector.vec3[1];
    vec3[2] *= vector.vec3[2];
    return *this;
}

Vector3& Vector3::operator/=(double divisor){
    vec3[0] /= divisor;
    vec3[1] /= divisor;
    vec3[2] /= divisor;
    return *this;
}


double Vector3::dotProduct(const Vector3& v1, const Vector3& v2){
    return v1.vec3[0] * v2.vec3[0] + v1.vec3[1] * v2.vec3[1] +
           v1.vec3[2] * v2.vec3[2];
}

Vector3 Vector3::crossProduct(const Vector3& v1, const Vector3& v2)
{
    return Vector3(v1[1] * v2[2] - v1[2] * v2[1],
                    v1[2] * v2[0] - v1[0] * v2[2],
                    v1[0] * v2[1] - v1[1] * v2[0]);
}

Vector3 Vector3::normal(const Vector3& v1, const Vector3& v2)
{
    return crossProduct(v1, v2).normalized();
}


Vector3 Vector3::normal(const Vector3& v1, const Vector3& v2, const Vector3& v3)
{
    return crossProduct((v2 - v1), (v3 - v1)).normalized();
}


//friend
bool operator==(const Vector3 &v1, const Vector3 &v2){
    return v1.vec3[0] == v2.vec3[0] && v1.vec3[1] == v2.vec3[1] &&
           v1.vec3[2] == v2.vec3[2];
}

bool operator!=(const Vector3 &v1, const Vector3 &v2){
    return v1.vec3[0] != v2.vec3[0] || v1.vec3[1] != v2.vec3[1] ||
           v1.vec3[2] != v2.vec3[2];
}

const Vector3 operator+(const Vector3 &v1, const Vector3 &v2){
    return Vector3(v1.vec3[0] + v2.vec3[0] , v1.vec3[1] + v2.vec3[1] ,
                   v1.vec3[2] + v2.vec3[2]);
}

const Vector3 operator-(const Vector3 &v1, const Vector3 &v2){
    return Vector3(v1.vec3[0] - v2.vec3[0] , v1.vec3[1] - v2.vec3[1] ,
                   v1.vec3[2] - v2.vec3[2]);
}

const Vector3 operator*(double factor, const Vector3 &vector){
    return Vector3(vector.vec3[0] * factor , vector.vec3[1] * factor ,
                   vector.vec3[2] * factor);
}

const Vector3 operator*(const Vector3 &vector, double factor){
    return Vector3(vector.vec3[0] * factor , vector.vec3[1] * factor ,
                   vector.vec3[2] * factor);
}

const Vector3 operator*(const Vector3 &v1, const Vector3& v2){
    return Vector3(v1.vec3[0] * v2.vec3[0] , v1.vec3[1] * v2.vec3[1] ,
                   v1.vec3[2] * v2.vec3[2]);
}

const Vector3 operator-(const Vector3 &vector){
    return Vector3(-vector.vec3[0], -vector.vec3[1],
                   -vector.vec3[2]);
}

const Vector3 operator/(const Vector3 &vector, double divisor){
    return Vector3(vector.vec3[0] / divisor , vector.vec3[1] / divisor ,
                   vector.vec3[2] / divisor);
}
