#include "vector4.h"

Vector4::Vector4()
{
    set_to_null();
}

Vector4::Vector4(double f){
    vec4[0] = f;
    vec4[1] = f;
    vec4[2] = f;
    vec4[3] = f;
}

Vector4::Vector4(double x, double y, double z, double w){
    vec4[0] = x;
    vec4[1] = y;
    vec4[2] = z;
    vec4[3] = w;
}

Vector4::Vector4(const double *vec4){
    this->vec4[0] = vec4[0];
    this->vec4[1] = vec4[1];
    this->vec4[2] = vec4[2];
    this->vec4[3] = vec4[3];
}

Vector4::Vector4(Vector4 const& vec){
    this->vec4[0] = vec.vec4[0];
    this->vec4[1] = vec.vec4[1];
    this->vec4[2] = vec.vec4[2];
    this->vec4[3] = vec.vec4[3];
}

void Vector4::set_to_null(){
    vec4[0] = 0.f;
    vec4[1] = 0.f;
    vec4[2] = 0.f;
    vec4[3] = 0.f;
}

bool Vector4::is_null() const{
    return vec4[0] == 0.0 && vec4[1] == 0.0 && vec4[2] == 0.0 && vec4[3] == 0.0;
}

void Vector4::set_value(int index, double value){
    if(index < 0 && index > 3){
        //qDebug("void Vector4::set_value(int index, double value) has a wrong index: %i", index);
        return;
    }
    vec4[index] = value;
}

double Vector4::get_value(int index){
    if(index < 0 && index > 3){
        //qDebug("double Vector4::get_value(int index) has a wrong index: %i", index);
        return 0;
    }
    return vec4[index];

}

double Vector4::x() const{
    return vec4[0];
}

double Vector4::y() const{
    return vec4[1];
}

double Vector4::z() const{
    return vec4[2];
}

double Vector4::w() const{
    return vec4[3];
}

void Vector4::set_x(double x){
    vec4[0] = x;
}

void Vector4::set_y(double y){
    vec4[1] = y;
}

void Vector4::set_z(double z){
    vec4[2] = z;
}

void Vector4::set_w(double w){
    vec4[3] = w;
}


double Vector4::length() const{
    return sqrt(vec4[0] * vec4[0] + vec4[1] * vec4[1] + vec4[2] * vec4[2] + vec4[3] * vec4[3]);
}

double Vector4::lengthSquared() const{
    return vec4[0] * vec4[0] + vec4[1] * vec4[1] + vec4[2] * vec4[2] + vec4[3] * vec4[3];
}

Vector4 Vector4::normalized() const{
    // Need some extra precision if the length is very small.
    double len = double(vec4[0]) * double(vec4[0]) +
                 double(vec4[1]) * double(vec4[1]) +
                 double(vec4[2]) * double(vec4[2]) +
                 double(vec4[3]) * double(vec4[3]);

    return *this / sqrt(len);

}

void Vector4::normalize(){
    // Need some extra precision if the length is very small.
    double len = double(vec4[0]) * double(vec4[0]) +
                 double(vec4[1]) * double(vec4[1]) +
                 double(vec4[2]) * double(vec4[2]) +
                 double(vec4[3]) * double(vec4[3]);

     len = sqrt(len);

     vec4[0] /= len;
     vec4[1] /= len;
     vec4[2] /= len;
     vec4[3] /= len;
}


const double& Vector4::operator[](int index) const{
    if(index < 0 && index > 3){
        //qDebug("const double& Vector4::operator[](int index) const has a wrong index: %i", index);
    }
    return vec4[index];
}

double& Vector4::operator[](int index){
    if(index < 0 && index > 3){
        //qDebug("double& Vector4::operator[](int index) has a wrong index: %i", index);
    }
    //var can be modified after return...
    return vec4[index];
}


const double& Vector4::operator()(int index) const{
    if(index < 0 && index > 3){
        //qDebug("const double& Vector4::operator()(int index) const has a wrong index: %i", index);
    }
    return vec4[index];
}

double& Vector4::operator()(int index){
    if(index < 0 && index > 3){
        //qDebug("double& Vector4::operator()(int index) has a wrong index: %i", index);
    }
    //var can be modified after return...
    return vec4[index];
}


Vector4& Vector4::operator+=(const Vector4 &vector){
    vec4[0] += vector.vec4[0];
    vec4[1] += vector.vec4[1];
    vec4[2] += vector.vec4[2];
    vec4[3] += vector.vec4[3];
    return *this;
}

Vector4& Vector4::operator-=(const Vector4 &vector){
    vec4[0] -= vector.vec4[0];
    vec4[1] -= vector.vec4[1];
    vec4[2] -= vector.vec4[2];
    vec4[3] -= vector.vec4[3];
    return *this;
}

Vector4& Vector4::operator*=(double factor){
    vec4[0] *= factor;
    vec4[1] *= factor;
    vec4[2] *= factor;
    vec4[3] *= factor;
    return *this;
}

Vector4& Vector4::operator*=(const Vector4 &vector){
    vec4[0] *= vector.vec4[0];
    vec4[1] *= vector.vec4[1];
    vec4[2] *= vector.vec4[2];
    vec4[3] *= vector.vec4[3];
    return *this;
}

Vector4& Vector4::operator/=(double divisor){
    vec4[0] /= divisor;
    vec4[1] /= divisor;
    vec4[2] /= divisor;
    vec4[3] /= divisor;
    return *this;
}


double Vector4::dotProduct(const Vector4& v1, const Vector4& v2){
    return v1.vec4[0] * v2.vec4[0] + v1.vec4[1] * v2.vec4[1] +
           v1.vec4[2] * v2.vec4[2] + v1.vec4[3] * v2.vec4[3];
}


//friend
bool operator==(const Vector4 &v1, const Vector4 &v2){
    return v1.vec4[0] == v2.vec4[0] && v1.vec4[1] == v2.vec4[1] &&
           v1.vec4[2] == v2.vec4[2] && v1.vec4[3] == v2.vec4[3];
}

bool operator!=(const Vector4 &v1, const Vector4 &v2){
    return v1.vec4[0] != v2.vec4[0] || v1.vec4[1] != v2.vec4[1] ||
           v1.vec4[2] != v2.vec4[2] || v1.vec4[3] != v2.vec4[3];
}

const Vector4 operator+(const Vector4 &v1, const Vector4 &v2){
    return Vector4(v1.vec4[0] + v2.vec4[0] , v1.vec4[1] + v2.vec4[1] ,
                   v1.vec4[2] + v2.vec4[2] , v1.vec4[3] + v2.vec4[3]);
}

const Vector4 operator-(const Vector4 &v1, const Vector4 &v2){
    return Vector4(v1.vec4[0] - v2.vec4[0] , v1.vec4[1] - v2.vec4[1] ,
                   v1.vec4[2] - v2.vec4[2] , v1.vec4[3] - v2.vec4[3]);
}

const Vector4 operator*(double factor, const Vector4 &vector){
    return Vector4(vector.vec4[0] * factor , vector.vec4[1] * factor ,
                   vector.vec4[2] * factor , vector.vec4[3] * factor);
}

const Vector4 operator*(const Vector4 &vector, double factor){
    return Vector4(vector.vec4[0] * factor , vector.vec4[1] * factor ,
                   vector.vec4[2] * factor , vector.vec4[3] * factor);
}

const Vector4 operator*(const Vector4 &v1, const Vector4& v2){
    return Vector4(v1.vec4[0] * v2.vec4[0] , v1.vec4[1] * v2.vec4[1] ,
                   v1.vec4[2] * v2.vec4[2] , v1.vec4[3] * v2.vec4[3]);
}

const Vector4 operator-(const Vector4 &vector){
    return Vector4(-vector.vec4[0], -vector.vec4[1],
                   -vector.vec4[2], -vector.vec4[3]);
}

const Vector4 operator/(const Vector4 &vector, double divisor){
    return Vector4(vector.vec4[0] / divisor , vector.vec4[1] / divisor ,
                   vector.vec4[2] / divisor , vector.vec4[3] / divisor);
}
