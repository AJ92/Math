#ifndef VECTOR3_H
#define VECTOR3_H

#include <QDebug>
#include "Math/Math.h"

/*

  Indizes
  |     0       1       2        |

  Pos
  |     x       y       z        |

*/
class Matrix4x4;

class Vector3
{
public:
    Vector3();
    Vector3(float f);
    Vector3(float x, float y, float z);
    Vector3(const float *vec3);
    Vector3(Vector3 const& vec);

    void set_to_null();

    bool is_null() const;

    void set_value(int index, float value);
    float get_value(int index);

    float x() const;
    float y() const;
    float z() const;


    void set_x(float x);
    void set_y(float y);
    void set_z(float z);




    float length() const;
    float lengthSquared() const;

    Vector3 normalized() const;
    void normalize();





    //operators...
    const float& operator[](int index) const;
    float& operator[](int index);

    const float& operator()(int index) const;
    float& operator()(int index);


    Vector3 &operator+=(const Vector3 &vector);
    Vector3 &operator-=(const Vector3 &vector);
    Vector3 &operator*=(float factor);
    Vector3 &operator*=(const Vector3 &vector);
    Vector3 &operator/=(float divisor);

    static float dotProduct(const Vector3& v1, const Vector3& v2);

    friend bool operator==(const Vector3 &v1, const Vector3 &v2);
    friend bool operator!=(const Vector3 &v1, const Vector3 &v2);
    friend const Vector3 operator+(const Vector3 &v1, const Vector3 &v2);
    friend const Vector3 operator-(const Vector3 &v1, const Vector3 &v2);
    friend const Vector3 operator*(float factor, const Vector3 &vector);
    friend const Vector3 operator*(const Vector3 &vector, float factor);
    friend const Vector3 operator*(const Vector3 &v1, const Vector3& v2);
    friend const Vector3 operator-(const Vector3 &vector);
    friend const Vector3 operator/(const Vector3 &vector, float divisor);


    //Matrix stuff
    friend Vector3 operator*(const Vector3& vector, const Matrix4x4& matrix);
    friend Vector3 operator*(const Matrix4x4& matrix, const Vector3& vector);


private:
    float vec3[3];

    //dirty inline hack ... i don't like it but its short :D
    Vector3(int) { }
};

#endif // VECTOR3_H
