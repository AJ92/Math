#ifndef VECTOR4_H
#define VECTOR4_H

#include <QDebug>
#include "Math/Math.h"

/*

  Indizes
  |     0       1       2       3       |

  Pos
  |     x       y       z       w       |

*/
class Matrix4x4;

class Vector4
{
public:
    Vector4();
    Vector4(float f);
    Vector4(float x, float y, float z, float w);
    Vector4(const float *vec4);
    Vector4(Vector4 const& vec);

    void set_to_null();

    bool is_null() const;

    void set_value(int index, float value);
    float get_value(int index);

    float x() const;
    float y() const;
    float z() const;
    float w() const;

    void set_x(float x);
    void set_y(float y);
    void set_z(float z);
    void set_w(float w);



    float length() const;
    float lengthSquared() const;

    Vector4 normalized() const;
    void normalize();





    //operators...
    const float& operator[](int index) const;
    float& operator[](int index);

    const float& operator()(int index) const;
    float& operator()(int index);


    Vector4 &operator+=(const Vector4 &vector);
    Vector4 &operator-=(const Vector4 &vector);
    Vector4 &operator*=(float factor);
    Vector4 &operator*=(const Vector4 &vector);
    Vector4 &operator/=(float divisor);

    static float dotProduct(const Vector4& v1, const Vector4& v2);

    friend bool operator==(const Vector4 &v1, const Vector4 &v2);
    friend bool operator!=(const Vector4 &v1, const Vector4 &v2);
    friend const Vector4 operator+(const Vector4 &v1, const Vector4 &v2);
    friend const Vector4 operator-(const Vector4 &v1, const Vector4 &v2);
    friend const Vector4 operator*(float factor, const Vector4 &vector);
    friend const Vector4 operator*(const Vector4 &vector, float factor);
    friend const Vector4 operator*(const Vector4 &v1, const Vector4& v2);
    friend const Vector4 operator-(const Vector4 &vector);
    friend const Vector4 operator/(const Vector4 &vector, float divisor);


    //Matrix stuff
    friend Vector4 operator*(const Vector4& vector, const Matrix4x4& matrix);
    friend Vector4 operator*(const Matrix4x4& matrix, const Vector4& vector);


private:
    float vec4[4];

    //dirty inline hack ... i don't like it but its short :D
    Vector4(int) { }

};

#endif // VECTOR4_H
