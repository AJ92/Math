#ifndef MATRIX3X3_H
#define MATRIX3X3_H

#include <QDebug>

#include "Math/Vector/vector3.h"
#include "Math/Vector/vector4.h"
#include "Math/Matrix/Matrix4x4.h"
#include "Math/mathematics.h"

/*column-major notation for openGL use...
  but using in code row-major notation for indizes...

  following stuff is colomn-major notation


    Indizes
    |    0        3        6    |
    |                           |
    |    1        4        7    |
    |                           |
    |    2        5        8    |

*/


class Matrix3x3
{
public:
    //Multi-purpose constructors:
    Matrix3x3();
    Matrix3x3(float f1, float f2, float f3,
              float f4, float f5, float f6,
              float f7, float f8, float f9);
    Matrix3x3(const float *mat3);
    Matrix3x3(Matrix3x3 const& mat);

    void set_to_identity();
    bool is_identity();



    void set_value(int index, float value);
    float get_value(int index);

    float* get_array();
    void set_array(float mat3[]);

    void debug();

    //operators...
    const float& operator[](int index) const;
    float& operator[](int index);

    const float& operator()(int index) const;
    float& operator()(int index);

    const float& operator()(int row, int column) const;
    float& operator()(int row, int column);

    Matrix3x3& operator+=(const Matrix3x3& other);
    friend Matrix3x3 operator+(const Matrix3x3& m1, const Matrix3x3& m2);

    Matrix3x3& operator+=(const float& value);
    friend Matrix3x3 operator+(const Matrix3x3& m1, const float& value);
    friend Matrix3x3 operator+(const float& value, const Matrix3x3& m1);

    Matrix3x3& operator-=(const Matrix3x3& other);
    friend Matrix3x3 operator-(const Matrix3x3& m1, const Matrix3x3& m2);

    Matrix3x3& operator*=(const Matrix3x3& other);
    friend Matrix3x3 operator*(const Matrix3x3& m1, const Matrix3x3& m2);

    Matrix3x3& operator*=(const float& multiplier);
    friend Matrix3x3 operator*(const Matrix3x3& m1, const float& multiplier);
    friend Matrix3x3 operator*(const float& multiplier, const Matrix3x3& m1);

    Matrix3x3& operator/=(const float& divisor);
    friend Matrix3x3 operator/(const Matrix3x3& m1, const float& divisor);

private:
    float mat3[9];

    int flagBits;           // Flag bits from the enum below.

    enum {
        Identity        = 0x0001,   // Identity matrix
        General         = 0x0002,   // General matrix, unknown contents
        Translation     = 0x0004,   // Contains a simple translation
        Scale           = 0x0008,   // Contains a simple scale
        Rotation        = 0x0010    // Contains a simple rotation
    };

    Matrix3x3(int) { flagBits = General; }
};

#endif // MATRIX3X3_H
