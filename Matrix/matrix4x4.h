#ifndef MATRIX4X4_H
#define MATRIX4X4_H

#include <QDebug>

#include "Math/Math.h"



class Matrix3x3;
class Vector4;
class Vector3;

/*column-major notation for openGL use...
  but using in code row-major notation for indizes...

  following stuff is colomn-major notation


    Indizes
    |    0        4        8        12   |
    |                                    |
    |    1        5        9        13   |
    |                                    |
    |    2        6        10       14   |
    |                                    |
    |    3        7        11       15   |

    Translation
    |    1        0        0        X    |
    |                                    |
    |    0        1        0        Y    |
    |                                    |
    |    0        0        1        Z    |
    |                                    |
    |    0        0        0        1    |

    Scale
    |    SX       0        0        0    |
    |                                    |
    |    0        SY       0        0    |
    |                                    |
    |    0        0        SZ       0    |
    |                                    |
    |    0        0        0        1    |

    Rotation X
    |    1        0        0        0    |
    |                                    |
    |    0      cos(?)   sin(?)     0    |
    |                                    |
    |    0     -sin(?)   cos(?)     0    |
    |                                    |
    |    0        0        0        1    |

    Rotation Y
    |  cos(?)     0    -sin(?)      0    |
    |                                    |
    |    0        1        0        0    |
    |                                    |
    |  sin(?)     0     cos(?)      0    |
    |                                    |
    |    0        0        0        1    |

    Rotation Z
    |  cos(?)  -sin(?)     0        0    |
    |                                    |
    |  sin(?)   cos(?)     0        0    |
    |                                    |
    |    0        0        1        0    |
    |                                    |
    |    0        0        0        1    |

    Directions
    |   RX        RY       RZ       12   |  R = right
    |                                    |
    |   UX        UY       UZ       13   |  U = up
    |                                    |
    |   LX        LY       LZ       14   |  L = look at (front)
    |                                    |
    |    3        7        11       15   |

*/

class Matrix4x4
{
public:

    //Multi-purpose constructors:
    Matrix4x4();
    Matrix4x4(float f1, float f2, float f3, float f4,
              float f5, float f6, float f7, float f8,
              float f9, float f10, float f11, float f12,
              float f13, float f14, float f15, float f16);
    Matrix4x4(const float *mat4);
    Matrix4x4(Matrix4x4 const& mat);



    void set_to_identity();
    bool is_identity();

    void translate(float x, float y, float z);
    void scale(float sx, float sy, float sz);

    void rotate_x(float degrees);
    void rotate_y(float degrees);
    void rotate_z(float degrees);

    //Quaternion Class is missing...
    //void rotate(const Quaternion& quaternion);


    void pre_multiply(Matrix4x4 mat);
    void post_multiply(Matrix4x4 mat);



    float determinant() const;
    Matrix4x4 inverted(bool *invertible = 0) const;
    Matrix4x4 transposed() const;

    Matrix3x3 normalMatrix() const;

    void ortho(float left, float right, float bottom, float top, float nearPlane, float farPlane);
    void frustum(float left, float right, float bottom, float top, float nearPlane, float farPlane);
    void perspective(float angle, float aspect, float nearPlane, float farPlane);

    //Vector classes not done yet...
    //void lookAt(const Vector3& eye, const Vector3& center, const Vector3& up);

    Matrix3x3 rotationMatrix() const;

    void set_value(int index, float value);
    float get_value(int index);

    float* get_array();
    void set_array(float mat4[]);


    //Vector ops...
    Vector3 get_vector_up();
    Vector3 get_vector_right();
    Vector3 get_vector_look_at();

    Vector3 get_vector_pos();

    Vector3 get_vector_scale();


    void debug();

    //operators...
    const float& operator[](int index) const;
    float& operator[](int index);

    const float& operator()(int index) const;
    float& operator()(int index);

    const float& operator()(int row, int column) const;
    float& operator()(int row, int column);

    Matrix4x4& operator+=(const Matrix4x4& other);
    friend Matrix4x4 operator+(const Matrix4x4& m1, const Matrix4x4& m2);

    Matrix4x4& operator+=(const float& value);
    friend Matrix4x4 operator+(const Matrix4x4& m1, const float& value);
    friend Matrix4x4 operator+(const float& value, const Matrix4x4& m1);

    Matrix4x4& operator-=(const Matrix4x4& other);
    friend Matrix4x4 operator-(const Matrix4x4& m1, const Matrix4x4& m2);

    Matrix4x4& operator*=(const Matrix4x4& other);
    friend Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2);

    Matrix4x4& operator*=(const float& multiplier);
    friend Matrix4x4 operator*(const Matrix4x4& m1, const float& multiplier);
    friend Matrix4x4 operator*(const float& multiplier, const Matrix4x4& m1);

    Matrix4x4& operator/=(const float& divisor);
    friend Matrix4x4 operator/(const Matrix4x4& m1, const float& divisor);


    //Vector stuff
    friend Vector4 operator*(const Vector4& vector, const Matrix4x4& matrix);
    friend Vector4 operator*(const Matrix4x4& matrix, const Vector4& vector);

    friend Vector3 operator*(const Vector3& vector, const Matrix4x4& matrix);
    friend Vector3 operator*(const Matrix4x4& matrix, const Vector3& vector);

private:
    float mat4[16];

    int flagBits;           // Flag bits from the enum below.

    enum {
        Identity        = 0x0001,   // Identity matrix
        General         = 0x0002,   // General matrix, unknown contents
        Translation     = 0x0004,   // Contains a simple translation
        Scale           = 0x0008,   // Contains a simple scale
        Rotation        = 0x0010    // Contains a simple rotation
    };

    Matrix4x4 orthonormalInverse() const;

    //dirty inline hack ... i don't like it but its short :D
    Matrix4x4(int) { flagBits = General; }

};

#endif // MATRIX4X4_H
