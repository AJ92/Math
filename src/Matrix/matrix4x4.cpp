#include "Matrix/matrix4x4.h"

/*!
    Constructs the 4x4 Matrix and sets it to it's Identity.

    \sa set_to_identity()
*/
Matrix4x4::Matrix4x4()
{
    set_to_identity();
}

/*!
    \overload

    Constructs the 4x4 Matrix and copies the
    components from \a mat in it's own 4x4 Matrix.

    \sa set_array()

*/
Matrix4x4::Matrix4x4(const Matrix4x4 &mat){
    this->mat4[0] = mat.mat4[0];
    this->mat4[1] = mat.mat4[1];
    this->mat4[2] = mat.mat4[2];
    this->mat4[3] = mat.mat4[3];
    this->mat4[4] = mat.mat4[4];
    this->mat4[5] = mat.mat4[5];
    this->mat4[6] = mat.mat4[6];
    this->mat4[7] = mat.mat4[7];
    this->mat4[8] = mat.mat4[8];
    this->mat4[9] = mat.mat4[9];
    this->mat4[10] = mat.mat4[10];
    this->mat4[11] = mat.mat4[11];
    this->mat4[12] = mat.mat4[12];
    this->mat4[13] = mat.mat4[13];
    this->mat4[14] = mat.mat4[14];
    this->mat4[15] = mat.mat4[15];
    flagBits = General;
}

/*!
    \overload

    Constructs the 4x4 Matrix and sets the components to \a f1 - \a f16.
*/
Matrix4x4::Matrix4x4(double f1, double f2, double f3, double f4,
                    double f5, double f6, double f7, double f8,
                    double f9, double f10, double f11, double f12,
                    double f13, double f14, double f15, double f16)
{
    mat4[0] = f1;
    mat4[1] = f2;
    mat4[2] = f3;
    mat4[3] = f4;
    mat4[4] = f5;
    mat4[5] = f6;
    mat4[6] = f7;
    mat4[7] = f8;
    mat4[8] = f9;
    mat4[9] = f10;
    mat4[10] = f11;
    mat4[11] = f12;
    mat4[12] = f13;
    mat4[13] = f14;
    mat4[14] = f15;
    mat4[15] = f16;
    flagBits = General;
}

/*!
    \overload

    Constructs the 4x4 Matrix and sets the
    components from \a mat in it's own 4x4 Matrix.

    \sa set_array()

*/
Matrix4x4::Matrix4x4(const double *mat)
{
    for (int i = 0; i < 16; i++){
        mat4[i] = mat[i];
    }
    flagBits = General;
}

/*!

    Sets the matrix to it's Identity.

*/
void Matrix4x4::set_to_identity(){
    //setting matrix to it's Identity
    mat4[0]    = mat4[5]  = mat4[10] = mat4[15] = 1.0;
    mat4[1]    = mat4[2]  = mat4[3]  = mat4[4]  = 0.0;
    mat4[6]    = mat4[7]  = mat4[8]  = mat4[9]  = 0.0;
    mat4[11]   = mat4[12] = mat4[13] = mat4[14] = 0.0;
    flagBits = Identity;
}

/*!
    Checks if the matrix is an identity matrix.
  */
bool Matrix4x4::is_identity(){
    if (flagBits == Identity)
        return true;
    if (mat4[0] != 1.0f || mat4[1] != 0.0f || mat4[2] != 0.0f)
        return false;
    if (mat4[3] != 0.0f || mat4[4] != 0.0f || mat4[5] != 1.0f)
        return false;
    if (mat4[6] != 0.0f || mat4[7] != 0.0f || mat4[8] != 0.0f)
        return false;
    if (mat4[9] != 0.0f || mat4[10] != 1.0f || mat4[11] != 0.0f)
        return false;
    if (mat4[12] != 0.0f || mat4[13] != 0.0f || mat4[14] != 0.0f)
        return false;
    return (mat4[15] == 1.0f);
}

/*!
    Loads the Identity and translates it to \a x, \a y and \a z coordinates.
*/
void Matrix4x4::translate(double x, double y, double z){
    set_to_identity();
    // Translate slots.
    mat4[12] = x;
    mat4[13] = y;
    mat4[14] = z;
    flagBits = Translation;
}

/*!
    Loads the Identity and scales it to \a sx, \a sy and \a sz.
*/
void Matrix4x4::scale(double sx, double sy, double sz){
    set_to_identity();
    // Scale slots.
    mat4[0]   = sx;
    mat4[5]   = sy;
    mat4[10]  = sz;
    flagBits = Scale;
}


void Matrix4x4::translate(Vector3 translate_vector){
    translate(translate_vector[0],translate_vector[1],translate_vector[2]);
}

void Matrix4x4::scale(Vector3 scale_vector){
    scale(scale_vector[0],scale_vector[1],scale_vector[2]);
}




//
//
//      Might need to combine all rotations in one function ...
//
//      and make all transformations be applied to the current matrix...
//
//


/*!
    Loads the Identity and rotates it around it's X-axis by \a degrees.
*/
void Matrix4x4::rotate_x(double degrees){
    set_to_identity();

    double radians = degreesToRadians(degrees);

    // Rotate X formula.
    mat4[5] = cos(radians);
    mat4[6] = -sin(radians);
    mat4[9] = -mat4[6];
    mat4[10] = mat4[5];
    flagBits = Rotation;
}

/*!
    Loads the Identity and rotates it around it's Y-axis by \a degrees.
*/
void Matrix4x4::rotate_y(double degrees){
    set_to_identity();

    double radians = degreesToRadians(degrees);

    // Rotate Y formula.
    mat4[0] = cos(radians);
    mat4[2] = sin(radians);
    mat4[8] = -mat4[2];
    mat4[10] = mat4[0];
    flagBits = Rotation;
}

/*!
    Loads the Identity and rotates it around it's Z-axis by \a degrees.
*/
void Matrix4x4::rotate_z(double degrees){
    set_to_identity();

    double radians = degreesToRadians(degrees);

    // Rotate Z formula.
    mat4[0] = cos(radians);
    mat4[1] = sin(radians);
    mat4[4] = -mat4[1];
    mat4[5] = mat4[0];
    flagBits = Rotation;
}




void Matrix4x4::rotate(double angle, const Vector3& vector)
{
    rotate(angle, vector.x(), vector.y(), vector.z());
}

void Matrix4x4::rotate(double angle, double x, double y, double z)
{
    if (angle == 0.0f)
        return;
    Matrix4x4 m(1); // The "1" says to not load the identity.
    double c, s, ic;
    if (angle == 90.0 || angle == -270.0) {
        s = 1.0;
        c = 0.0;
    } else if (angle == -90.0 || angle == 270.0) {
        s = -1.0;
        c = 0.0;
    } else if (angle == 180.0 || angle == -180.0) {
        s = 0.0;
        c = -1.0;
    } else {
        double a = angle * kPI180;
        c = cos(a);
        s = sin(a);
    }
    bool quick = false;
    if (x == 0.0) {
        if (y == 0.0) {
            if (z != 0.0) {
                // Rotate around the Z axis.
                m.set_to_identity();
                m.mat4[0] = c;
                m.mat4[5] = c;
                if (z < 0.0) {
                    m.mat4[4] = s;
                    m.mat4[1] = -s;
                } else {
                    m.mat4[4] = -s;
                    m.mat4[1] = s;
                }
                m.flagBits = General;
                quick = true;
            }
        } else if (z == 0.0) {
            // Rotate around the Y axis.
            m.set_to_identity();
            m.mat4[0] = c;
            m.mat4[10] = c;
            if (y < 0.0) {
                m.mat4[8] = -s;
                m.mat4[2] = s;
            } else {
                m.mat4[8] = s;
                m.mat4[2] = -s;
            }
            m.flagBits = General;
            quick = true;
        }
    } else if (y == 0.0 && z == 0.0) {
        // Rotate around the X axis.
        m.set_to_identity();
        m.mat4[5] = c;
        m.mat4[10] = c;
        if (x < 0.0) {
            m.mat4[9] = s;
            m.mat4[6] = -s;
        } else {
            m.mat4[9] = -s;
            m.mat4[6] = s;
        }
        m.flagBits = General;
        quick = true;
    }
    if (!quick) {
        double len = x * x + y * y + z * z;
        if (!Helper::fuzzyIsNull(len - 1.0) && !Helper::fuzzyIsNull(len)) {
            len = sqrt(len);
            x /= len;
            y /= len;
            z /= len;
        }
        ic = 1.0 - c;
        m.mat4[0] = x * x * ic + c;
        m.mat4[4] = x * y * ic - z * s;
        m.mat4[8] = x * z * ic + y * s;
        m.mat4[12] = 0.0;
        m.mat4[1] = y * x * ic + z * s;
        m.mat4[5] = y * y * ic + c;
        m.mat4[9] = y * z * ic - x * s;
        m.mat4[13] = 0.0;
        m.mat4[2] = x * z * ic - y * s;
        m.mat4[6] = y * z * ic + x * s;
        m.mat4[10] = z * z * ic + c;
        m.mat4[14] = 0.0;
        m.mat4[3] = 0.0;
        m.mat4[7] = 0.0;
        m.mat4[11] = 0.0;
        m.mat4[15] = 1.0;
    }
    int flags = flagBits;
    *this *= m;
    if (flags != Identity)
        flagBits = flags | Rotation;
    else
        flagBits = Rotation;
}


/*!
    Multiplyies the matrix with \a mat in prefix order.
*/
void Matrix4x4::pre_multiply(Matrix4x4 mat){
    *this = mat * *this;
}

/*!
    Multiplyies the matrix with \a mat in postfix order.
*/
void Matrix4x4::post_multiply(Matrix4x4 mat){
    *this = *this * mat;
}



// Calculate the determinant of a 3x3 sub-matrix.
//     | A B C |
// M = | D E F |   det(M) = A * (EI - HF) - B * (DI - GF) + C * (DH - GE)
//     | G H I |
static inline double matrixDet3 (const double mat4[16], int row0, int row1, int row2,
     int col0, int col1, int col2)
{
    return mat4[4*row0+col0] *
                (mat4[4*row1+col1] * mat4[4*row2+col2] -
                 mat4[4*row1+col2] * mat4[4*row2+col1]) -
           mat4[4*row1+col0] *
                (mat4[4*row0+col1] * mat4[4*row2+col2] -
                 mat4[4*row0+col2] * mat4[4*row2+col1]) +
           mat4[4*row2+col0] *
                (mat4[4*row0+col1] * mat4[4*row1+col2] -
                 mat4[4*row0+col2] * mat4[4*row1+col1]);
}

// Calculate the determinant of a 4x4 matrix.
static inline double matrixDet4(const double mat4[16])
{
    double det;
    det  = mat4[0] * matrixDet3(mat4, 1, 2, 3, 1, 2, 3);
    det -= mat4[4] * matrixDet3(mat4, 0, 2, 3, 1, 2, 3);
    det += mat4[8] * matrixDet3(mat4, 0, 1, 3, 1, 2, 3);
    det -= mat4[12] * matrixDet3(mat4, 0, 1, 2, 1, 2, 3);
    return det;
}

/*!
    Calculates the determinant.
    The 4x4 matrix inverse algorithm is based on that described at:
    http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q24
    Some optimization has been done to avoid making copies of 3x3
    sub-matrices and to unroll the loops.
*/
double Matrix4x4::determinant() const{
    return double(matrixDet4(mat4));
}

// Helper routine for inverting orthonormal matrices that consist
// of just rotations and translations.
Matrix4x4 Matrix4x4::orthonormalInverse() const
{
    Matrix4x4 result(1);  // The '1' says not to load identity

    result.mat4[0] = mat4[0];
    result.mat4[4] = mat4[1];
    result.mat4[8] = mat4[2];

    result.mat4[1] = mat4[4];
    result.mat4[5] = mat4[5];
    result.mat4[9] = mat4[6];

    result.mat4[2] = mat4[8];
    result.mat4[6] = mat4[9];
    result.mat4[10] = mat4[10];

    result.mat4[3] = 0.0f;
    result.mat4[7] = 0.0f;
    result.mat4[11] = 0.0f;

    //might have wrong values
    result.mat4[12] = -(result.mat4[0] * mat4[12] + result.mat4[4] * mat4[13] + result.mat4[8] * mat4[14]);
    result.mat4[13] = -(result.mat4[1] * mat4[12] + result.mat4[5] * mat4[13] + result.mat4[9] * mat4[14]);
    result.mat4[14] = -(result.mat4[2] * mat4[12] + result.mat4[6] * mat4[13] + result.mat4[10] * mat4[14]);
    result.mat4[15] = 1.0f;

    return result;
}

/*!
    Returns the inverse of this matrix.  Returns the identity if
    this matrix cannot be inverted; i.e. determinant() is zero.
    If \a invertible is not null, then true will be written to
    that location if the matrix can be inverted; false otherwise.

    If the matrix is recognized as the identity or an orthonormal
    matrix, then this function will quickly invert the matrix
    using optimized routines.

    \sa determinant(), normalMatrix()
*/
Matrix4x4 Matrix4x4::inverted(bool *invertible) const{

    // Handle some of the easy cases first.
    if (flagBits == Identity) {
        if (invertible)
            *invertible = true;
        return Matrix4x4();
    } else if (flagBits == Translation) {
        Matrix4x4 inv;
        inv.mat4[12] = -mat4[12];
        inv.mat4[13] = -mat4[13];
        inv.mat4[14] = -mat4[14];
        inv.flagBits = Translation;
        if (invertible)
            *invertible = true;
        return inv;
    } else if (flagBits == Rotation || flagBits == (Rotation | Translation)) {
        if (invertible)
            *invertible = true;
        return orthonormalInverse();
    }

    Matrix4x4 inv(1); // The "1" says to not load the identity.

    float det = matrixDet4(mat4);
    if (det == 0.0f) {
        if (invertible)
            *invertible = false;
        return Matrix4x4();
    }
    det = 1.0f / det;

    inv.mat4[0] =  matrixDet3(mat4, 1, 2, 3, 1, 2, 3) * det;
    inv.mat4[1] = -matrixDet3(mat4, 0, 2, 3, 1, 2, 3) * det;
    inv.mat4[2] =  matrixDet3(mat4, 0, 1, 3, 1, 2, 3) * det;
    inv.mat4[3] = -matrixDet3(mat4, 0, 1, 2, 1, 2, 3) * det;
    inv.mat4[4] = -matrixDet3(mat4, 1, 2, 3, 0, 2, 3) * det;
    inv.mat4[5] =  matrixDet3(mat4, 0, 2, 3, 0, 2, 3) * det;
    inv.mat4[6] = -matrixDet3(mat4, 0, 1, 3, 0, 2, 3) * det;
    inv.mat4[7] =  matrixDet3(mat4, 0, 1, 2, 0, 2, 3) * det;
    inv.mat4[8] =  matrixDet3(mat4, 1, 2, 3, 0, 1, 3) * det;
    inv.mat4[9] = -matrixDet3(mat4, 0, 2, 3, 0, 1, 3) * det;
    inv.mat4[10] =  matrixDet3(mat4, 0, 1, 3, 0, 1, 3) * det;
    inv.mat4[11] = -matrixDet3(mat4, 0, 1, 2, 0, 1, 3) * det;
    inv.mat4[12] = -matrixDet3(mat4, 1, 2, 3, 0, 1, 2) * det;
    inv.mat4[13] =  matrixDet3(mat4, 0, 2, 3, 0, 1, 2) * det;
    inv.mat4[14] = -matrixDet3(mat4, 0, 1, 3, 0, 1, 2) * det;
    inv.mat4[15] =  matrixDet3(mat4, 0, 1, 2, 0, 1, 2) * det;

    if (invertible)
        *invertible = true;
    return inv;
}

/*!

*/
Matrix4x4 Matrix4x4::transposed() const{
    Matrix4x4 result(1); // The "1" says to not load the identity.
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            result.mat4[col+row*4] = mat4[4*row+col];
        }
    }
    return result;
}

/*!

*/
Matrix3x3 Matrix4x4::normalMatrix() const{

    Matrix3x3 inv;

    // Handle the simple cases first.
    if (flagBits == Identity || flagBits == Translation) {
        return inv;
    } else if (flagBits == Scale || flagBits == (Translation | Scale)) {
        if (mat4[0] == 0.0 || mat4[5] == 0.0 || mat4[10] == 0.0)
            return inv;
        inv[0] = 1.0 / mat4[0];
        inv[4] = 1.0 / mat4[5];
        inv[8] = 1.0 / mat4[10];
        return inv;
    }

    double det = matrixDet3(mat4, 0, 1, 2, 0, 1, 2);
    if (det == 0.0)
        return inv;
    det = 1.0 / det;

    double *invm = inv.get_array();

    // Invert and transpose in a single step.
    invm[0 + 0 * 3] =  (mat4[5] * mat4[10] - mat4[9] * mat4[6]) * det;
    invm[1 + 0 * 3] = -(mat4[4] * mat4[10] - mat4[6] * mat4[8]) * det;
    invm[2 + 0 * 3] =  (mat4[4] * mat4[9] - mat4[5] * mat4[8]) * det;
    invm[0 + 1 * 3] = -(mat4[1] * mat4[10] - mat4[9] * mat4[2]) * det;
    invm[1 + 1 * 3] =  (mat4[0] * mat4[10] - mat4[2] * mat4[8]) * det;
    invm[2 + 1 * 3] = -(mat4[0] * mat4[9] - mat4[1] * mat4[8]) * det;
    invm[0 + 2 * 3] =  (mat4[1] * mat4[6] - mat4[2] * mat4[5]) * det;
    invm[1 + 2 * 3] = -(mat4[0] * mat4[6] - mat4[2] * mat4[4]) * det;
    invm[2 + 2 * 3] =  (mat4[0] * mat4[5] - mat4[4] * mat4[1]) * det;

    return inv;
}



/*!
    Multiples this matrix by another that rotates coordinates according
    to a specified \a quaternion.  The \a quaternion is assumed to have
    been normalized.

    \sa scale(), translate(), QQuaternion
*/

/*
void Matrix4x4::rotate(const Quaternion& quaternion)
{
    // Algorithm from:
    // http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q54
    Matrix4x4 m(1);
    double xx = quaternion.x() * quaternion.x();
    double xy = quaternion.x() * quaternion.y();
    double xz = quaternion.x() * quaternion.z();
    double xw = quaternion.x() * quaternion.scalar();
    double yy = quaternion.y() * quaternion.y();
    double yz = quaternion.y() * quaternion.z();
    double yw = quaternion.y() * quaternion.scalar();
    double zz = quaternion.z() * quaternion.z();
    double zw = quaternion.z() * quaternion.scalar();
    m(0,0) = 1.0f - 2 * (yy + zz);
    m(1,0) =        2 * (xy - zw);
    m(2,0) =        2 * (xz + yw);
    m(3,0) = 0.0f;
    m(0,1) =        2 * (xy + zw);
    m(1,1) = 1.0f - 2 * (xx + zz);
    m(2,1) =        2 * (yz - xw);
    m(3,1) = 0.0f;
    m(0,2) =        2 * (xz - yw);
    m(1,2) =        2 * (yz + xw);
    m(2,2) = 1.0f - 2 * (xx + yy);
    m(3,2) = 0.0f;
    m(0,3) = 0.0f;
    m(1,3) = 0.0f;
    m(2,3) = 0.0f;
    m(3,3) = 1.0f;
    int flags = flagBits;
    *this *= m;
    if (flags != Identity)
        flagBits = flags | Rotation;
    else
        flagBits = Rotation;
}
*/



/*!
    Multiplies this matrix by another that applies an orthographic
    projection for a window with lower-left corner (\a left, \a bottom),
    upper-right corner (\a right, \a top), and the specified \a nearPlane
    and \a farPlane clipping planes.

    \sa frustum(), perspective()
*/
void Matrix4x4::ortho(double left, double right, double bottom, double top, double nearPlane, double farPlane){
    // Bail out if the projection volume is zero-sized.
    if (left == right || bottom == top || nearPlane == farPlane)
        return;

    // Construct the projection.
    double width = right - left;
    double invheight = top - bottom;
    double clip = farPlane - nearPlane;

    //Vector class is incomplete and thus functions are missing...
    //A more efficient way of ortho projection...

    /*
    if (clip == 2.0f && (nearPlane + farPlane) == 0.0f) {
        // We can express this projection as a translate and scale
        // which will be more efficient to modify with further
        // transformations than producing a "General" matrix.
        translate(Vector3
            (-(left + right) / width,
             -(top + bottom) / invheight,
             0.0f));
        scale(Vector3
            (2.0f / width,
             2.0f / invheight,
             -1.0f));
        return;
    }
    */

    Matrix4x4 m(1);
    m(0,0) = 2.0 / width;
    m(1,0) = 0.0;
    m(2,0) = 0.0;
    m(3,0) = -(left + right) / width;
    m(0,1) = 0.0;
    m(1,1) = 2.0 / invheight;
    m(2,1) = 0.0;
    m(3,1) = -(top + bottom) / invheight;
    m(0,2) = 0.0;
    m(1,2) = 0.0;
    m(2,2) = -2.0 / clip;
    m(3,2) = -(nearPlane + farPlane) / clip;
    m(0,3) = 0.0;
    m(1,3) = 0.0;
    m(2,3) = 0.0;
    m(3,3) = 1.0;

    // Apply the projection.
    *this *= m;
    return;
}

/*!
    Multiplies this matrix by another that applies a perspective
    frustum projection for a window with lower-left corner (\a left, \a bottom),
    upper-right corner (\a right, \a top), and the specified \a nearPlane
    and \a farPlane clipping planes.

    \sa ortho(), perspective()
*/
void Matrix4x4::frustum(double left, double right, double bottom, double top, double nearPlane, double farPlane){
    // Bail out if the projection volume is zero-sized.
    if (left == right || bottom == top || nearPlane == farPlane)
        return;

    // Construct the projection.
    Matrix4x4 m(1);
    double width = right - left;
    double invheight = top - bottom;
    double clip = farPlane - nearPlane;
    m(0,0) = 2.0 * nearPlane / width;
    m(1,0) = 0.0;
    m(2,0) = (left + right) / width;
    m(3,0) = 0.0;
    m(0,1) = 0.0;
    m(1,1) = 2.0 * nearPlane / invheight;
    m(2,1) = (top + bottom) / invheight;
    m(3,1) = 0.0;
    m(0,2) = 0.0;
    m(1,2) = 0.0;
    m(2,2) = -(nearPlane + farPlane) / clip;
    m(3,2) = -2.0 * nearPlane * farPlane / clip;
    m(0,3) = 0.0;
    m(1,3) = 0.0;
    m(2,3) = -1.0;
    m(3,3) = 0.0;

    // Apply the projection.
    *this *= m;
}

/*!
    Multiplies this matrix by another that applies a perspective
    projection.  The field of view will be \a angle degrees within
    a window with a given \a aspect ratio.  The projection will
    have the specified \a nearPlane and \a farPlane clipping planes.

    \sa ortho(), frustum()
*/
void Matrix4x4::perspective(double angle, double aspect, double nearPlane, double farPlane){
    // Bail out if the projection volume is zero-sized.
    if (nearPlane == farPlane || aspect == 0.0)
        return;

    // Construct the projection.
    Matrix4x4 m(1);
    double radians = (angle / 2.0) * kPI180;
    double sine = sin(radians);
    if (sine == 0.0)
        return;
    double cotan = cos(radians) / sine;
    double clip = farPlane - nearPlane;
    m(0,0) = cotan / aspect;
    m(1,0) = 0.0;
    m(2,0) = 0.0;
    m(3,0) = 0.0;
    m(0,1) = 0.0;
    m(1,1) = cotan;
    m(2,1) = 0.0;
    m(3,1) = 0.0;
    m(0,2) = 0.0;
    m(1,2) = 0.0;
    m(2,2) = -(nearPlane + farPlane) / clip;
    m(3,2) = -(2.0 * nearPlane * farPlane) / clip;
    m(0,3) = 0.0;
    m(1,3) = 0.0;
    m(2,3) = -1.0;
    m(3,3) = 0.0;

    // Apply the projection.
    *this *= m;
}


/*!
    Multiplies this matrix by another that applies an \a eye position
    transformation.  The \a center value indicates the center of the
    view that the \a eye is looking at.  The \a up value indicates
    which direction should be considered up with respect to the \a eye.
*/

/*
void Matrix4x4::lookAt(const Vector3& eye, const Vector3& center, const Vector3& up){
    QVector3D forward = (center - eye).normalized();
    QVector3D side = QVector3D::crossProduct(forward, up).normalized();
    QVector3D upVector = QVector3D::crossProduct(side, forward);

    QMatrix4x4 m(1);

    m.m[0][0] = side.x();
    m.m[1][0] = side.y();
    m.m[2][0] = side.z();
    m.m[3][0] = 0.0f;
    m.m[0][1] = upVector.x();
    m.m[1][1] = upVector.y();
    m.m[2][1] = upVector.z();
    m.m[3][1] = 0.0f;
    m.m[0][2] = -forward.x();
    m.m[1][2] = -forward.y();
    m.m[2][2] = -forward.z();
    m.m[3][2] = 0.0f;
    m.m[0][3] = 0.0f;
    m.m[1][3] = 0.0f;
    m.m[2][3] = 0.0f;
    m.m[3][3] = 1.0f;

    *this *= m;
    translate(-eye);
}
*/


Matrix3x3 Matrix4x4::rotationMatrix() const{
    return Matrix3x3(mat4[0],mat4[1],mat4[2],
                     mat4[4],mat4[5],mat4[6],
                     mat4[8],mat4[9],mat4[10]);
}


/*!
    Sets the given \a value at the \a index in the matrix.
*/
void Matrix4x4::set_value(int index, double value){
    if(index < 0 || index > 15){
        //qDebug("void Matrix4x4::set_value(int index, double value) has a wrong index: %i", index);
        return;
    }
    mat4[index] = value;
    flagBits = General;
}

/*!
    Returns the value at the given \a index.
*/
double Matrix4x4::get_value(int index){
    if(index < 0 || index > 15){
        //qDebug("double Matrix4x4::get_value(int index) has a wrong index: %i", index);
        return 0;
    }
    return mat4[index];
}

/*!
    Returns the matrix as an array of doubles.
*/
double* Matrix4x4::get_array(){
    return mat4;
}

/*!
    Sets the components of the matrix to the values stored int the \a mat4 array.
*/
void Matrix4x4::set_array(double mat4[]){
    this->mat4[0] = mat4[0];
    this->mat4[1] = mat4[1];
    this->mat4[2] = mat4[2];
    this->mat4[3] = mat4[3];
    this->mat4[4] = mat4[4];
    this->mat4[5] = mat4[5];
    this->mat4[6] = mat4[6];
    this->mat4[7] = mat4[7];
    this->mat4[8] = mat4[8];
    this->mat4[9] = mat4[9];
    this->mat4[10] = mat4[10];
    this->mat4[11] = mat4[11];
    this->mat4[12] = mat4[12];
    this->mat4[13] = mat4[13];
    this->mat4[14] = mat4[14];
    this->mat4[15] = mat4[15];
    flagBits = General;
}



//Vector ops:

Vector3 Matrix4x4::get_vector_up(){
    return Vector3(mat4[1],mat4[5],mat4[9]);
}

Vector3 Matrix4x4::get_vector_right(){
    return Vector3(mat4[0],mat4[4],mat4[8]);
}

Vector3 Matrix4x4::get_vector_look_at(){
    return Vector3(mat4[2],mat4[6],mat4[10]);
}

Vector3 Matrix4x4::get_vector_pos(){
    return Vector3(mat4[12],mat4[13],mat4[14]);
}

Vector3 Matrix4x4::get_vector_scale(){
    return Vector3(mat4[0],mat4[5],mat4[10]);
}





/*!
    Outputs the matrix in the console by using Qt's qDebug() function.
*/
void Matrix4x4::debug(){
    //qDebug("Matrix:");
    //qDebug("%f  %f  %f  %f",mat4[0],mat4[1],mat4[2],mat4[3]);
    //qDebug("%f  %f  %f  %f",mat4[4],mat4[5],mat4[6],mat4[7]);
    //qDebug("%f  %f  %f  %f",mat4[8],mat4[9],mat4[10],mat4[11]);
    //qDebug("%f  %f  %f  %f",mat4[12],mat4[13],mat4[14],mat4[15]);
}


//operators


/*!
    For easier access to the values of the matrix you can use this function.
*/
const double& Matrix4x4::operator[](int index) const
{
    if(index < 0 && index > 15){
        //qDebug("const double& Matrix4x4::operator()(int index) has a wrong index: %i", index);
    }
    return mat4[index];
}

/*!
    For easier access to the values of the matrix you can use this function.
*/
double& Matrix4x4::operator[](int index)
{
    if(index < 0 && index > 15){
        //qDebug("double& Matrix4x4::operator()(int index) has a wrong index: %i", index);
    }
    flagBits = General;
    return mat4[index];
}

/*!
    For easier access to the values of the matrix you can use this function.
*/
const double& Matrix4x4::operator()(int index) const
{
    if(index < 0 && index > 15){
        //qDebug("const double& Matrix4x4::operator()(int index) has a wrong index: %i", index);
    }
    return mat4[index];
}

/*!
    For easier access to the values of the matrix you can use this function.
*/
double& Matrix4x4::operator()(int index)
{
    if(index < 0 && index > 15){
        //qDebug("double& Matrix4x4::operator()(int index) has a wrong index: %i", index);
    }
    flagBits = General;
    return mat4[index];
}

/*!
    For easier access to the values of the matrix you can use this function.
*/
const double& Matrix4x4::operator()(int aRow, int aColumn) const
{
    if(aRow < 0 && aRow > 3 && aColumn < 0 && aColumn > 3){
        //qDebug("const double& Matrix4x4::operator()(int aRow, int aColumn) has wrong indizes: %i  %i", aRow, aColumn);
    }
    return mat4[aColumn+aRow*4];
}

/*!
    For easier access to the values of the matrix you can use this function.
*/
double& Matrix4x4::operator()(int aRow, int aColumn)
{
    if(aRow < 0 && aRow > 3 && aColumn < 0 && aColumn > 3){
        //qDebug("double& Matrix4x4::operator()(int aRow, int aColumn) has wrong indizes: %i  %i", aRow, aColumn);
    }
    flagBits = General;
    return mat4[aColumn+aRow*4];
}


/*!
    Adds the \a other matrix to itself and returns itself.
*/
Matrix4x4& Matrix4x4::operator+=(const Matrix4x4& other){
    *this = *this + other;
    return *this;
}

/*!
    Adds the \a m1 and \a m2 matrix and returns the result.
*/
Matrix4x4 operator+(const Matrix4x4& m1, const Matrix4x4& m2){
    Matrix4x4 m(1);
    m.mat4[0] = m1.mat4[0] + m2.mat4[0];
    m.mat4[1] = m1.mat4[1] + m2.mat4[1];
    m.mat4[2] = m1.mat4[2] + m2.mat4[2];
    m.mat4[3] = m1.mat4[3] + m2.mat4[3];
    m.mat4[4] = m1.mat4[4] + m2.mat4[4];
    m.mat4[5] = m1.mat4[5] + m2.mat4[5];
    m.mat4[6] = m1.mat4[6] + m2.mat4[6];
    m.mat4[7] = m1.mat4[7] + m2.mat4[7];
    m.mat4[8] = m1.mat4[8] + m2.mat4[8];
    m.mat4[9] = m1.mat4[9] + m2.mat4[9];
    m.mat4[10] = m1.mat4[10] + m2.mat4[10];
    m.mat4[11] = m1.mat4[11] + m2.mat4[11];
    m.mat4[12] = m1.mat4[12] + m2.mat4[12];
    m.mat4[13] = m1.mat4[13] + m2.mat4[13];
    m.mat4[14] = m1.mat4[14] + m2.mat4[14];
    m.mat4[15] = m1.mat4[15] + m2.mat4[15];
    //m is General matrix...
    return m;
}



/*!
    Adds the \a value to itself and returns itself.
*/
Matrix4x4& Matrix4x4::operator+=(const double& value){
    *this = *this + value;
    return *this;
}

/*!
    Adds the \a value to \a m1 matrix and returns the result.
*/
Matrix4x4 operator+(const Matrix4x4& m1, const double& value){
    Matrix4x4 m(1);
    m.mat4[0] = m1.mat4[0] + value;
    m.mat4[1] = m1.mat4[1] + value;
    m.mat4[2] = m1.mat4[2] + value;
    m.mat4[3] = m1.mat4[3] + value;
    m.mat4[4] = m1.mat4[4] + value;
    m.mat4[5] = m1.mat4[5] + value;
    m.mat4[6] = m1.mat4[6] + value;
    m.mat4[7] = m1.mat4[7] + value;
    m.mat4[8] = m1.mat4[8] + value;
    m.mat4[9] = m1.mat4[9] + value;
    m.mat4[10] = m1.mat4[10] + value;
    m.mat4[11] = m1.mat4[11] + value;
    m.mat4[12] = m1.mat4[12] + value;
    m.mat4[13] = m1.mat4[13] + value;
    m.mat4[14] = m1.mat4[14] + value;
    m.mat4[15] = m1.mat4[15] + value;
    //m is General matrix...
    return m;
}

/*!
    Adds the \a m1 and \a m2 matrix and returns the result.
*/
Matrix4x4 operator+(const double& value, const Matrix4x4& m1){
    return m1 + value;
}


/*!
    Subtracts the \a other matrix from itself and returns itself.
*/
Matrix4x4& Matrix4x4::operator-=(const Matrix4x4& other){
    *this = *this - other;
    return *this;
}

/*!
    Subtracts \a m2 from \a m1 and returns the result.
*/
Matrix4x4 operator-(const Matrix4x4& m1, const Matrix4x4& m2){
    Matrix4x4 m(1);
    m.mat4[0] = m1.mat4[0] - m2.mat4[0];
    m.mat4[1] = m1.mat4[1] - m2.mat4[1];
    m.mat4[2] = m1.mat4[2] - m2.mat4[2];
    m.mat4[3] = m1.mat4[3] - m2.mat4[3];
    m.mat4[4] = m1.mat4[4] - m2.mat4[4];
    m.mat4[5] = m1.mat4[5] - m2.mat4[5];
    m.mat4[6] = m1.mat4[6] - m2.mat4[6];
    m.mat4[7] = m1.mat4[7] - m2.mat4[7];
    m.mat4[8] = m1.mat4[8] - m2.mat4[8];
    m.mat4[9] = m1.mat4[9] - m2.mat4[9];
    m.mat4[10] = m1.mat4[10] - m2.mat4[10];
    m.mat4[11] = m1.mat4[11] - m2.mat4[11];
    m.mat4[12] = m1.mat4[12] - m2.mat4[12];
    m.mat4[13] = m1.mat4[13] - m2.mat4[13];
    m.mat4[14] = m1.mat4[14] - m2.mat4[14];
    m.mat4[15] = m1.mat4[15] - m2.mat4[15];
    //m is General matrix...
    return m;
}


/*!
    Multiplys the \a other matrix to itself and returns itself.
*/
Matrix4x4& Matrix4x4::operator*=(const Matrix4x4& other){
    *this = *this * other;
    return *this;
}

/*!
    Multiplys \a m1 with \a m2 and returns the result.
*/
Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2){

    if(m1.flagBits == Matrix4x4::Identity)
        return m2;
    else if(m2.flagBits == Matrix4x4::Identity)
        return m1;

    Matrix4x4 m(1);
    // Fisrt Column
    m.mat4[0] = m1.mat4[0]*m2.mat4[0] + m1.mat4[4]*m2.mat4[1] + m1.mat4[8]*m2.mat4[2] + m1.mat4[12]*m2.mat4[3];
    m.mat4[1] = m1.mat4[1]*m2.mat4[0] + m1.mat4[5]*m2.mat4[1] + m1.mat4[9]*m2.mat4[2] + m1.mat4[13]*m2.mat4[3];
    m.mat4[2] = m1.mat4[2]*m2.mat4[0] + m1.mat4[6]*m2.mat4[1] + m1.mat4[10]*m2.mat4[2] + m1.mat4[14]*m2.mat4[3];
    m.mat4[3] = m1.mat4[3]*m2.mat4[0] + m1.mat4[7]*m2.mat4[1] + m1.mat4[11]*m2.mat4[2] + m1.mat4[15]*m2.mat4[3];

    // Second Column
    m.mat4[4] = m1.mat4[0]*m2.mat4[4] + m1.mat4[4]*m2.mat4[5] + m1.mat4[8]*m2.mat4[6] + m1.mat4[12]*m2.mat4[7];
    m.mat4[5] = m1.mat4[1]*m2.mat4[4] + m1.mat4[5]*m2.mat4[5] + m1.mat4[9]*m2.mat4[6] + m1.mat4[13]*m2.mat4[7];
    m.mat4[6] = m1.mat4[2]*m2.mat4[4] + m1.mat4[6]*m2.mat4[5] + m1.mat4[10]*m2.mat4[6] + m1.mat4[14]*m2.mat4[7];
    m.mat4[7] = m1.mat4[3]*m2.mat4[4] + m1.mat4[7]*m2.mat4[5] + m1.mat4[11]*m2.mat4[6] + m1.mat4[15]*m2.mat4[7];

    // Third Column
    m.mat4[8] = m1.mat4[0]*m2.mat4[8] + m1.mat4[4]*m2.mat4[9] + m1.mat4[8]*m2.mat4[10] + m1.mat4[12]*m2.mat4[11];
    m.mat4[9] = m1.mat4[1]*m2.mat4[8] + m1.mat4[5]*m2.mat4[9] + m1.mat4[9]*m2.mat4[10] + m1.mat4[13]*m2.mat4[11];
    m.mat4[10] = m1.mat4[2]*m2.mat4[8] + m1.mat4[6]*m2.mat4[9] + m1.mat4[10]*m2.mat4[10] + m1.mat4[14]*m2.mat4[11];
    m.mat4[11] = m1.mat4[3]*m2.mat4[8] + m1.mat4[7]*m2.mat4[9] + m1.mat4[11]*m2.mat4[10] + m1.mat4[15]*m2.mat4[11];

    // Fourth Column
    m.mat4[12] = m1.mat4[0]*m2.mat4[12] + m1.mat4[4]*m2.mat4[13] + m1.mat4[8]*m2.mat4[14] + m1.mat4[12]*m2.mat4[15];
    m.mat4[13] = m1.mat4[1]*m2.mat4[12] + m1.mat4[5]*m2.mat4[13] + m1.mat4[9]*m2.mat4[14] + m1.mat4[13]*m2.mat4[15];
    m.mat4[14] = m1.mat4[2]*m2.mat4[12] + m1.mat4[6]*m2.mat4[13] + m1.mat4[10]*m2.mat4[14] + m1.mat4[14]*m2.mat4[15];
    m.mat4[15] = m1.mat4[3]*m2.mat4[12] + m1.mat4[7]*m2.mat4[13] + m1.mat4[11]*m2.mat4[14] + m1.mat4[15]*m2.mat4[15];
    //m is General matrix...
    return m;
}

/*!
    Multiplys the matrix with a \a multiplier. Scales all matrix values.
*/
Matrix4x4& Matrix4x4::operator*=(const double& multiplier){
    *this = *this * multiplier;
    return *this;
}

/*!
    Multiplys \a m1 with a \a multiplier and returns the result. Scales all matrix values.
*/
Matrix4x4 operator*(const Matrix4x4& m1, const double& multiplier){


    if(m1.flagBits == Matrix4x4::Identity){
        Matrix4x4 m_i;
        m_i.mat4[0]    = m_i.mat4[5]  = m_i.mat4[10] = m_i.mat4[15] = multiplier;
        return m_i;
    }


    Matrix4x4 m(1);

    // Fisrt Column
    m.mat4[0] = m1.mat4[0] * multiplier;
    m.mat4[1] = m1.mat4[1] * multiplier;
    m.mat4[2] = m1.mat4[2] * multiplier;
    m.mat4[3] = m1.mat4[3] * multiplier;

    // Second Column
    m.mat4[4] = m1.mat4[4] * multiplier;
    m.mat4[5] = m1.mat4[5] * multiplier;
    m.mat4[6] = m1.mat4[6] * multiplier;
    m.mat4[7] = m1.mat4[7] * multiplier;

    // Third Column
    m.mat4[8] = m1.mat4[8] * multiplier;
    m.mat4[9] = m1.mat4[9] * multiplier;
    m.mat4[10] = m1.mat4[10] * multiplier;
    m.mat4[11] = m1.mat4[11] * multiplier;

    // Fourth Column
    m.mat4[12] = m1.mat4[12] * multiplier;
    m.mat4[13] = m1.mat4[13] * multiplier;
    m.mat4[14] = m1.mat4[14] * multiplier;
    m.mat4[15] = m1.mat4[15] * multiplier;
    //m is General matrix...
    return m;
}

/*!
    Multiplys \a m1 with a \a multiplier and returns the result. Scales all matrix values.
    Actually just a kommutative way of the above function.

    \sa operator*(const Matrix4x4& m1, const double& multiplier)
*/
Matrix4x4 operator*(const double& multiplier, const Matrix4x4& m1){
    return m1 * multiplier;
}


/*!
    Divides the matrix with a \a divisor. Scales all matrix values.
*/
Matrix4x4& Matrix4x4::operator/=(const double& divisor){
    *this = *this / divisor;
    return *this;
}

/*!
    Divides \a m1 with a \a divisor and returns the result. Scales all matrix values.
*/
Matrix4x4 operator/(const Matrix4x4& m1, const double& divisor){


    if(m1.flagBits == Matrix4x4::Identity){
        Matrix4x4 m_i;
        m_i.mat4[0]    = m_i.mat4[5]  = m_i.mat4[10] = m_i.mat4[15] = divisor;
        return m_i;
    }


    Matrix4x4 m(1);

    // Fisrt Column
    m.mat4[0] = m1.mat4[0] / divisor;
    m.mat4[1] = m1.mat4[1] / divisor;
    m.mat4[2] = m1.mat4[2] / divisor;
    m.mat4[3] = m1.mat4[3] / divisor;

    // Second Column
    m.mat4[4] = m1.mat4[4] / divisor;
    m.mat4[5] = m1.mat4[5] / divisor;
    m.mat4[6] = m1.mat4[6] / divisor;
    m.mat4[7] = m1.mat4[7] / divisor;

    // Third Column
    m.mat4[8] = m1.mat4[8] / divisor;
    m.mat4[9] = m1.mat4[9] / divisor;
    m.mat4[10] = m1.mat4[10] / divisor;
    m.mat4[11] = m1.mat4[11] / divisor;

    // Fourth Column
    m.mat4[12] = m1.mat4[12] / divisor;
    m.mat4[13] = m1.mat4[13] / divisor;
    m.mat4[14] = m1.mat4[14] / divisor;
    m.mat4[15] = m1.mat4[15] / divisor;
    //m is General matrix...
    return m;
}



//Vector4
Vector4 operator*(const Vector4& vector, const Matrix4x4& matrix)
{
    double x, y, z, w;
    x = vector.x() * matrix.mat4[0] +
        vector.y() * matrix.mat4[1] +
        vector.z() * matrix.mat4[2] +
        vector.w() * matrix.mat4[3];
    y = vector.x() * matrix.mat4[4] +
        vector.y() * matrix.mat4[5] +
        vector.z() * matrix.mat4[6] +
        vector.w() * matrix.mat4[7];
    z = vector.x() * matrix.mat4[8] +
        vector.y() * matrix.mat4[9] +
        vector.z() * matrix.mat4[10] +
        vector.w() * matrix.mat4[11];
    w = vector.x() * matrix.mat4[12] +
        vector.y() * matrix.mat4[13] +
        vector.z() * matrix.mat4[14] +
        vector.w() * matrix.mat4[15];
    return Vector4(x, y, z, w);
}

Vector4 operator*(const Matrix4x4& matrix, const Vector4& vector)
{
    double x, y, z, w;
    x = vector.x() * matrix.mat4[0] +
        vector.y() * matrix.mat4[4] +
        vector.z() * matrix.mat4[8] +
        vector.w() * matrix.mat4[12];
    y = vector.x() * matrix.mat4[1] +
        vector.y() * matrix.mat4[5] +
        vector.z() * matrix.mat4[9] +
        vector.w() * matrix.mat4[13];
    z = vector.x() * matrix.mat4[2] +
        vector.y() * matrix.mat4[6] +
        vector.z() * matrix.mat4[10] +
        vector.w() * matrix.mat4[14];
    w = vector.x() * matrix.mat4[3] +
        vector.y() * matrix.mat4[7] +
        vector.z() * matrix.mat4[11] +
        vector.w() * matrix.mat4[15];
    return Vector4(x, y, z, w);
}

//Vector3
Vector3 operator*(const Vector3& vector, const Matrix4x4& matrix)
{
    double x, y, z, w;
    x = vector.x() * matrix.mat4[0] +
        vector.y() * matrix.mat4[1] +
        vector.z() * matrix.mat4[2] +
        matrix.mat4[3];
    y = vector.x() * matrix.mat4[4] +
        vector.y() * matrix.mat4[5] +
        vector.z() * matrix.mat4[6] +
        matrix.mat4[7];
    z = vector.x() * matrix.mat4[8] +
        vector.y() * matrix.mat4[9] +
        vector.z() * matrix.mat4[10] +
        matrix.mat4[11];
    w = vector.x() * matrix.mat4[12] +
        vector.y() * matrix.mat4[13] +
        vector.z() * matrix.mat4[14] +
        matrix.mat4[15];
    if (w == 1.0f)
        return Vector3(x, y, z);
    else
        return Vector3(x / w, y / w, z / w);
}

Vector3 operator*(const Matrix4x4& matrix, const Vector3& vector)
{
    double x, y, z, w;
    if (matrix.flagBits == Matrix4x4::Identity) {
        return vector;
    } else if (matrix.flagBits == Matrix4x4::Translation) {
        return Vector3(vector.x() + matrix.mat4[12],
                       vector.y() + matrix.mat4[13],
                       vector.z() + matrix.mat4[14]);
    } else if (matrix.flagBits ==
                    (Matrix4x4::Translation | Matrix4x4::Scale)) {
        return Vector3(vector.x() * matrix.mat4[0] + matrix.mat4[12],
                       vector.y() * matrix.mat4[5] + matrix.mat4[13],
                       vector.z() * matrix.mat4[10] + matrix.mat4[14]);
    } else if (matrix.flagBits == Matrix4x4::Scale) {
        return Vector3(vector.x() * matrix.mat4[0],
                       vector.y() * matrix.mat4[5],
                       vector.z() * matrix.mat4[10]);
    } else {
        x = vector.x() * matrix.mat4[0] +
            vector.y() * matrix.mat4[4] +
            vector.z() * matrix.mat4[8] +
            matrix.mat4[12];
        y = vector.x() * matrix.mat4[1] +
            vector.y() * matrix.mat4[5] +
            vector.z() * matrix.mat4[9] +
            matrix.mat4[13];
        z = vector.x() * matrix.mat4[2] +
            vector.y() * matrix.mat4[6] +
            vector.z() * matrix.mat4[10] +
            matrix.mat4[14];
        w = vector.x() * matrix.mat4[3] +
            vector.y() * matrix.mat4[7] +
            vector.z() * matrix.mat4[11] +
            matrix.mat4[15];
        if (w == 1.0f)
            return Vector3(x, y, z);
        else
            return Vector3(x / w, y / w, z / w);
    }
}
