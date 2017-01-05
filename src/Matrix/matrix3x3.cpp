#include "Matrix/matrix3x3.h"

/*!
    Constructs the 3x3 Matrix and sets it to it's Identity.

    \sa set_to_identity()
*/
Matrix3x3::Matrix3x3()
{
    set_to_identity();
}

/*!
    \overload

    Constructs the 3x3 Matrix and copies the
    components from \a mat in it's own 3x3 Matrix.

    \sa set_array()

*/
Matrix3x3::Matrix3x3(const Matrix3x3 &mat){
    this->mat3[0] = mat.mat3[0];
    this->mat3[1] = mat.mat3[1];
    this->mat3[2] = mat.mat3[2];
    this->mat3[3] = mat.mat3[3];
    this->mat3[4] = mat.mat3[4];
    this->mat3[5] = mat.mat3[5];
    this->mat3[6] = mat.mat3[6];
    this->mat3[7] = mat.mat3[7];
    this->mat3[8] = mat.mat3[8];
    flagBits = General;
}

/*!
    \overload

    Constructs the 3x3 Matrix and sets the components to \a f1 - \a f9.
*/
Matrix3x3::Matrix3x3(double f1, double f2, double f3,
                     double f4, double f5, double f6,
                     double f7, double f8, double f9)
{
    mat3[0] = f1;
    mat3[1] = f2;
    mat3[2] = f3;
    mat3[3] = f4;
    mat3[4] = f5;
    mat3[5] = f6;
    mat3[6] = f7;
    mat3[7] = f8;
    mat3[8] = f9;

    flagBits = General;
}

/*!
    \overload

    Constructs the 4x4 Matrix and sets the
    components from \a mat in it's own 3x3 Matrix.

    \sa set_array()

*/
Matrix3x3::Matrix3x3(const double *mat)
{
    for (int i = 0; i < 9; i++){
        mat3[i] = mat[i];
    }
    flagBits = General;
}

/*!

    Sets the matrix to it's Identity.

*/
void Matrix3x3::set_to_identity(){
    //setting matrix to it's Identity
    mat3[0]    = mat3[4]  = mat3[8] = 1.0;
    mat3[1]    = mat3[2]  = mat3[3] = 0.0;
    mat3[5]    = mat3[6]  = mat3[7]  = 0.0;
    flagBits = Identity;
}

/*!
    Checks if the matrix is an identity matrix.
  */
bool Matrix3x3::is_identity(){
    if (flagBits == Identity)
        return true;
    if (mat3[0] != 1.0 || mat3[1] != 0.0 || mat3[2] != 0.0)
        return false;
    if (mat3[3] != 0.0 || mat3[4] != 1.0 || mat3[5] != 0.0)
        return false;
    if (mat3[6] != 0.0 || mat3[7] != 0.0)
        return false;
    return (mat3[8] == 1.0f);
}









/*!
    Sets the given \a value at the \a index in the matrix.
*/
void Matrix3x3::set_value(int index, double value){
    if(index < 0 || index > 8){
        return;
    }
    mat3[index] = value;
    flagBits = General;
}

/*!
    Returns the value at the given \a index.
*/
double Matrix3x3::get_value(int index){
    if(index < 0 || index > 8){
        return 0;
    }
    return mat3[index];
}

/*!
    Returns the matrix as an array of doubles.
*/
double* Matrix3x3::get_array(){
    return mat3;
}

/*!
    Sets the components of the matrix to the values stored int the \a mat3 array.
*/
void Matrix3x3::set_array(double mat3[]){
    this->mat3[0] = mat3[0];
    this->mat3[1] = mat3[1];
    this->mat3[2] = mat3[2];
    this->mat3[3] = mat3[3];
    this->mat3[4] = mat3[4];
    this->mat3[5] = mat3[5];
    this->mat3[6] = mat3[6];
    this->mat3[7] = mat3[7];
    this->mat3[8] = mat3[8];
    flagBits = General;
}

/*!
    Outputs the matrix in the console by using Qt's qDebug() function.
*/
void Matrix3x3::debug(){
    //qDebug("Matrix:");
    //qDebug("%f  %f  %f",mat3[0],mat3[1],mat3[2]);
    //qDebug("%f  %f  %f",mat3[3],mat3[4],mat3[5]);
    //qDebug("%f  %f  %f",mat3[6],mat3[7],mat3[8]);

}


//operators


/*!
    For easier access to the values of the matrix you can use this function.
*/
const double& Matrix3x3::operator[](int index) const
{
    if(index < 0 && index > 8){
        //qDebug("const double& Matrix3x3::operator()(int index) has a wrong index: %i", index);
    }
    return mat3[index];
}

/*!
    For easier access to the values of the matrix you can use this function.
*/
double& Matrix3x3::operator[](int index)
{
    if(index < 0 && index > 8){
        //qDebug("double& Matrix3x3::operator()(int index) has a wrong index: %i", index);
    }
    flagBits = General;
    return mat3[index];
}

/*!
    For easier access to the values of the matrix you can use this function.
*/
const double& Matrix3x3::operator()(int index) const
{
    if(index < 0 && index > 8){
        //qDebug("const double& Matrix3x3::operator()(int index) has a wrong index: %i", index);
    }
    return mat3[index];
}

/*!
    For easier access to the values of the matrix you can use this function.
*/
double& Matrix3x3::operator()(int index)
{
    if(index < 0 && index > 8){
        //qDebug("double& Matrix3x3::operator()(int index) has a wrong index: %i", index);
    }
    flagBits = General;
    return mat3[index];
}

/*!
    For easier access to the values of the matrix you can use this function.
*/
const double& Matrix3x3::operator()(int aRow, int aColumn) const
{
    if(aRow < 0 && aRow > 2 && aColumn < 0 && aColumn > 2){
        //qDebug("const double& Matrix3x3::operator()(int aRow, int aColumn) has wrong indizes: %i  %i", aRow, aColumn);
    }
    return mat3[aColumn+aRow*3];
}

/*!
    For easier access to the values of the matrix you can use this function.
*/
double& Matrix3x3::operator()(int aRow, int aColumn)
{
    if(aRow < 0 && aRow > 2 && aColumn < 0 && aColumn > 2){
        //qDebug("double& Matrix3x3::operator()(int aRow, int aColumn) has wrong indizes: %i  %i", aRow, aColumn);
    }
    flagBits = General;
    return mat3[aColumn+aRow*3];
}


/*!
    Adds the \a other matrix to itself and returns itself.
*/
Matrix3x3& Matrix3x3::operator+=(const Matrix3x3& other){
    *this = *this + other;
    return *this;
}

/*!
    Adds the \a m1 and \a m2 matrix and returns the result.
*/
Matrix3x3 operator+(const Matrix3x3& m1, const Matrix3x3& m2){
    Matrix3x3 m(1);
    m.mat3[0] = m1.mat3[0] + m2.mat3[0];
    m.mat3[1] = m1.mat3[1] + m2.mat3[1];
    m.mat3[2] = m1.mat3[2] + m2.mat3[2];
    m.mat3[3] = m1.mat3[3] + m2.mat3[3];
    m.mat3[4] = m1.mat3[4] + m2.mat3[4];
    m.mat3[5] = m1.mat3[5] + m2.mat3[5];
    m.mat3[6] = m1.mat3[6] + m2.mat3[6];
    m.mat3[7] = m1.mat3[7] + m2.mat3[7];
    m.mat3[8] = m1.mat3[8] + m2.mat3[8];
    //m is General matrix...
    return m;
}



/*!
    Adds the \a value to itself and returns itself.
*/
Matrix3x3& Matrix3x3::operator+=(const double& value){
    *this = *this + value;
    return *this;
}

/*!
    Adds the \a value to \a m1 matrix and returns the result.
*/
Matrix3x3 operator+(const Matrix3x3& m1, const double& value){
    Matrix3x3 m(1);
    m.mat3[0] = m1.mat3[0] + value;
    m.mat3[1] = m1.mat3[1] + value;
    m.mat3[2] = m1.mat3[2] + value;
    m.mat3[3] = m1.mat3[3] + value;
    m.mat3[4] = m1.mat3[4] + value;
    m.mat3[5] = m1.mat3[5] + value;
    m.mat3[6] = m1.mat3[6] + value;
    m.mat3[7] = m1.mat3[7] + value;
    m.mat3[8] = m1.mat3[8] + value;
    //m is General matrix...
    return m;
}

/*!
    Adds the \a m1 and \a m2 matrix and returns the result.
*/
Matrix3x3 operator+(const double& value, const Matrix3x3& m1){
    return m1 + value;
}


/*!
    Subtracts the \a other matrix from itself and returns itself.
*/
Matrix3x3& Matrix3x3::operator-=(const Matrix3x3& other){
    *this = *this - other;
    return *this;
}

/*!
    Subtracts \a m2 from \a m1 and returns the result.
*/
Matrix3x3 operator-(const Matrix3x3& m1, const Matrix3x3& m2){
    Matrix3x3 m(1);
    m.mat3[0] = m1.mat3[0] - m2.mat3[0];
    m.mat3[1] = m1.mat3[1] - m2.mat3[1];
    m.mat3[2] = m1.mat3[2] - m2.mat3[2];
    m.mat3[3] = m1.mat3[3] - m2.mat3[3];
    m.mat3[4] = m1.mat3[4] - m2.mat3[4];
    m.mat3[5] = m1.mat3[5] - m2.mat3[5];
    m.mat3[6] = m1.mat3[6] - m2.mat3[6];
    m.mat3[7] = m1.mat3[7] - m2.mat3[7];
    m.mat3[8] = m1.mat3[8] - m2.mat3[8];
    //m is General matrix...
    return m;
}


/*!
    Multiplys the \a other matrix to itself and returns itself.
*/
Matrix3x3& Matrix3x3::operator*=(const Matrix3x3& other){
    *this = *this * other;
    return *this;
}

/*!
    Multiplys \a m1 with \a m2 and returns the result.
*/
Matrix3x3 operator*(const Matrix3x3& m1, const Matrix3x3& m2){

    if(m1.flagBits == Matrix3x3::Identity)
        return m2;
    else if(m2.flagBits == Matrix3x3::Identity)
        return m1;

    Matrix3x3 m(1);
    //seems to be transposed...
    /*
    // Fisrt Column
    m.mat3[0] = m1.mat3[0]*m2.mat3[0] + m1.mat3[1]*m2.mat3[3] + m1.mat3[2]*m2.mat3[6];
    m.mat3[1] = m1.mat3[0]*m2.mat3[1] + m1.mat3[1]*m2.mat3[4] + m1.mat3[2]*m2.mat3[7];
    m.mat3[2] = m1.mat3[0]*m2.mat3[2] + m1.mat3[1]*m2.mat3[5] + m1.mat3[2]*m2.mat3[8];

    // Second Column
    m.mat3[3] = m1.mat3[3]*m2.mat3[0] + m1.mat3[4]*m2.mat3[3] + m1.mat3[5]*m2.mat3[6];
    m.mat3[4] = m1.mat3[3]*m2.mat3[1] + m1.mat3[4]*m2.mat3[4] + m1.mat3[5]*m2.mat3[7];
    m.mat3[5] = m1.mat3[3]*m2.mat3[2] + m1.mat3[4]*m2.mat3[5] + m1.mat3[5]*m2.mat3[8];

    // Third Column
    m.mat3[6] = m1.mat3[6]*m2.mat3[0] + m1.mat3[7]*m2.mat3[3] + m1.mat3[8]*m2.mat3[6];
    m.mat3[7] = m1.mat3[6]*m2.mat3[1] + m1.mat3[7]*m2.mat3[4] + m1.mat3[8]*m2.mat3[7];
    m.mat3[8] = m1.mat3[6]*m2.mat3[2] + m1.mat3[7]*m2.mat3[5] + m1.mat3[8]*m2.mat3[8];
    */


    // Fisrt Column
    m.mat3[0] = m1.mat3[0]*m2.mat3[0] + m1.mat3[3]*m2.mat3[1] + m1.mat3[6]*m2.mat3[2];
    m.mat3[1] = m1.mat3[1]*m2.mat3[0] + m1.mat3[4]*m2.mat3[1] + m1.mat3[7]*m2.mat3[2];
    m.mat3[2] = m1.mat3[2]*m2.mat3[0] + m1.mat3[5]*m2.mat3[1] + m1.mat3[8]*m2.mat3[2];

    // Second Column
    m.mat3[3] = m1.mat3[0]*m2.mat3[3] + m1.mat3[3]*m2.mat3[4] + m1.mat3[6]*m2.mat3[5];
    m.mat3[4] = m1.mat3[1]*m2.mat3[3] + m1.mat3[4]*m2.mat3[4] + m1.mat3[7]*m2.mat3[5];
    m.mat3[5] = m1.mat3[2]*m2.mat3[3] + m1.mat3[5]*m2.mat3[4] + m1.mat3[8]*m2.mat3[5];

    // Third Column
    m.mat3[6] = m1.mat3[0]*m2.mat3[6] + m1.mat3[3]*m2.mat3[7] + m1.mat3[6]*m2.mat3[8];
    m.mat3[7] = m1.mat3[1]*m2.mat3[6] + m1.mat3[4]*m2.mat3[7] + m1.mat3[7]*m2.mat3[8];
    m.mat3[8] = m1.mat3[2]*m2.mat3[6] + m1.mat3[5]*m2.mat3[7] + m1.mat3[8]*m2.mat3[8];

    //m is General matrix...
    return m;
}

/*!
    Multiplys the matrix with a \a multiplier. Scales all matrix values.
*/
Matrix3x3& Matrix3x3::operator*=(const double& multiplier){
    *this = *this * multiplier;
    return *this;
}

/*!
    Multiplys \a m1 with a \a multiplier and returns the result. Scales all matrix values.
*/
Matrix3x3 operator*(const Matrix3x3& m1, const double& multiplier){


    if(m1.flagBits == Matrix3x3::Identity){
        Matrix3x3 m_i;
        m_i.mat3[0]    = m_i.mat3[4]  = m_i.mat3[8] = multiplier;
        return m_i;
    }


    Matrix3x3 m(1);

    // Fisrt Column
    m.mat3[0] = m1.mat3[0] * multiplier;
    m.mat3[1] = m1.mat3[1] * multiplier;
    m.mat3[2] = m1.mat3[2] * multiplier;

    // Second Column
    m.mat3[3] = m1.mat3[3] * multiplier;
    m.mat3[4] = m1.mat3[4] * multiplier;
    m.mat3[5] = m1.mat3[5] * multiplier;

    // Third Column
    m.mat3[6] = m1.mat3[6] * multiplier;
    m.mat3[7] = m1.mat3[7] * multiplier;
    m.mat3[8] = m1.mat3[8] * multiplier;

    //m is General matrix...
    return m;
}

/*!
    Multiplys \a m1 with a \a multiplier and returns the result. Scales all matrix values.
    Actually just a kommutative way of the above function.

    \sa operator*(const Matrix3x3& m1, const double& multiplier)
*/
Matrix3x3 operator*(const double& multiplier, const Matrix3x3& m1){
    return m1 * multiplier;
}


/*!
    Divides the matrix with a \a divisor. Scales all matrix values.
*/
Matrix3x3& Matrix3x3::operator/=(const double& divisor){
    *this = *this / divisor;
    return *this;
}

/*!
    Divides \a m1 with a \a divisor and returns the result. Scales all matrix values.
*/
Matrix3x3 operator/(const Matrix3x3& m1, const double& divisor){


    if(m1.flagBits == Matrix3x3::Identity){
        Matrix3x3 m_i;
        m_i.mat3[0] /= divisor;
        m_i.mat3[4] /= divisor;
        m_i.mat3[8] /= divisor;
        return m_i;
    }


    Matrix3x3 m(1);

    // Fisrt Column
    m.mat3[0] = m1.mat3[0] / divisor;
    m.mat3[1] = m1.mat3[1] / divisor;
    m.mat3[2] = m1.mat3[2] / divisor;

    // Second Column
    m.mat3[3] = m1.mat3[3] / divisor;
    m.mat3[4] = m1.mat3[4] / divisor;
    m.mat3[5] = m1.mat3[5] / divisor;

    // Third Column
    m.mat3[6] = m1.mat3[6] / divisor;
    m.mat3[7] = m1.mat3[7] / divisor;
    m.mat3[8] = m1.mat3[8] / divisor;

    //m is General matrix...
    return m;
}
