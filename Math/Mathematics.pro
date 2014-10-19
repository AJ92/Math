#-------------------------------------------------
#
# Project created by QtCreator 2014-10-19T19:50:51
#
#-------------------------------------------------

QT       -= gui

TARGET = Mathematics
TEMPLATE = lib

DEFINES += MATHEMATICS_LIBRARY

SOURCES += \
    Matrix/matrix3x3.cpp \
    Matrix/matrix4x4.cpp \
    Vector/vector3.cpp \
    Vector/vector4.cpp \
    Intersections/intersections.cpp \
    Geometry/aabb.cpp \
    Geometry/plane.cpp \
    Geometry/sphere.cpp

HEADERS +=\
    mathematics.h \
    Matrix/matrix3x3.h \
    Matrix/matrix4x4.h \
    Vector/vector3.h \
    Vector/vector4.h \
    Intersections/intersections.h \
    Geometry/aabb.h \
    Geometry/plane.h \
    Geometry/sphere.h \
    mathematics_global.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
