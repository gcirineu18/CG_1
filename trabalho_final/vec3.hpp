#ifndef VEC3_HPP
#define VEC3_HPP

#include <cmath>

struct Vec3 {
    double x, y, z;

    Vec3(double x=0, double y=0, double z=0): x(x),y(y),z(z) {}

    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }

    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }

    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

    Vec3 operator*(const Vec3& v) const { return {x * v.x, y * v.y, z * v.z}; }

    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

    Vec3 cross(const Vec3& v) const {
        return {y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x};
    }

    double length() const { return sqrt(dot(*this)); }

    Vec3 normalize() const { return (*this) * (1.0 / length()); }
};

#endif
