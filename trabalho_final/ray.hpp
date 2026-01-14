#ifndef RAY_HPP
#define RAY_HPP

#include "vec3.hpp"

struct Ray {
    Vec3 origin;
    Vec3 direction;

    Ray(const Vec3& o, const Vec3& d) : origin(o), direction(d) {}

    Vec3 pt(double t) const { return origin + direction * t; }
};

#endif
