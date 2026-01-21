#ifndef LIGHT_HPP
#define LIGHT_HPP

#include "vec3.hpp"

struct Light {
    Vec3 position;
    Vec3 color;

    Light(Vec3 p, Vec3 c) : position(p), color(c) {}
};

#endif // LIGHT_HPP
