#ifndef OBJECT_HPP
#define OBJECT_HPP

#include "ray.hpp"
#include "commons.hpp"

class Object {
public:
    Material mat;
    int id;

    virtual bool intersect(const Ray& r, HitRecord& rec) = 0;

    virtual ~Object() {}
};

#endif
