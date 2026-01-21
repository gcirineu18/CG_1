#ifndef OBJECT_HPP
#define OBJECT_HPP

#include "ray.hpp"
#include "commons.hpp"

class Object {
public:
    Material mat;
    int id;   
    const char* name;

    static int nextID;

    virtual bool intersect(const Ray& r, HitRecord& rec) = 0;

    virtual ~Object() {}

    Object(bool auto_increment = true){
        if (auto_increment) {
            id = nextID++;
        }
    }
};

inline int Object::nextID = 0;
#endif

