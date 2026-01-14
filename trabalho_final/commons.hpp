#ifndef COMMONS_HPP
#define COMMONS_HPP

#include "vec3.hpp"

struct Material {
    Vec3 color;
    double kd, ks, ka;
    double shininess;
};

struct HitRecord {
    double t;         // Distancia da interseção
    Vec3 p;           // Ponto exato no espaco (x,y,z)
    Vec3 normal;      // Vetor normal
    int objectID;
    Material mat;
};

#endif
