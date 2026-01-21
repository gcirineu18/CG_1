#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "object.hpp"

class Sphere : public Object {
public:
    Vec3 centro;
    float raio;

    Sphere(Vec3 c, float r, Material m, const char* nome){
        centro = c;
        raio = r;
        mat = m;
        name = nome;
    }

    bool intersect(const Ray& r, HitRecord& rec) override {
        Vec3 w = r.origin - centro;

        float a = r.direction.dot(r.direction);
        float b = 2.0f * w.dot(r.direction);
        float c = w.dot(w) - (raio * raio);

        float delta = b * b - (4.0f * a * c);

        if (delta < 0) return false;

        float sqrt_delta = sqrt(delta);
        float t1 = (-b - sqrt_delta) / (2.0f * a);
        float t2 = (-b + sqrt_delta) / (2.0f * a);

        float t_hit = -1.0f;
        if (t1 > 0.001f) t_hit = t1;
        else if (t2 > 0.001f) t_hit = t2;

        if (t_hit > 0) {
            rec.t = t_hit;
            rec.p = r.pt(t_hit);
            rec.normal = (rec.p - centro).normalize();
            rec.objectID = id;
            rec.mat = mat;
            return true;
        }
        return false;
    }
};

#endif
