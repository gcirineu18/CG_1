#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP

#include "object.hpp"
#include "vec3.hpp"

class Triangle : public Object {
public:
    Vec3 v0, v1, v2;  // Vértices do triângulo
    Vec3 normal;

    Triangle(Vec3 _v0, Vec3 _v1, Vec3 _v2, Material m) 
        :Object(false), v0(_v0), v1(_v1), v2(_v2) {
        mat = m;
        
        // Calcular normal usando produto cruzado
        Vec3 e1 = v1 - v0;
        Vec3 e2 = v2 - v0;
        normal = e1.cross(e2).normalize();
    }

    bool intersect(const Ray& r, HitRecord& rec) override {
        const float EPSILON = 0.0001f;
        
        Vec3 e1 = v1 - v0;
        Vec3 e2 = v2 - v0;
        Vec3 h = r.direction.cross(e2);
        float a = e1.dot(h);

        if (fabs(a) < EPSILON) return false;

        float f = 1.0f / a;
        Vec3 s = r.origin - v0;
        float u = f * s.dot(h);

        if (u < 0.0f || u > 1.0f) return false;

        Vec3 q = s.cross(e1);
        float v = f * r.direction.dot(q);

        if (v < 0.0f || u + v > 1.0f) return false;

        float t = f * e2.dot(q);

        if (t > 0.001f) {
            rec.t = t;
            rec.p = r.pt(t);
            rec.normal = normal;
            rec.objectID = id;
            rec.mat = mat;
            rec.u = u;
            rec.v = v;
            return true;
        }
        return false;
    }
};

#endif // TRIANGLE_HPP