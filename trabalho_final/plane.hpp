#ifndef PLANE_HPP
#define PLANE_HPP

#include "object.hpp"

class Plane : public Object {
public:
    Vec3 p_pi;
    Vec3 n_bar;

    Plane(Vec3 p, Vec3 n, Material m, int _id){
        p_pi = p;
        n_bar = n;
        mat = m;
        id = _id;
    }

    bool intersect(const Ray& r, HitRecord& rec) override {
        float drn = r.direction.dot(n_bar);

        if (std::abs(drn) < 1e-6) return false;

        Vec3 w_p = r.origin - p_pi;
        float t_hit = -(w_p.dot(n_bar)) / drn;

        if (t_hit > 0.001f) {
            rec.t = t_hit;
            rec.p = r.pt(t_hit);
            rec.normal = n_bar;

            if (drn > 0) {
                rec.normal = n_bar * -1.0f;
            }

            rec.objectID = id;
            rec.mat = mat;
            return true;
        }

        return false;
    }
};

#endif // PLANE_HPP
