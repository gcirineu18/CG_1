#ifndef CAMERA_HPP
#define CAMERA_HPP

#include "vec3.hpp"
#include "ray.hpp"

class Camera {
public:
    Vec3 eye, up, para;
    float d;
    float xmin, xmax, ymin, ymax;

    Vec3 u, v, w;

    Camera(Vec3 _eye, Vec3 _para, Vec3 _up, float _d, float x0, float x1, float y0, float y1) {
        eye = _eye;
        para = _para;
        up = _up;
        d = _d;
        xmin = x0;
        xmax = x1;
        ymin = y0;
        ymax = y1;

        Vec3 view = (para - eye);
        if (view.length() < 1e-6) view = Vec3(0, 0, -1);

        w = view.normalize();

        Vec3 right = w.cross(up);
        if (right.length() < 1e-6) {
            // Se o 'up' for paralelo ao 'w', forçamos um vetor lateral
            right = Vec3(1, 0, 0);
        }

        u = right.normalize();
        v = u.cross(w).normalize();
    }

    Ray gerarRaio(float s, float t) {
        float u_coord = xmin + s * (xmax - xmin);
        float v_coord = ymin + t * (ymax - ymin);

        Vec3 dr = (w * d) + (u * u_coord) + (v * v_coord);
        return Ray{eye, dr.normalize()};
    }
};

#endif // CAMERA_HPP
