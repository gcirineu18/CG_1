#ifndef CONE_HPP
#define CONE_HPP

#include "object.hpp"

class Cone : public Object {
public:
    Vec3 cBase, V, n;
    float raio, H, theta;
    bool base;

    Cone(Vec3 base, Vec3 eixo, float r_param, float altura, Material m, bool baseOn = true, const char* nome = "") {
        cBase = base;
        n = eixo.normalize();
        raio = r_param;
        H = altura;
        theta = std::atan2(raio, altura);
        V = cBase + n * H;
        mat = m;
        base = baseOn;
        name = nome;
    }

    float intersectBase(const Ray& r) {
        float rBase = raio;
        float denom = r.direction.dot(n);

        if(std::abs(denom) < 1e-6) return -1.0f;

        float t = (cBase - r.origin).dot(n) / denom;
        if (t < 0.001) return -1.0f;

        Vec3 P = r.pt(t);
        if ((P - cBase).length() <= rBase) return t;
        return -1.0f;
    }

    bool intersect(const Ray& r, HitRecord& rec) override {
        Vec3 d = r.direction;
        Vec3 v = V - r.origin;
        float cos2 = cos(theta) * cos(theta);

        float dn = d.dot(n);
        float vn = v.dot(n);
        float vd = v.dot(d);
        float dd = d.dot(d);
        float vv = v.dot(v);

        float a = (dn * dn) - (dd * cos2);
        float b = (vd * cos2) - (vn * dn);
        float c = (vn * vn) - (vv * cos2);

        float t_hit = -1.0f;

        if (std::abs(a) < 1e-6) {
            if (std::abs(b) > 1e-6) t_hit = -c / (2.0 * b);
        } else {
            float delta = b * b - a * c;
            if (delta >= 0) {
                float sqrt_delta = std::sqrt(delta);
                float t1 = (-b - sqrt_delta)/ a;
                float t2 = (-b + sqrt_delta)/ a;

                if (t1 > t2) std::swap(t1, t2);

                for (float t_teste: {t1, t2}){
                    if (t_teste > 0.001){
                        Vec3 P = r.pt(t_teste);

                        float h_check = (V - P).dot(n);
                        if (h_check >= 0 && h_check <= H) {
                            t_hit = t_teste;
                            break;
                        }
                    }
                }
            }
        }

        float t_base = -1.0f;
        if (base) {
            t_base = intersectBase(r);
        }

        if (t_hit > 0.001 && (t_base <= 0.001 || t_hit < t_base)) {
            hitLateral(rec, r, t_hit);
            return true;
        } else if (t_base > 0.001) {
            hitBase(rec, r, t_base);
            return true;
        }

        return false;
    }

    void rotateX(float angle){
        cBase = cBase.rotX(angle);
        V = V.rotX(angle);
        n = (V - cBase).normalize();
        // recalcular ângulo de abertura do cone
        theta = std::atan2(raio, H);
    }

    void rotateY(float angle){
        cBase = cBase.rotY(angle);
        V = V.rotY(angle);
        n = (V - cBase).normalize();
        theta = std::atan2(raio, H);
    }

    void rotateZ(float angle){
        cBase = cBase.rotZ(angle);
        V = V.rotZ(angle);
        n = (V - cBase).normalize();
        theta = std::atan2(raio, H);
    }

    void translate( float tx, float ty, float tz){
        cBase = cBase.translate(tx, ty, tz);
        V = V.translate(tx, ty, tz);
    }

    void translate(const Vec3& t){
        translate(t.x, t.y, t.z);
    }

    void scaleTransform(float sx, float sy, float sz) {
        // Transformar base e vértice
        cBase = cBase.scale(sx, sy, sz);
        V = V.scale(sx, sy, sz);
        
        // Recalcular direção e altura
        Vec3 newDir = V - cBase;
        H = newDir.length();
        n = newDir.normalize();
        
        // O raio também precisa ser escalado
        raio = raio * sqrt(sx * sx + sz * sz) / 2.0;  // média das escalas perpendiculares
        theta = std::atan2(raio, H);
    }

    void scaleTransform(float s) {
        scaleTransform(s, s, s);
    }

    void shearTransform(float xy, float xz, float yx, float yz, float zx, float zy) {
        cBase = cBase.shear(xy, xz, yx, yz, zx, zy);
        V = V.shear(xy, xz, yx, yz, zx, zy);
        
        // Recalcular direção e altura
        Vec3 newDir = V - cBase;
        H = newDir.length();
        n = newDir.normalize();
        theta = std::atan2(raio, H);
    }


private:
    void hitLateral(HitRecord& rec, const Ray& r, float t){
        rec.t = t;
        rec.p = r.pt(t);

        Vec3 VP = rec.p - V;
        float h = VP.dot(n);
        float k = (raio*raio)/(H*H);

        Vec3 n_dir = VP - n * (h * (1.0 + k));
        rec.normal = n_dir.normalize();

        if (rec.normal.dot(r.direction) > 0) {
            rec.normal = rec.normal * -1.0;
        }

        rec.mat = mat;
        rec.objectID = id;
    }

    void hitBase(HitRecord& rec, const Ray& r, float t){
        rec.t = t;
        rec.p = r.pt(t);
        rec.normal = n * -1.0;
        rec.mat = mat;
        rec.objectID = id;
    }
};

#endif // CONE_HPP
