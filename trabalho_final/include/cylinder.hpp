#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include "object.hpp"

class Cylinder : public Object {
public:
    Vec3 B;
    Vec3 u;
    float R;
    float H;
    bool tampaBase;
    bool tampaTopo;

    Cylinder(Vec3 base, Vec3 d_eixo, float raio, float altura, Material m, 
        bool baseOn = true, bool topoOn = true, const char* nome = "") {
        B = base;
        u = d_eixo.normalize();
        R = raio;
        H = altura;
        mat = m;
        tampaBase = baseOn;
        tampaTopo = topoOn;
        name = nome;
    }

    bool intersect(const Ray& r, HitRecord& rec) override {
        Vec3 dp = r.origin - B;

        Vec3 w = r.direction - u * (r.direction.dot(u));
        Vec3 v = dp - u * (dp.dot(u));

        float a = w.dot(w);
        float b = 2.0f * v.dot(w);
        float c = v.dot(v) - (R * R);

        float delta = b * b - 4.0f * a * c;
        float t_corpo = -1.0f;

        // Se a for 0, raio paralelo ao eixo
        if ( std::abs(a) > 1e-6 && delta >= 0) {
            float sqrt_delta = std::sqrt(delta);
            float t1 = (-b - sqrt_delta)/(2.0f * a);
            float t2 = (-b + sqrt_delta)/(2.0f * a);

            for(float t_teste : {t1, t2}) {
                if (t_teste > 0.001f) {
                    Vec3 P = r.pt(t_teste);

                    float h_check = (P - B).dot(u);

                    if (h_check >= 0 && h_check <= H) {
                        t_corpo = t_teste;
                        rec.t = t_corpo;
                        rec.p = P;

                        Vec3 p_eixo = B + u * h_check;
                        rec.normal = (rec.p - p_eixo).normalize();
                        break;
                    }
                }
            }
        }

        // Checar para as tampas - mesmo atingindo algum ponto no corpo, a tampa pode estar mais perto
        float t_tampas = -1.0f;
        HitRecord rec_tampas;
        if (checkCaps(r, rec_tampas)) {
            t_tampas = rec_tampas.t;
        }

        // Se tiver interseção no corpo e na tampa, pegamos o menor
        if (t_corpo > 0.001f && (t_tampas < 0 || t_corpo < t_tampas)){
            rec.mat = mat;
            rec.objectID = id;
            return true;
        } else if (t_tampas > 0.001f) {
            rec = rec_tampas;
            return true;
        }

        return false;
    }

    float intersectDisc(const Ray& r, Vec3 centro_disco, Vec3 n_disco) {
        float denom = r.direction.dot(n_disco);
        if (std::abs(denom) < 1e-6) return -1.0f;

        float t = (centro_disco - r.origin).dot(n_disco) / denom;
        if (t < 0.001f) return -1.0f;

        Vec3 p = r.pt(t);
        Vec3 dist_vec = p - centro_disco;
        if (dist_vec.dot(dist_vec) <= R*R) {
            return t;
        }
        return -1.0f;
    }

private:
    bool checkCaps(const Ray& r, HitRecord& rec_cap) {
        float t_base = -1.0f;
        float t_topo = -1.0f;

        if (tampaBase) {
            t_base = intersectDisc(r, B, u * -1.0f);
        }

        if (tampaTopo) {
            t_topo = intersectDisc(r, B + u * H, u);
        }

        float t_hit = -1.0f;
        Vec3 n_escolhida;

        if (t_base > 0.001f && t_topo > 0.001f) {
            if (t_base < t_topo) {
                t_hit = t_base;
                n_escolhida = u * -1.0f;
            } else {
                t_hit = t_topo;
                n_escolhida = u;
            }
        } else if (t_base > 0.001f) {
            t_hit = t_base;
            n_escolhida = u * -1.0f;
        } else if (t_topo > 0.001f) {
            t_hit = t_topo;
            n_escolhida = u;
        }

        if (t_hit > 0.001f) {
            rec_cap.t = t_hit;
            rec_cap.p = r.pt(t_hit);
            rec_cap.normal = n_escolhida;
            rec_cap.mat = mat;
            rec_cap.objectID = id;
            return true;
        }

        return false;
    }
};

#endif // CYLINDER_HPP
