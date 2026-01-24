#ifndef LIGHT_HPP
#define LIGHT_HPP
#include "vec3.hpp"
#include "textura.hpp" 

using namespace std;

struct Light {
    Vec3 position;
    Vec3 color;

    Light(Vec3 p, Vec3 c) : position(p), color(c) {}


    static Vec3 calcula_luzes(vector<Light> luzes, HitRecord hit, Cena cena, Camera cam){
        Vec3 baseColor = hit.mat.color;
                
        // Aplicar textura se disponível
        if (hit.mat.useTexture && hit.mat.texture != nullptr) {
            baseColor = hit.mat.texture->sample(hit.u, hit.v);
        }
                
        Vec3 corFinal = baseColor * hit.mat.ka;
        for (const auto& luz : luzes){
            Vec3 L = (luz.position - hit.p).normalize();
            float distL = (luz.position - hit.p).length();

            Ray shadowRay(hit.p + hit.normal * 0.01, L);
            HitRecord shadowHit;

            bool emSombra =  false;
            if (cena.trace(shadowRay, shadowHit)){
                if (shadowHit.objectID != hit.objectID && shadowHit.t < distL) emSombra = true;
            }

            if(!emSombra){
                float dotNL = std::max(0.0f, hit.normal.dot(L));
                corFinal = corFinal + (baseColor * luz.color)* (hit.mat.kd * dotNL);

                Vec3 V_d = (cam.eye - hit.p).normalize();
                Vec3 R_d = (hit.normal * (2.0 * hit.normal.dot(L))) - L;
                float dotRV = std::max(0.0f, R_d.dot(V_d));

                corFinal = corFinal + (luz.color) * (hit.mat.ks * std::pow(dotRV, hit.mat.shininess));
            }
        }
        return corFinal;
    }
};

#endif // LIGHT_HPP
