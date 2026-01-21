
#ifndef SPOTLIGHT_HPP
#define SPOTLIGHT_HPP
#include "vec3.hpp"
#include "light.hpp"
#include "commons.hpp"
#include "cena.hpp"
#include "camera.hpp"

class SpotLight : public Light{
public:
    Vec3 direction;
    float angle;

    SpotLight(Vec3 p, Vec3 c, Vec3 dir, float ang)
    : Light(p, c),  direction(dir.normalize()), angle(ang) 
    {
    }

    float getIntensity(const Vec3& pointToLight) const {

        float theta = direction.dot(pointToLight);
        if( theta > cos(angle)){
            return 1.0f;
        }
        return 0.0f;
    }

    static Vec3 calcula_luzSpot(vector<SpotLight> spotLights, HitRecord hit, Cena cena, Camera cam){
        Vec3 baseColor = hit.mat.color;

        if (hit.mat.useTexture && hit.mat.texture != nullptr) {
            baseColor = hit.mat.texture->sample(hit.u, hit.v);
        }
                
        Vec3 corFinal = Vec3(0, 0, 0);

        for (const auto& spot : spotLights) {
            Vec3 L = (spot.position - hit.p).normalize();
            double distL = (spot.position - hit.p).length();
            
            // Verificar se está dentro do cone
            float intensity = spot.getIntensity(L * -1.0); 
            
            if (intensity > 0.01f) {
                // Shadow ray
                Ray shadowRay(hit.p + hit.normal * 0.01, L);
                HitRecord shadowHit;
                
                bool emSombra = false;
                if (cena.trace(shadowRay, shadowHit)) {
                    if (shadowHit.objectID != hit.objectID && shadowHit.t < distL) 
                        emSombra = true;
                }
                
                if (!emSombra) {
                    double dotNL = std::max(0.0, hit.normal.dot(L));
                    
                    // Aplicar intensidade do spotlight
                    Vec3 lightContrib = (baseColor * spot.color) * 
                                    (hit.mat.kd * dotNL );
                    corFinal = corFinal + lightContrib;
                    
                    // Specular
                    Vec3 V_d = (cam.eye - hit.p).normalize();
                    Vec3 R_d = (hit.normal * (2.0 * hit.normal.dot(L))) - L;
                    double dotRV = std::max(0.0, R_d.dot(V_d));
                    
                    corFinal = corFinal + (spot.color ) * 
                            (hit.mat.ks * std::pow(dotRV, hit.mat.shininess));
                }
            }
        }
        return corFinal;
    }
};


#endif 