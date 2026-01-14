#include <iostream>
#include <fstream>
#include "vec3.hpp"
#include "ray.hpp"
#include "camera.hpp"
#include "cena.hpp"
#include "sphere.hpp"
#include "plane.hpp"
#include "cylinder.hpp"
#include "cone.hpp"
#include "light.hpp"

// Função para limitar os valores entre 0 e 1
double clamp(double x) {
    if (x < 0) return 0;
    if (x > 1) return 1;
    return x;
}

int main() {
    const int width = 500;
    const int height = 500;

    Camera cam(Vec3(0,0,0),Vec3(0,0,-1),Vec3(0,1,0),1,-1.0,1.0,-1.0,1.0);

    std::vector<Light> luzes;
    //luzes.push_back(Light(origin, color)
    luzes.push_back(Light(Vec3(-1,1.4,0.2),Vec3(0.7,0.7,0.7)));
    
    //Material material = {cor, ka, kd, ke, m}
    Material chao = {Vec3(0.2,0.2,0.2), 0.2,0.8,0.0,1.0};
    Material parede = {Vec3(0.2,0.2,0.8), 0.686, 0.933, 0.933,1.0};
    Material fundo = {Vec3(0.2,0.2,0.4), 0.686, 0.933, 0.933,1.0};
    Material teto = {Vec3(0.5,0.5,0.5), 0.933, 0.933, 0.933, 1.0};
    Material cil = {Vec3(0.39,0.12,0),0.824, 0.706, 0.549, 2.0};
    Material con = {Vec3(0.22,0.73,0.37),0.0, 1.0, 0.498, 2.0};
    Material esf = {Vec3(1.0, 0.73, 0.37), 0.854, 0.647, 0.125, 2.0};

    Cena cena;
    //cena.adicionar
    cena.adicionar(new Plane(Vec3(0,-1.5,0),Vec3(0,1,0),chao,1));
    cena.adicionar(new Plane(Vec3(2,-1.5,0),Vec3(-1,0,0),parede,2));
    cena.adicionar(new Plane(Vec3(2,-1.5,-4),Vec3(0,0,1),fundo,3));
    cena.adicionar(new Plane(Vec3(-2,-1.5,0),Vec3(1,0,0),parede,4));
    cena.adicionar(new Plane(Vec3(0,1.5,0),Vec3(0,-1,0),teto,5));

    cena.adicionar(new Cylinder(Vec3(0, -1.5,-2.0),Vec3(0,1,0),0.05, 0.9, cil,6));
    cena.adicionar(new Cone(Vec3(0, -0.6, -2), Vec3(0,1,0),0.9,1.5,con,7));
    cena.adicionar(new Sphere(Vec3(0, 0.95,-2.0),0.05,esf,8));

    std::ofstream file("resultado.ppm");
    file << "P3\n" << width << " " << height << "\n255\n";

    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            double u = (double)x /(width - 1);
            double v = (double)y /(height - 1);

            Ray r = cam.gerarRaio(u, v);
            HitRecord hit;
            Vec3 corPixel;

            if (cena.trace(r, hit)) {
                Vec3 corFinal = hit.mat.color * hit.mat.ka;     //Começamos com Ka

                for (const auto& luz : luzes){
                    Vec3 L = (luz.position - hit.p).normalize();
                    double distL = (luz.position - hit.p).length();

                    Ray shadowRay(hit.p + hit.normal * 0.01, L);
                    HitRecord shadowHit;

                    bool emSombra =  false;
                    if (cena.trace(shadowRay, shadowHit)){
                        if (shadowHit.objectID != hit.objectID && shadowHit.t < distL) emSombra = true;
                    }

                    if(!emSombra){

                        double dotNL = std::max(0.0, hit.normal.dot(L));
                        corFinal = corFinal + (hit.mat.color * luz.color)* (hit.mat.kd * dotNL);     //Soma com Kd

                        Vec3 V_d = (cam.eye - hit.p).normalize();
                        Vec3 R_d = (hit.normal * (2.0 * hit.normal.dot(L))) - L;
                        double dotRV = std::max(0.0, R_d.dot(V_d));

                        corFinal = corFinal + (luz.color) * (hit.mat.ks * std::pow(dotRV, hit.mat.shininess)); //Soma com Ke
                    }

                }
                corPixel = corFinal;

            } else {
                corPixel = cena.bg_color;
            }


            file << (int)(clamp(corPixel.x) * 255) << " "
                 << (int)(clamp(corPixel.y) * 255) << " "
                 << (int)(clamp(corPixel.z) * 255) << "\n";
        }
    }

    std::cout << "Imagem gerada com sucesso!" << std::endl;
    return 0;
}
