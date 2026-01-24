#ifndef CENA_HPP
#define CENA_HPP
using namespace std;

#include <vector>
#include "object.hpp"

class Cena{
public:
    ~Cena(){}
    
    std::vector<Object*> objetos;
    Vec3 bg_color = Vec3(0.2f, 0.2f, 0.2f);

    void adicionar(Object* obj) {
        objetos.push_back(obj);
    }

    bool trace(const Ray& r, HitRecord& melhorHit) {
        bool atingiuAlgo = false;
        float t_menor = 1e30f;

        for (auto obj : objetos) {
            HitRecord tempHit;
            if (obj->intersect(r, tempHit)){
                if (tempHit.t < t_menor) {
                    t_menor = tempHit.t;
                    melhorHit = tempHit;
                    atingiuAlgo = true;
                }
            }
        }
        return atingiuAlgo;
    }
};

#endif // CENA_HPP
