#ifndef CUBE_HPP
#define CUBE_HPP

#include "object.hpp"
#include "triangle.hpp"
#include <vector>

class Cube : public Object {
public:
    Vec3 min_corner;
    Vec3 max_corner;
    std::vector<Triangle> triangles;

    // Construtor: posição e tamanho do cubo
    Cube(Vec3 center, float size, Material m, int _id) {
        mat = m;
        id = _id;

        float half = size / 3.0f;
        min_corner = center - Vec3(half, half, half);
        max_corner = center + Vec3(half, half, half);

        createMesh(m, _id);
    }

    void createMesh(Material m, int _id) {
        // Define os 8 vértices do cubo
        Vec3 v0 = Vec3(min_corner.x, min_corner.y, min_corner.z);
        Vec3 v1 = Vec3(max_corner.x, min_corner.y, min_corner.z);
        Vec3 v2 = Vec3(max_corner.x, max_corner.y, min_corner.z);
        Vec3 v3 = Vec3(min_corner.x, max_corner.y, min_corner.z);
        Vec3 v4 = Vec3(min_corner.x, min_corner.y, max_corner.z);
        Vec3 v5 = Vec3(max_corner.x, min_corner.y, max_corner.z);
        Vec3 v6 = Vec3(max_corner.x, max_corner.y, max_corner.z);
        Vec3 v7 = Vec3(min_corner.x, max_corner.y, max_corner.z);

        // Face frontal (z mínimo) - 2 triângulos
        triangles.push_back(Triangle(v0, v1, v2, m, _id));
        triangles.push_back(Triangle(v0, v2, v3, m, _id));

        // Face traseira (z máximo) - 2 triângulos
        triangles.push_back(Triangle(v5, v4, v7, m, _id));
        triangles.push_back(Triangle(v5, v7, v6, m, _id));

        // Face esquerda (x mínimo) - 2 triângulos
        triangles.push_back(Triangle(v4, v0, v3, m, _id));
        triangles.push_back(Triangle(v4, v3, v7, m, _id));

        // Face direita (x máximo) - 2 triângulos
        triangles.push_back(Triangle(v1, v5, v6, m, _id));
        triangles.push_back(Triangle(v1, v6, v2, m, _id));

        // Face inferior (y mínimo) - 2 triângulos
        triangles.push_back(Triangle(v4, v5, v1, m, _id));
        triangles.push_back(Triangle(v4, v1, v0, m, _id));

        // Face superior (y máximo) - 2 triângulos
        triangles.push_back(Triangle(v3, v2, v6, m, _id));
        triangles.push_back(Triangle(v3, v6, v7, m, _id));
    }

    bool intersect(const Ray& r, HitRecord& rec) override {
        bool hit = false;
        float closest_t = 1e9f;

        // Testar interseção com todos os triângulos
        for (auto& tri : triangles) {
            HitRecord temp_rec;
            if (tri.intersect(r, temp_rec)) {
                if (temp_rec.t < closest_t) {
                    closest_t = temp_rec.t;
                    rec = temp_rec;
                    hit = true;
                }
            }
        }
        return hit;
    }
};

#endif // CUBE_HPP