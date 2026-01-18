#ifndef CUBE_HPP
#define CUBE_HPP

#include "object.hpp"
#include "triangle.hpp"
#include <vector>

class Cube : public Object {
public:
    Vec3 min_corner;
    Vec3 max_corner; 
    Vec3 center;
    float size;
    std::vector<Triangle> triangles;

    // Construtor: posição e tamanho do cubo
    Cube(Vec3 center, float size, Material m, int _id) {
        mat = m;
        id = _id;

        this->center = center;
        this->size = size;

        float half = size / 3.0f;
        min_corner = center - Vec3(half, half, half);
        max_corner = center + Vec3(half, half, half);

        createMesh(m, _id);
    }

    void createMesh(Material m, int _id) {

        triangles.clear();

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

    void rotateX(double angle){
        float half = size / 3.0f;

        // Calcular os 8 vértices relativos ao centro
        Vec3 v0 = Vec3(-half, -half, -half).rotX(angle);
        Vec3 v1 = Vec3( half, -half, -half).rotX(angle);
        Vec3 v2 = Vec3( half,  half, -half).rotX(angle);
        Vec3 v3 = Vec3(-half,  half, -half).rotX(angle);
        Vec3 v4 = Vec3(-half, -half,  half).rotX(angle);
        Vec3 v5 = Vec3( half, -half,  half).rotX(angle);
        Vec3 v6 = Vec3( half,  half,  half).rotX(angle);
        Vec3 v7 = Vec3(-half,  half,  half).rotX(angle);

        updateMeshWithVertices(v0, v1, v2, v3, v4, v5, v6, v7);
    }

    void rotateY(double angle) {
        float half = size / 3.0f;
        
        Vec3 v0 = Vec3(-half, -half, -half).rotY(angle);
        Vec3 v1 = Vec3( half, -half, -half).rotY(angle);
        Vec3 v2 = Vec3( half,  half, -half).rotY(angle);
        Vec3 v3 = Vec3(-half,  half, -half).rotY(angle);
        Vec3 v4 = Vec3(-half, -half,  half).rotY(angle);
        Vec3 v5 = Vec3( half, -half,  half).rotY(angle);
        Vec3 v6 = Vec3( half,  half,  half).rotY(angle);
        Vec3 v7 = Vec3(-half,  half,  half).rotY(angle);

        updateMeshWithVertices(v0, v1, v2, v3, v4, v5, v6, v7);
    }
    
    void rotateZ(double angle) {
        float half = size / 3.0f;
        
        Vec3 v0 = Vec3(-half, -half, -half).rotZ(angle);
        Vec3 v1 = Vec3( half, -half, -half).rotZ(angle);
        Vec3 v2 = Vec3( half,  half, -half).rotZ(angle);
        Vec3 v3 = Vec3(-half,  half, -half).rotZ(angle);
        Vec3 v4 = Vec3(-half, -half,  half).rotZ(angle);
        Vec3 v5 = Vec3( half, -half,  half).rotZ(angle);
        Vec3 v6 = Vec3( half,  half,  half).rotZ(angle);
        Vec3 v7 = Vec3(-half,  half,  half).rotZ(angle);

        updateMeshWithVertices(v0, v1, v2, v3, v4, v5, v6, v7);
    }

    // Translação
    void translate(double tx, double ty, double tz) {
        center = center.translate(tx, ty, tz);
        min_corner = min_corner.translate(tx, ty, tz);
        max_corner = max_corner.translate(tx, ty, tz);
        createMesh(mat, id);
    }

    void translate(const Vec3& t) {
        translate(t.x, t.y, t.z);
    }

    // Escala
    void scaleTransform(double sx, double sy, double sz) {
        float half = size / 3.0f;
        
        Vec3 v0 = Vec3(-half, -half, -half).scale(sx, sy, sz);
        Vec3 v1 = Vec3( half, -half, -half).scale(sx, sy, sz);
        Vec3 v2 = Vec3( half,  half, -half).scale(sx, sy, sz);
        Vec3 v3 = Vec3(-half,  half, -half).scale(sx, sy, sz);
        Vec3 v4 = Vec3(-half, -half,  half).scale(sx, sy, sz);
        Vec3 v5 = Vec3( half, -half,  half).scale(sx, sy, sz);
        Vec3 v6 = Vec3( half,  half,  half).scale(sx, sy, sz);
        Vec3 v7 = Vec3(-half,  half,  half).scale(sx, sy, sz);

        updateMeshWithVertices(v0, v1, v2, v3, v4, v5, v6, v7);
    }

    void scaleTransform(double s) {
        scaleTransform(s, s, s);
    }

    // Cisalhamento
    void shearTransform(double xy, double xz, double yx, double yz, double zx, double zy) {
        float half = size / 3.0f;
        
        Vec3 v0 = Vec3(-half, -half, -half).shear(xy, xz, yx, yz, zx, zy);
        Vec3 v1 = Vec3( half, -half, -half).shear(xy, xz, yx, yz, zx, zy);
        Vec3 v2 = Vec3( half,  half, -half).shear(xy, xz, yx, yz, zx, zy);
        Vec3 v3 = Vec3(-half,  half, -half).shear(xy, xz, yx, yz, zx, zy);
        Vec3 v4 = Vec3(-half, -half,  half).shear(xy, xz, yx, yz, zx, zy);
        Vec3 v5 = Vec3( half, -half,  half).shear(xy, xz, yx, yz, zx, zy);
        Vec3 v6 = Vec3( half,  half,  half).shear(xy, xz, yx, yz, zx, zy);
        Vec3 v7 = Vec3(-half,  half,  half).shear(xy, xz, yx, yz, zx, zy);

        updateMeshWithVertices(v0, v1, v2, v3, v4, v5, v6, v7);
    }

    private:
        void updateMeshWithVertices(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3, 
            Vec3 v4, Vec3 v5, Vec3 v6, Vec3 v7)
        {
            triangles.clear();
            
            v0 = v0 + center;
            v1 = v1 + center;
            v2 = v2 + center;
            v3 = v3 + center;
            v4 = v4 + center;
            v5 = v5 + center;
            v6 = v6 + center;
            v7 = v7 + center;

            triangles.push_back(Triangle(v0, v1, v2, mat, id));
            triangles.push_back(Triangle(v0, v2, v3, mat, id));
            triangles.push_back(Triangle(v5, v4, v7, mat, id));
            triangles.push_back(Triangle(v5, v7, v6, mat, id));
            triangles.push_back(Triangle(v4, v0, v3, mat, id));
            triangles.push_back(Triangle(v4, v3, v7, mat, id));
            triangles.push_back(Triangle(v1, v5, v6, mat, id));
            triangles.push_back(Triangle(v1, v6, v2, mat, id));
            triangles.push_back(Triangle(v4, v5, v1, mat, id));
            triangles.push_back(Triangle(v4, v1, v0, mat, id));
            triangles.push_back(Triangle(v3, v2, v6, mat, id));
            triangles.push_back(Triangle(v3, v6, v7, mat, id));

        }

};

#endif // CUBE_HPP