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

    Vec3 vertices[8];

    Cube(Vec3 center, float size, Material m, const char* nome) {
        mat = m;
        name = nome;

        this->center = center;
        this->size = size;

        float half = size / 3.0f;
        min_corner = center - Vec3(half, half, half);
        max_corner = center + Vec3(half, half, half);

        initializeVertices();
        createMesh(m);
    }

    // Inicializar vértices no estado original
    void initializeVertices() {
        float half = size / 3.0f;
        vertices[0] = Vec3(-half, -half, -half);
        vertices[1] = Vec3( half, -half, -half);
        vertices[2] = Vec3( half,  half, -half);
        vertices[3] = Vec3(-half,  half, -half);
        vertices[4] = Vec3(-half, -half,  half);
        vertices[5] = Vec3( half, -half,  half);
        vertices[6] = Vec3( half,  half,  half);
        vertices[7] = Vec3(-half,  half,  half);
    }

    void createMesh(Material m) {
        triangles.clear();

        // Usar vértices armazenados + centro
        Vec3 v0 = vertices[0] + center;
        Vec3 v1 = vertices[1] + center;
        Vec3 v2 = vertices[2] + center;
        Vec3 v3 = vertices[3] + center;
        Vec3 v4 = vertices[4] + center;
        Vec3 v5 = vertices[5] + center;
        Vec3 v6 = vertices[6] + center;
        Vec3 v7 = vertices[7] + center;

        triangles.push_back(Triangle(v0, v1, v2, m));
        triangles.push_back(Triangle(v0, v2, v3, m));
        triangles.push_back(Triangle(v6, v7, v5, m));
        triangles.push_back(Triangle(v4, v5, v7, m));
        triangles.push_back(Triangle(v4, v0, v3, m));
        triangles.push_back(Triangle(v4, v3, v7, m));
        triangles.push_back(Triangle(v1, v5, v6, m));
        triangles.push_back(Triangle(v1, v6, v2, m));
        triangles.push_back(Triangle(v4, v5, v1, m));
        triangles.push_back(Triangle(v4, v1, v0, m));
        triangles.push_back(Triangle(v3, v2, v6, m));
        triangles.push_back(Triangle(v3, v6, v7, m));
    }

    void rotateY(double angle) {
        for (int i = 0; i < 8; i++) {
            vertices[i] = vertices[i].rotY(angle);
        }
        createMesh(mat);
    }

    void rotateX(double angle) {
        for (int i = 0; i < 8; i++) {
            vertices[i] = vertices[i].rotX(angle);
        }
        createMesh(mat);
    }

    void rotateZ(double angle) {
        for (int i = 0; i < 8; i++) {
            vertices[i] = vertices[i].rotZ(angle);
        }
        createMesh(mat);
    }

    void scaleTransform(double sx, double sy, double sz) {
        for (int i = 0; i < 8; i++) {
            vertices[i] = vertices[i].scale(sx, sy, sz);
        }
        createMesh(mat);
    }

    void scaleTransform(double s) {
        scaleTransform(s, s, s);
    }

    void shearTransform(double xy, double xz, double yx, double yz, double zx, double zy) {
        for (int i = 0; i < 8; i++) {
            vertices[i] = vertices[i].shear(xy, xz, yx, yz, zx, zy);
        }
        createMesh(mat);
    }

    void translate(double tx, double ty, double tz) {
        center = center.translate(tx, ty, tz);
        min_corner = min_corner.translate(tx, ty, tz);
        max_corner = max_corner.translate(tx, ty, tz);
        createMesh(mat);
    }

    void translate(const Vec3& t) {
        translate(t.x, t.y, t.z);
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