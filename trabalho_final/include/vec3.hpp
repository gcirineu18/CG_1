#ifndef VEC3_HPP
#define VEC3_HPP

#include <cmath>

struct Vec3 {
    float x, y, z;

    Vec3(float x=0, float y=0, float z=0): x(x),y(y),z(z) {}

    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }

    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }

    Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }

    Vec3 operator*(const Vec3& v) const { return {x * v.x, y * v.y, z * v.z}; }
    
    Vec3 operator-() const {
        return Vec3(-x, -y, -z);
    }
    //Produto escalar
    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

    //Produto vetorial
    Vec3 cross(const Vec3& v) const {
        return {y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x};
    }

    float length() const { return sqrt(dot(*this)); }

    Vec3 normalize() const { return (*this) * (1.0 / length()); }

    Vec3 rotX(float a) const{
        float matrix[4][4] = {
            {1, 0, 0, 0},
            {0, cos(a), -sin(a), 0},
            {0, sin(a), cos(a), 0},
            {0, 0, 0, 1}
        };

        float nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        float ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        float nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    Vec3 rotY(float a) const {
        float matrix[4][4] = {
            {cos(a), 0, sin(a), 0}, 
            {0, 1, 0, 0}, 
            {-sin(a), 0, cos(a), 0}, 
            {0, 0, 0, 1}
        };

        float nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        float ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        float nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    Vec3 rotZ(float a) const {
        float matrix[4][4] = {
            {cos(a), -sin(a), 0, 0}, 
            {sin(a), cos(a), 0, 0}, 
            {0, 0, 1, 0}, 
            {0, 0, 0, 1}
        };

        float nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        float ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        float nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    // Rotação Arbitrária
    // Método de Mudança de Sistemas de Coordenadas
    Vec3 rotArb(float a, Vec3 axis) const {
        Vec3 w = axis.normalize();

        Vec3 t = (std::abs(w.x) > 0.9f) ? Vec3(0, 1, 0) : Vec3(1, 0, 0);
        Vec3 u = t.cross(w).normalize();
        Vec3 v = w.cross(u);

        float B[4][4] = {
            {u.x, u.y, u.z, 0},
            {v.x, v.y, v.z, 0},
            {w.x, w.y, w.z, 0},
            {  0,   0,   0, 1}
        };

        float px = (B[0][0] * x) + (B[0][1] * y) + (B[0][2] * z);
        float py = (B[1][0] * x) + (B[1][1] * y) + (B[1][2] * z);
        float pz = (B[2][0] * x) + (B[2][1] * y) + (B[2][2] * z);
        Vec3 p_prime(px, py, pz);

        Vec3 p_rot = p_prime.rotZ(a);

        float nx = (B[0][0] * p_rot.x) + (B[1][0] * p_rot.y) + (B[2][0] * p_rot.z);
        float ny = (B[0][1] * p_rot.x) + (B[1][1] * p_rot.y) + (B[2][1] * p_rot.z);
        float nz = (B[0][2] * p_rot.x) + (B[1][2] * p_rot.y) + (B[2][2] * p_rot.z);

        return Vec3(nx, ny, nz);
    }

    // ----- TRANSLAÇÃO -----
    Vec3 translate(float tx, float ty, float tz) const {
        float matrix[4][4] = {
            {1, 0, 0, tx}, 
            {0, 1, 0, ty}, 
            {0, 0, 1, tz}, 
            {0, 0, 0, 1}
        };

        float nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        float ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        float nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    Vec3 translate(const Vec3& t) const {
        return translate(t.x, t.y, t.z);
    }

    // ----- ESCALA -----
    Vec3 scale(float sx, float sy, float sz) const {
        float matrix[4][4] = {
            {sx, 0, 0, 0}, 
            {0, sy, 0, 0}, 
            {0, 0, sz, 0}, 
            {0, 0, 0, 1}
        };

        float nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        float ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        float nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    Vec3 scale(float s) const {
        return scale(s, s, s);
    }

    // ----- CISALHAMENTO -----
    Vec3 shear(float xy, float xz, float yx, float yz, float zx, float zy) const {
        float matrix[4][4] = {
            {1, xy, xz, 0}, 
            {yx, 1, yz, 0}, 
            {zx, zy, 1, 0}, 
            {0, 0, 0, 1}
        };

        float nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        float ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        float nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    
};

#endif
