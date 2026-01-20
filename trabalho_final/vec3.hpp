#ifndef VEC3_HPP
#define VEC3_HPP

#include <cmath>

struct Vec3 {
    double x, y, z;

    Vec3(double x=0, double y=0, double z=0): x(x),y(y),z(z) {}

    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }

    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }

    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

    Vec3 operator*(const Vec3& v) const { return {x * v.x, y * v.y, z * v.z}; }

    //Produto escalar
    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

    //Produto vetorial
    Vec3 cross(const Vec3& v) const {
        return {y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x};
    }

    double length() const { return sqrt(dot(*this)); }

    Vec3 normalize() const { return (*this) * (1.0 / length()); }

    Vec3 rotX(double a) const{
        double matrix[4][4] = {
            {1, 0, 0, 0},
            {0, cos(a), -sin(a), 0},
            {0, sin(a), cos(a), 0},
            {0, 0, 0, 1}
        };

        double nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        double ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        double nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    Vec3 rotY(double a) const {
        double matrix[4][4] = {
            {cos(a), 0, sin(a), 0}, 
            {0, 1, 0, 0}, 
            {-sin(a), 0, cos(a), 0}, 
            {0, 0, 0, 1}
        };

        double nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        double ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        double nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    Vec3 rotZ(double a) const {
        double matrix[4][4] = {
            {cos(a), -sin(a), 0, 0}, 
            {sin(a), cos(a), 0, 0}, 
            {0, 0, 1, 0}, 
            {0, 0, 0, 1}
        };

        double nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        double ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        double nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    // ----- TRANSLAÇÃO -----
    Vec3 translate(double tx, double ty, double tz) const {
        double matrix[4][4] = {
            {1, 0, 0, tx}, 
            {0, 1, 0, ty}, 
            {0, 0, 1, tz}, 
            {0, 0, 0, 1}
        };

        double nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        double ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        double nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    Vec3 translate(const Vec3& t) const {
        return translate(t.x, t.y, t.z);
    }

    // ----- ESCALA -----
    Vec3 scale(double sx, double sy, double sz) const {
        double matrix[4][4] = {
            {sx, 0, 0, 0}, 
            {0, sy, 0, 0}, 
            {0, 0, sz, 0}, 
            {0, 0, 0, 1}
        };

        double nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        double ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        double nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    Vec3 scale(double s) const {
        return scale(s, s, s);
    }

    // ----- CISALHAMENTO -----
    Vec3 shear(double xy, double xz, double yx, double yz, double zx, double zy) const {
        double matrix[4][4] = {
            {1, xy, xz, 0}, 
            {yx, 1, yz, 0}, 
            {zx, zy, 1, 0}, 
            {0, 0, 0, 1}
        };

        double nx = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z) + (matrix[0][3] * 1);
        double ny = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z) + (matrix[1][3] * 1);
        double nz = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z) + (matrix[2][3] * 1);

        return Vec3(nx, ny, nz);
    }

    
};

#endif
