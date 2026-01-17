#ifndef VERTEX_HPP
#define VERTEX_HPP

#include "object.hpp"

class Vertex : public Object {
    private:
        Vec3 vertex;

    public:
         Vertex(double x, double y, double z){
            vertex.x = x;
            vertex.y = y;
            vertex.z = z;
         }    



};

#endif //VERTEX_HPP
