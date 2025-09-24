#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

/*
    Coordenadas:

    Janela:
        (0, 0, -dJanela)

    Observador:
        (0, 0, 0)

    Esfera:
        (0, 0, -(dJanela + rEsfera))

*/
class Point {
    public:
        float x;
        float y;
        float z;
        int p;

        // Construtor
        Point(float x,float y, float z){
            this->x = x;
            this->y = y;
            this->z = z;
            this->p = 1;
        }
};

class Vector {
    public:
        float x;
        float y;
        float z;
        int v;       

        // Construtor
        Vector(float x, float y, float z) {
            this->x = x;
            this->y = y;
            this->z = z;
            this->v = 0;
        };
};


float calculaNorma(Vector& v){  
 return sqrt(pow(v.x,2) + pow(v.y,2) + pow(v.z,2));
}

float calculaProdutoEscalar(Vector& v, Vector& p){
    return (v.x * p.x) + (v.y * p.y) + (v.z * p.z);
}



int main(){

    Point olhoPintor(0, 0, 0);
    
    float dJanela = 20.0;
    float wJanela = 4.0;
    float hJanela = 4.0;
    Point centroJanela(0, 0, -dJanela);

    vector<int> esfColor = {255, 0, 0};
    float rEsfera = 1.0;
    Point centroEsfera(0, 0, -(dJanela + rEsfera));

    vector<int> bgColor = {100, 100, 100};

    int nCol = 400;
    int nLin = 400;
                        
    FILE *fp = fopen("tela.ppm", "wb");
    
    fprintf(fp, "P6\n%d %d\n255\n", nCol, nLin);

    float Dx = wJanela/nCol;
    float Dy = hJanela/nLin;

    float z = -dJanela;
    for (int l = 0; l < nLin; l++) {
        float y = (hJanela/2) - (Dy/2) - l*Dy;

        for (int c = 0; c < nCol; c++) {
            float x = - (wJanela/2) + (Dx/2) + c*Dx;

            Point pj(x, y, -dJanela);
            
            Vector Dr(pj.x - olhoPintor.x, pj.y - olhoPintor.y, pj.z - olhoPintor.z);

            float normaDr = calculaNorma(Dr);

            Vector dr((Dr.x/normaDr), (Dr.y/normaDr), (Dr.z/normaDr));
            
            Vector w(olhoPintor.x - centroEsfera.x, olhoPintor.y - centroEsfera.y, olhoPintor.z - centroEsfera.z);

            float a = calculaProdutoEscalar(dr, dr); 

            float bDelta = 2 * calculaProdutoEscalar(w, dr);

            float cDelta = calculaProdutoEscalar(w, w) - pow(rEsfera, 2) ;

            float delta = pow(bDelta,2) - (4 * a * cDelta);

            static unsigned char color[3];
            if(delta >= 0) {
                color[0] = esfColor[0];  
                color[1] = esfColor[1];  
                color[2] = esfColor[2];  
            } else {
                color[0] = bgColor[0];   
                color[1] = bgColor[1];   
                color[2] = bgColor[2];   
            }

            fwrite(color, 1, 3, fp);
        }
    }

    fclose(fp);
    
    return 0;
}

