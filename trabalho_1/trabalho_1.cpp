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

class Ray{
    public:
       Point p0; 
       int t;
       
    //    Ray(Point p0, int t){
    //     this->p0 = new Point(p0.x, p0.y, p0.z);
    //     this->t = t;
    //    };
       
};

/*
float calcularIntersecaoRaio() {                            - Talvez n precisa, pois como eh 2D 
                                                              o rEsfera eh suficiente para saber 
}                                                             se esta dentro ou fora da area da esfera.
*/

int calculaNorma(Vector& v){  
 return sqrt(pow(v.x,2) + pow(v.y,2) + pow(v.z,2));
}

int calculaProdutoEscalar(Vector& v, Vector& p){
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

    int nCol = 100;
    int nLin = 100;

    //int matCanvas[100][100];

    vector<vector<vector<int>>> matCanvas(nLin, vector<vector<int>>(nCol, vector<int>(3, 0)));

    float Dx = wJanela/nCol;
    float Dy = hJanela/nLin;

    float z = -dJanela;
    for (int l = 0; l < nLin; l++) {
        float y = (hJanela/2) - (Dy/2) - l*Dy;

        for (int c = 0; c < nCol; c++) {
            float x = - (wJanela/2) + (Dx/2) + c*Dx;

            Point pj(x, y, -dJanela);
            
            Vector Dr(pj.x - olhoPintor.x, pj.y - olhoPintor.y, pj.z - olhoPintor.z);

            int normaDr = calculaNorma(Dr);

            Vector dr((Dr.x/normaDr), (Dr.y/normaDr), (Dr.z/normaDr));
            
            Vector w(olhoPintor.x - centroEsfera.x, olhoPintor.y - centroEsfera.y, olhoPintor.z - centroEsfera.z);

            int bDelta = 2 * calculaProdutoEscalar(w, dr);

            int cDelta = pow(rEsfera, 2) * calculaProdutoEscalar(w, w);

            int delta = pow(bDelta,2) + 4 * cDelta;

            if(delta >= 0) {
                matCanvas[l][c][0] = esfColor[0];  
                matCanvas[l][c][1] = esfColor[1];  
                matCanvas[l][c][2] = esfColor[2];  
            } else {
                matCanvas[l][c][0] = bgColor[0];   
                matCanvas[l][c][1] = bgColor[1];   
                matCanvas[l][c][2] = bgColor[2];   
            }
            // Cast um raio de olhoPintor para P - talvez n (?)

            // calcura distancia entre P e centroEsfera

            // Se distancia P-centroEsfera > rEsfera => bgColor
            // Senao => esfColor
        }
    }

    return 0;
}

