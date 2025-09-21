#include <iostream>

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
/*
float calcularIntersecaoRaio() {                            - Talvez n precisa, pois como eh 2D 
                                                              o rEsfera eh suficiente para saber 
}                                                             se esta dentro ou fora da area da esfera.
*/
int main(){

    Point olhoPintor(0, 0, 0);
    
    float dJanela = 20.0;
    float wJanela = 4.0;
    float hJanela = 4.0;
    Point centroJanela(0, 0, -dJanela);

    float esfColor[3] = {255, 0, 0};
    float rEsfera = 1.0;
    Point centroEsfera(0, 0, -(dJanela + rEsfera));

    float bgColor[3] = {100, 100, 100};

    int nCol = 100;
    int nLin = 100;

    float Dx = wJanela/nCol;
    float Dy = hJanela/nLin;

    for (int l = 0; l < nLin; l++) {
        float y = (hJanela/2) - (Dy/2) - l*Dy;

        for (int c = 0; c < nCol; c++) {
            float x = - (wJanela/2) + (Dx/2) + c*Dx;

            Point P(x, y, -dJanela);

            // Cast um raio de olhoPintor para P - talvez n (?)

            // calcura distancia entre P e centroEsfera

            // Se distancia P-centroEsfera > rEsfera => bgColor
            // Senao => esfColor
        }
    }

    return 0;
}

