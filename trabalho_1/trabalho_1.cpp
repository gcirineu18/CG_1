#include <iostream>
#include <fstream>
#include <cmath>

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
        Point(float x, float y, float z)
        : x(x), y(y), z(z), p(1) {}
};

class Vector {
    public:
        float x;
        float y;
        float z;
        int v;       

        // Construtor
        Vector(float x, float y, float z)
        : x(x), y(y), z(z), v(1) {}
};

class Color {
public:
    int r;
    int g;
    int b;

    // Construtor
    Color(int r, int g, int b)
        : r(r), g(g), b(b) {}
};

class Janela {
public:
    float largura;
    float altura;
    Point centro;

    // Construtor
    Janela(float largura, float altura, Point centro)
        : largura(largura), altura(altura), centro(centro) {}
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

    Color esfColor(255, 0, 0);
    float rEsfera = 1.0;
    Point centroEsfera(0, 0, -(dJanela + rEsfera));

    Color bgColor(100, 100, 100);

    int nCol = 400;
    int nLin = 400;
                        
    ofstream fp("tela.ppm", ios::binary);

    fp << "P6\n" << nCol << " " << nLin << "\n255\n";    // Gera cabeçalho PPM 

    // Definir deltaX e deltaY - Tamanho de 1 quadrante do Canvas
    float Dx = wJanela/nCol;
    float Dy = hJanela/nLin;

    float z = -dJanela;
    for (int l = 0; l < nLin; l++) {
        float y = (hJanela/2) - (Dy/2) - l*Dy;

        for (int c = 0; c < nCol; c++) {
            float x = - (wJanela/2) + (Dx/2) + c*Dx;

            Point pj(x, y, -dJanela);    // Ponto correspondente ao centro do quadrante
            
            Vector Dr(pj.x - olhoPintor.x, pj.y - olhoPintor.y, pj.z - olhoPintor.z); // Vetor que corresponde ao raio que saio do olho e passa no centro do quadrante

            float normaDr = calculaNorma(Dr);

            Vector dr((Dr.x/normaDr), (Dr.y/normaDr), (Dr.z/normaDr)); // Vetor Unitário
            
            // De acordo com o processo matemático calculado, criamos o vetor w, que corresponde ao ponto olhoPintor - centroEsfera
            // P(t) = Po + t*dr -> Pj = Po + t*dr -> Pj - C = Po + t*dr - C => Po - C + t*dr => w + t*dr
            
            Vector w(olhoPintor.x - centroEsfera.x, olhoPintor.y - centroEsfera.y, olhoPintor.z - centroEsfera.z);

            // Com o w podemos calcular os componentes a, b e c da nossa equação do segundo grau
            // No entando, por agora apenas precisamos saber se o discriminante é maior ou menor que 0
            // Quando formos implementar a iluminação, iremos precisar dos valores de t
            
            float a = calculaProdutoEscalar(dr, dr); 

            float bDelta = 2 * calculaProdutoEscalar(w, dr);

            float cDelta = calculaProdutoEscalar(w, w) - pow(rEsfera, 2) ;

            float delta = pow(bDelta,2) - (4 * a * cDelta);

            // Se o delta é positivo ou 0, encontramos a esfera, se não pintamos com a cor de background
            if(delta >= 0) {
                fp.put(esfColor.r);
                fp.put(esfColor.g);
                fp.put(esfColor.b);  
            } else {
                fp.put(bgColor.r);
                fp.put(bgColor.r);
                fp.put(bgColor.r);   
            }
        }
    }

    fp.close();
    cout << "Imagem gerada." << endl;
    return 0;
}

