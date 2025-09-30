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

Vector produtoVetorEscalar(float a, Vector& v) {
    Vector p(v.x*a, v.y*a, v.z*a);
    return p;
}

Vector somaPontoVetor(Point& v, Vector& p) {
    Vector s(v.x + p.x, v.y + p.y, v.z + p.z);
    return s;
}

Vector subtracaoVetorPonto(Vector& v, Point& p) {
    Vector d(v.x - p.x, v.y - p.y, v.z - p.z);
    return d;
}

Vector subtracaoPontoVetor(Point& p, Vector& v) {
    Vector d(p.x - v.x, p.y - v.y, p.z - v.z);
    return d;
}

Vector divisaoVetorEscalar(Vector& v, float a) {
    Vector d(v.x/a, v.y/a, v.z/a);
    return d;
}

float minimo(float x, float y) {
    if (x <= y) {
        return x;
    } else {
        return y;
    }
}

Vector arroba(Vector& v, Vector& p) {
    Vector w((v.x * p.x),(v.y * p.y),(v.z * p.z));
    return w;
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

    Vector i_f(0.7, 0.7, 0.7);
    Point p_f(0, 5, 0);

    Vector K(1, 0, 0);
    // Vamos tratar i_f e K como Vetores para fazermos os cálculos e após tudo vamo transformá-los para Color

    float m = 10.0;

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

            float aDelta = calculaProdutoEscalar(dr, dr);

            float bDelta = 2 * calculaProdutoEscalar(w, dr);

            float cDelta = calculaProdutoEscalar(w, w) - pow(rEsfera, 2) ;

            float delta = pow(bDelta,2) - (4 * aDelta * cDelta);

            // Se o delta é positivo ou 0, encontramos a esfera, se não pintamos com a cor de background
            if(delta < 0) {
                fp.put(bgColor.r);
                fp.put(bgColor.r);
                fp.put(bgColor.r);
            } else {
                // Calcular t
                float t1 = (- bDelta + sqrt(delta))/2 * aDelta;
                float t2 = (- bDelta - sqrt(delta))/2 * aDelta;

                float t = minimo(t1, t2);

                // Calcular vetores n, l, v e r

                // Calcula P(t) de acordo com a equação do Ray
                Vector tdr = produtoVetorEscalar(t, dr);
                Vector Pi = somaPontoVetor(olhoPintor, tdr);

                // Calcula n
                Vector PiC = subtracaoVetorPonto(Pi, centroEsfera);
                Vector n = divisaoVetorEscalar(PiC,rEsfera);

                // Calcula l
                Vector pfPI = subtracaoPontoVetor(p_f, Pi);
                float norm_pfPI = calculaNorma(pfPI);
                Vector l = divisaoVetorEscalar(pfPI, norm_pfPI);

                // Calcula x
                Vector x(-dr.x, -dr.y, -dr.z);

                // Calcula r
                float ln = calculaProdutoEscalar(l, n);
                Vector lnn = produtoVetorEscalar(ln ,n);
                Vector r(lnn.x - l.x, lnn.y- l.y, lnn.z - l.z);

                // Calcular Id e Ie
                Vector Ifk = arroba(i_f, K);
                Vector Id = produtoVetorEscalar(ln, Ifk);

                float rx = calculaProdutoEscalar(r, x);
                float rxm = pow(rx, m);
                Vector Ie = produtoVetorEscalar(rxm, Ifk);

                Vector I_E(Ie.x + Id.x, Ie.y + Id.y, Ie.z + Id.z);

                // Imprimir as cores no arquivo
                fp.put(I_E.x * 255);
                fp.put(I_E.y * 255);
                fp.put(I_E.z * 255);
            }
        }
    }

    fp.close();
    cout << "Imagem gerada." << endl;
    return 0;
}
