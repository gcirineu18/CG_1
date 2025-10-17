#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

/* ============= FUNÇÕES ============= */

vector<float> criarPonto( float x, float y, float z){
    vector<float> p = {x,y,z,1};
    return p;
}

vector<float> criarVetor(float x, float y, float z){
    vector<float> v {x,y,z,0};
    return v;
}

float norma(vector<float>& v){
    return sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2));
}

float prod_escalar(vector<float>& v, vector<float>& p){
    return (v[0] * p[0]) + (v[1] * p[1]) + (v[2] * p[2]);
}

vector<float> mult_escalar(vector<float>& v, float x) {
    vector<float> q;
    if (v[3] == 1)
        q = criarPonto(v[0] * x, v[1] * x, v[2] * x);
    else
        q = criarVetor(v[0] * x, v[1] * x, v[2] * x);

    return q;
}

vector<float> soma(vector<float>& v, vector<float>& p) {
    vector<float> s = criarPonto(v[0] + p[0], v[1] + p[1], v[2] + p[2]);
    return s;
}

vector<float> subt(vector<float>& v, vector<float>& p) {
    vector<float> s = criarVetor(v[0] - p[0], v[1] - p[1], v[2] - p[2]);
    return s;
}

vector<float> div_escalar(vector<float>& v, float x) {
    vector<float> d;
    if (v[3] == 1)
        d = criarPonto(v[0]/x, v[1]/x, v[2]/x);
    else
        d = criarVetor(v[0]/x, v[1]/x, v[2]/x);

    return d;
}

vector<float> arroba(vector<float>& v, vector<float>& p) {
    vector<float> arr = criarVetor(v[0] * p[0], v[1] * p[1], v[2] * p[2]);
    return arr;
}

float minimo(float x, float y) {
    if (x <= y) {
        return x;
    } else {
        return y;
    }
}

int main() {

    // Localização do olho do Pintor
    vector<float> olhoPintor = criarPonto(0.0f, 0.0f, 0.0f);

    // Definir características da Janela
    float dJanela = 4.0f;
    float wJanela = 4.0f;
    float hJanela = 4.0f;
    vector<float> centroJanela = criarPonto(0.0f, 0.0f, -dJanela);

    // Definir características da Tela de Mosquito
    int nCol = 400;
    int nLin = 400;
    float Dx = wJanela/nCol;    /* Medidas de  1 quadrante da tela de mosquito */
    float Dy = hJanela/nLin;

    // Definir cor do Background - Cinza
    vector<float> bgColor = criarPonto(100,100,100);

    // Definir características da Fonte Luminosa
    vector<float> fonte_int = criarPonto(0.7f, 0.7f, 0.7f);
    vector<float> fonte_coord = criarPonto(0.0f, 4.0f, 0.0f);
    vector<float> K = criarPonto(1.0f, 0.0f, 0.0f);

    // m
    float m = 10.0f;

    /* ============= OBJETOS DA CENA ============= */
    // Iniciar Esfera
    float rEsfera = 1.0f;
    vector<float> esfColor = criarPonto(255, 0, 0);
    vector<float> centroEsfera = criarPonto(0.0f, 0.0f, -(dJanela + rEsfera));

    /* ============= ARQUIVO PPM - Cabeçalho ============= */
    ofstream fp("tela.ppm", ios::binary);
    fp << "P6\n" << nCol << " " << nLin << "\n255\n";

    /* ============= PERCORRER QUADRANTES DA TELA DE MOSQUITO ============= */
    float z = -dJanela;
    for (int l = 0; l < nLin; l++) {
        float y = (hJanela/2) - (Dy/2) - l*Dy;

        for (int c = 0; c < nCol; c++) {
            float x = - (wJanela/2) + (Dx/2) + c*Dx;

            // Vetor com coordenadas do centro do quadrante
            vector<float> pj = criarPonto(x, y, z);

            // Vetor correspondente ao raio lançado do olho do pintor para o centro do quadrante
            vector<float> dr = subt(pj, olhoPintor);

            // Calculando o vetor unitário de dr
            float dr_norma = norma(dr);
            vector<float> dr_u = div_escalar(dr, dr_norma);

            /* ============= VERIFICAÇÃO PARA ESFERA ============= */
            // De acordo com o processo matemático calculado, criamos o vetor w, que corresponde ao ponto olhoPintor - centroEsfera
            // P(t) = Po + t*dr -> Pj = Po + t*dr -> Pj - C = Po + t*dr - C => Po - C + t*dr => w + t*dr

            vector<float> w = subt(olhoPintor, centroEsfera);

            // Com o w podemos calcular os componentes a, b e c da nossa equação do segundo grau
            // No entando, por agora apenas precisamos saber se o discriminante é maior ou menor que 0
            // Quando formos implementar a iluminação, iremos precisar dos valores de t

            float aDelta = prod_escalar(dr_u, dr_u);

            float bDelta = 2 * prod_escalar(w, dr_u);

            float cDelta = prod_escalar(w, w) - pow(rEsfera, 2) ;

            float delta = pow(bDelta,2) - (4 * aDelta * cDelta);

            // Se o delta é positivo ou 0, encontramos a esfera, se não vamos testar para outros objetos utilizados na cena
            if (delta >= 0) {
                // calcular menor valor de t
                float t1 = (- bDelta + sqrt(delta))/2 * aDelta;
                float t2 = (- bDelta - sqrt(delta))/2 * aDelta;

                float t = minimo(t1, t2);

                /* Calcular vetores n, l, v e r */

                // Calcula P(t) de acordo com a equação do Ray
                vector<float> tdr = mult_escalar(dr_u, t);
                vector<float> Pi = soma(olhoPintor, tdr);

                // Calcula n
                vector<float> PiC = subt(Pi, centroEsfera);
                vector<float> n = div_escalar(PiC, rEsfera);

                // Calcula l
                vector<float> pfPI = subt(fonte_coord, Pi);
                float pfPI_norma = norma(pfPI);
                vector<float> l = div_escalar(pfPI, pfPI_norma);

                // Calcula x
                vector<float> x = {-dr_u[0], -dr_u[1], -dr_u[2]};

                // Calcula r
                float ln = prod_escalar(l, n);
                vector<float> r = {(2*ln*n[0]) - l[0], (2*ln*n[1]) - l[1], (2 * ln*n[2]) - l[2]};

                // Calcular Id e Ie
                float ln_limitado = max(0.0f, prod_escalar(l, n));
                vector<float> Ifk = arroba(fonte_int, K);
                vector<float> Id = mult_escalar(Ifk, ln_limitado);

                float rx = max(0.0f, prod_escalar(r, x));
                float rxm = pow(rx, m);
                vector<float> Ie = mult_escalar(Ifk, rxm);

                vector<float> I_E = {Ie[0] + Id[0], Ie[1] + Id[1], Ie[2] + Id[2]};

                // Imprimir as cores no arquivo
                if (I_E[0] > 1.0f){
                    fp.put(255);
                } else{
                    fp.put(I_E[0] * 255);
                }
                fp.put(I_E[1] * 255);
                fp.put(I_E[2] * 255);
            }
            /* else {

                adicionar verificação para plano

            } */
            else {
                fp.put(bgColor[0]);
                fp.put(bgColor[1]);
                fp.put(bgColor[2]);
            }
        }
    }

    /* ============= MENSAGEM DE SUCESSO =============  */
    fp.close();
    cout << "Imagem gerada com sucesso!" << endl;
    return 0;
}