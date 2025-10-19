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

vector<float> prod_vetorial(vector<float>& v, vector<float>& p){
    
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
    float dJanela = 30.0f;
    float wJanela = 60.0f;
    float hJanela = 60.0f;
    vector<float> centroJanela = criarPonto(0.0f, 0.0f, -dJanela);

    // Definir características da Tela de Mosquito
    int nCol = 500;
    int nLin = 500;
    float Dx = wJanela/nCol;    /* Medidas de  1 quadrante da tela de mosquito */
    float Dy = hJanela/nLin;

    // Definir cor do Background - Cinza
    vector<float> bgColor = criarVetor(100,100,100);

    // Definir características da Fonte Pontual
    vector<float> fonte_p_int = criarVetor(0.7f, 0.7f, 0.7f);
    vector<float> fonte_p_coord = criarPonto(0.0f, 60.0f, -30.0f);

    // Definir características da Fonte Ambiente
    vector<float> I_A = criarVetor(0.3f, 0.3f, 0.3f);

    /* ============= OBJETOS DA CENA ============= */
    // Iniciar Esfera
    float rEsfera = 40.0f;
    vector<float> esfColor = criarVetor(255, 0, 0);
    vector<float> centroEsfera = criarPonto(0.0f, 0.0f, -100.0f);
    vector<float> K_esf = criarVetor(0.7f, 0.2f, 0.2f);
    float m_esf = 10.0f;

    // Plano de Chão
    vector<float> P_pi_c = criarPonto(0.0f, -rEsfera, 0.0f);
    vector<float> n_bar_c = criarVetor(0.0f, 1.0f, 0.0f);
    vector<float> Kd_c = criarVetor(0.2f, 0.7f, 0.2f);
    vector<float> Ke_c = criarVetor(0.0f, 0.0f, 0.0f);
    float m_c = 1.0f;

    // Plano de Fundo
    vector<float> P_pi_f = criarPonto(0.0f, 0.0f, -200.0f);
    vector<float> n_bar_f = criarVetor(0.0f, 0.0f, 1.0f);
    vector<float> Kd_f = criarVetor(0.3f, 0.3f, 0.7f);
    vector<float> Ke_f = criarVetor(0.0f, 0.0f, 0.0f);
    float m_f = 1.0f;


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

            vector<float> w_esf = subt(olhoPintor, centroEsfera);

            // Com o w podemos calcular os componentes a, b e c da nossa equação do segundo grau
            float aDelta = prod_escalar(dr_u, dr_u);

            float bDelta = 2 * prod_escalar(w_esf, dr_u);

            float cDelta = prod_escalar(w_esf, w_esf) - pow(rEsfera, 2) ;

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
                vector<float> pfPI = subt(fonte_p_coord, Pi);
                float pfPI_norma = norma(pfPI);
                vector<float> l = div_escalar(pfPI, pfPI_norma);

                // Calcula x
                vector<float> x = {-dr_u[0], -dr_u[1], -dr_u[2]};

                // Calcula r
                float ln = prod_escalar(l, n);
                vector<float> r = {(2*ln*n[0]) - l[0], (2*ln*n[1]) - l[1], (2 * ln*n[2]) - l[2]};

                // Calcular Id e Ie e Iamb
                float ln_limitado = max(0.0f, prod_escalar(l, n));
                vector<float> Ifk = arroba(fonte_p_int, K_esf);
                vector<float> Id = mult_escalar(Ifk, ln_limitado);

                float rx = max(0.0f, prod_escalar(r, x));
                float rxm = pow(rx, m_esf);
                vector<float> Ie = mult_escalar(Ifk, rxm);

                vector<float> Iamb = arroba(I_A, K_esf);

                vector<float> I_E = {Ie[0] + Id[0] + Iamb[0], Ie[1] + Id[1] + Iamb[1], Ie[2] + Id[2] + Iamb[2]};

                // Imprimir as cores no arquivo
                if (I_E[0] > 1.0f){
                    fp.put(255);
                } else{
                    fp.put(I_E[0] * 255);
                }
                fp.put(I_E[1] * 255);
                fp.put(I_E[2] * 255);

            } else {
                /* ============= VERIFICAÇÃO PARA PLANO CHÃO ============= */
                vector<float> w_p_chao = subt(olhoPintor, P_pi_c);
                float wpcnc = - prod_escalar(w_p_chao, n_bar_c);
                float drnc = prod_escalar(dr_u, n_bar_c);
                float ti_p_c = wpcnc / drnc;

                /* ============= VERIFICAÇÃO PARA PLANO FUNDO ============= */
                vector<float> w_p_fundo = subt(olhoPintor, P_pi_f);
                float wpfnf = - prod_escalar(w_p_fundo, n_bar_f);
                float drnf = prod_escalar(dr_u, n_bar_f);
                float ti_p_f = wpfnf / drnf;

                if (ti_p_c <= ti_p_f && ti_p_c > 0) {
                    // Calcula P(t) de acordo com a equação do Ray
                    vector<float> tdr = mult_escalar(dr_u, ti_p_c);
                    vector<float> Pi = soma(olhoPintor, tdr);

                    // Calcula l
                    vector<float> pfPI = subt(fonte_p_coord, Pi);
                    float pfPI_norma = norma(pfPI);
                    vector<float> l = div_escalar(pfPI, pfPI_norma);

                    // P(s) = Pi + s*l
                    vector<float> s_esf = subt(Pi, centroEsfera);

                    // Com o w podemos calcular os componentes a, b e c da nossa equação do segundo grau
                    float aDelta = prod_escalar(l, l);

                    float bDelta = 2 * prod_escalar(s_esf, l);

                    float cDelta = prod_escalar(s_esf, s_esf) - pow(rEsfera, 2) ;

                    float delta = pow(bDelta,2) - (4 * aDelta * cDelta);

                    if (delta >= 0) {
                        l = criarVetor(0.0f, 0.0f, 0.0f);
                    }

                    // Calcula x
                    vector<float> x = {-dr_u[0], -dr_u[1], -dr_u[2]};

                    // Calcula r
                    float ln = prod_escalar(l, n_bar_c);
                    vector<float> r = {(2*ln*n_bar_c[0]) - l[0], (2*ln*n_bar_c[1]) - l[1], (2 * ln*n_bar_c[2]) - l[2]};

                    // Calcular Id e Ie e Iamb
                    float ln_limitado = max(0.0f, ln);
                    vector<float> Ifkd = arroba(fonte_p_int, Kd_c);
                    vector<float> Id = mult_escalar(Ifkd, ln_limitado);

                    float rx = max(0.0f, prod_escalar(r, x));
                    float rxm = pow(rx, m_c);
                    vector<float> Ifke = arroba(fonte_p_int, Ke_c);
                    vector<float> Ie = mult_escalar(Ifke, rxm);

                    vector<float> Iamb = arroba(I_A, Kd_c);

                    vector<float> I_E_c = {Ie[0] + Id[0] + Iamb[0], Ie[1] + Id[1] + Iamb[1], Ie[2] + Id[2] + Iamb[2]};
                    
                    // Imprimir as cores no arquivo
                    if (I_E_c[0] > 1.0f){
                        fp.put(255);
                    } else{
                        fp.put(I_E_c[0] * 255);
                    }
                    if (I_E_c[1] > 1.0f){
                        fp.put(255);
                    } else{
                        fp.put(I_E_c[1] * 255);
                    }
                    if (I_E_c[2] > 1.0f){
                        fp.put(255);
                    } else{
                        fp.put(I_E_c[2] * 255);
                    }

                } else if (ti_p_f > 0) {
                    // Calcula P(t) de acordo com a equação do Ray
                    vector<float> tdr = mult_escalar(dr_u, ti_p_f);
                    vector<float> Pi = soma(olhoPintor, tdr);

                    // Calcula l
                    vector<float> pfPI = subt(fonte_p_coord, Pi);
                    float pfPI_norma = norma(pfPI);
                    vector<float> l = div_escalar(pfPI, pfPI_norma);

                    // P(s) = Pi + s*l
                    vector<float> s_esf = subt(Pi, centroEsfera);

                    // Com o w podemos calcular os componentes a, b e c da nossa equação do segundo grau
                    float aDelta = prod_escalar(l, l);

                    float bDelta = 2 * prod_escalar(s_esf, l);

                    float cDelta = prod_escalar(s_esf, s_esf) - pow(rEsfera, 2) ;

                    float delta = pow(bDelta,2) - (4 * aDelta * cDelta);

                    if (delta >= 0) {
                        l = criarVetor(0.0f, 0.0f, 0.0f);
                    }

                    // Calcula x
                    vector<float> x = {-dr_u[0], -dr_u[1], -dr_u[2]};

                    // Calcula r
                    float ln = prod_escalar(l, n_bar_f);
                    vector<float> r = {(2*ln*n_bar_f[0]) - l[0], (2*ln*n_bar_f[1]) - l[1], (2 * ln*n_bar_c[2]) - l[2]};

                    // Calcular Id e Ie e Iamb
                    float ln_limitado = max(0.0f, ln);
                    vector<float> Ifkd = arroba(fonte_p_int, Kd_f);
                    vector<float> Id = mult_escalar(Ifkd, ln_limitado);

                    float rx = max(0.0f, prod_escalar(r, x));
                    float rxm = pow(rx, m_f);
                    vector<float> Ifke = arroba(fonte_p_int, Ke_f);
                    vector<float> Ie = mult_escalar(Ifke, rxm);

                    vector<float> Iamb = arroba(I_A, Kd_f);

                    vector<float> I_E_f = {Ie[0] + Id[0] + Iamb[0], Ie[1] + Id[1] + Iamb[1], Ie[2] + Id[2] + Iamb[2]};

                    if (I_E_f[0] > 1.0f){
                        fp.put(255);
                    } else{
                        fp.put(I_E_f[0] * 255);
                    }
                    if (I_E_f[1] > 1.0f){
                        fp.put(255);
                    } else{
                        fp.put(I_E_f[1] * 255);
                    }
                    if (I_E_f[2] > 1.0f){
                        fp.put(255);
                    } else{
                        fp.put(I_E_f[2] * 255);
                    }
                        
                } else {
                    fp.put(bgColor[0]);
                    fp.put(bgColor[1]);
                    fp.put(bgColor[2]);
                }
            }
        }
    }

    /* ============= MENSAGEM DE SUCESSO =============  */
    fp.close();
    cout << "Imagem gerada com sucesso!" << endl;
    return 0;
}