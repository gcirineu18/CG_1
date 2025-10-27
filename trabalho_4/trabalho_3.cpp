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

float calcula_ti_plano(vector<float> olhoPintor, vector<float> P_pi, vector<float> n_bar, vector<float> dr_u){
    vector<float> w_p = subt(olhoPintor, P_pi);
    float wpcn = - prod_escalar(w_p, n_bar);
    float drn = prod_escalar(dr_u, n_bar);
    return wpcn /drn;
}

vector<float> arroba(vector<float>& v, vector<float>& p) {
    vector<float> arr = criarVetor(v[0] * p[0], v[1] * p[1], v[2] * p[2]);
    return arr;
}

vector<float> calcular_intensidade_olho(vector<float> dr_u, vector<float> l, vector<float> n_bar, vector<float> fonte_p_int, vector<float> Kd, vector<float> Ke, float m, vector<float> I_A){
    // Calcula x
    vector<float> x = {-dr_u[0], -dr_u[1], -dr_u[2]};

    // Calcula r
    float ln = prod_escalar(l, n_bar);
    vector<float> r = {(2*ln*n_bar[0]) - l[0], (2*ln*n_bar[1]) - l[1], (2 * ln*n_bar[2]) - l[2]};

    // Calcular Id e Ie e Iamb
    float ln_limitado = max(0.0f, ln);
    vector<float> Ifkd = arroba(fonte_p_int, Kd);
    vector<float> Id = mult_escalar(Ifkd, ln_limitado);

    float rx = max(0.0f, prod_escalar(r, x));
    float rxm = pow(rx, m);
    vector<float> Ifke = arroba(fonte_p_int, Ke);
    vector<float> Ie = mult_escalar(Ifke, rxm);

    vector<float> Iamb = arroba(I_A, Kd);

    return {Ie[0] + Id[0] + Iamb[0], Ie[1] + Id[1] + Iamb[1], Ie[2] + Id[2] + Iamb[2]};
}

float calcula_delta_plano(vector<float>*l, vector<float> dr_u, float ti_p, vector<float> olhoPintor, vector<float> fonte_p_coord, vector<float> centroEsfera, float rEsfera){
    // Calcula P(t) de acordo com a equação do Ray
    vector<float> tdr = mult_escalar(dr_u, ti_p);
    vector<float> Pi = soma(olhoPintor, tdr);

    // Calcula l
    vector<float> pfPI = subt(fonte_p_coord, Pi);
    float pfPI_norma = norma(pfPI);
     *l = div_escalar(pfPI, pfPI_norma);

    // P(s) = Pi + s*l
    vector<float> s_esf = subt(Pi, centroEsfera);

    // Com o w podemos calcular os componentes a, b e c da nossa equação do segundo grau
    float aDelta = prod_escalar(*l, *l);

    float bDelta = 2 * prod_escalar(s_esf, *l);

    float cDelta = prod_escalar(s_esf, s_esf) - pow(rEsfera, 2) ;

    return pow(bDelta,2) - (4 * aDelta * cDelta);
}

float calcula_intersecao_cone(vector<float>& olhoPintor, vector<float>& dr_u, vector<float>& V, vector<float>& n, float H, float R, float cosTheta) {
    vector<float> v = subt(olhoPintor, V);  // Ponto de origem do raio - Vértice do cone

    float a = pow(prod_escalar(dr_u, n),2) - (prod_escalar(dr_u, dr_u)*cosTheta*cosTheta);
    float b = prod_escalar(v, dr_u)*cosTheta*cosTheta - (prod_escalar(v,n)*prod_escalar(dr_u, n));
    float c = pow(prod_escalar(v, n),2) - (prod_escalar(v, v)*cosTheta*cosTheta);

    float delta = b * b - a * c;

    if (delta < 0) {
        return -1;  // Não há interseção
    }

    // Calculando o menor valor de t (pode ser uma interseção no cone)
    float t1 = (-b + sqrt(delta)) / a;
    float t2 = (-b - sqrt(delta)) / a;

    return min(t1, t2);
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

    // Cilindro
    vector<float> CB_cilindro = centroEsfera;
    float h_cilindro = 3.0f*rEsfera;
    vector<float> d_cil = criarVetor(-1/sqrt(3), 1/sqrt(3), -1/sqrt(3));
    float d_cil_norma = norma(d_cil);
    vector<float> d_cil_u = div_escalar(d_cil, d_cil_norma);
    vector<float> hd_cilindro = mult_escalar(d_cil_u, h_cilindro);
    vector<float> CT_cilindro = soma(CB_cilindro, hd_cilindro);

    cout << CB_cilindro[0] << " " << CB_cilindro[1] << " " << CB_cilindro[2] << endl;
    cout << CT_cilindro[0] << " " << CT_cilindro[1] << " " << CT_cilindro[2] << endl;

    // Cone
    vector<float> CB_cone = CT_cilindro;
    float rB_cone = 1.5f * rEsfera;
    float h_cone = rB_cone/3.0f;
    vector<float> d_cone = d_cil_u;
    vector<float> hd_cone = mult_escalar(d_cone, h_cone);
    vector<float> V_cone = soma(CB_cone, hd_cone);
    vector<float> Kd_cone = criarVetor(0.8f, 0.3f, 0.2f);
    vector<float> Ke_cone = criarVetor(0.8f, 0.3f, 0.2f);
    vector<float> Ka_cone = criarVetor(0.8f, 0.3f, 0.2f);
    float m_cone = 1.0f;

    // Nos cálculos de interseção vamos precisar do valor de Cos(T) em que T é o angulo que a geratriz faz com o eixo do cone. Vamos definir essas variáveis também para facilitar cálculos seguintes.
    // Para calcular o vetor geratriz, vamos fazer a diferença entre um ponto na circunferência da base do cone e o vértice do cone.
    // Para calcular o ponto na circunferencia da base do cone, vamos pegar o ponto do centro da base do cone e somar com o raio - seguindo o pensamento de que P - CB = r
    vector<float> ortogonal = criarVetor(1.0f, 1.0f, 0.0f);
    float ortogonal_norma = norma(ortogonal);
    vector<float> ortogonal_u = div_escalar(ortogonal, ortogonal_norma);

    vector<float> Circ_cone = criarPonto(CB_cone[0] + rB_cone*ortogonal_u[0], CB_cone[1] + rB_cone*ortogonal_u[1], CB_cone[2] + rB_cone*ortogonal_u[2]);
    vector<float> g_cone = subt(Circ_cone, V_cone);
    float g_cone_norma = norma(g_cone);
    vector<float> g_cone_u = div_escalar(g_cone, g_cone_norma);
    // Vamos calcular o vetor que sai do vértice e vai para o centro da base (o eixo) - pois o vetor d_cone já está normalizado
    vector<float> eixo_cone = subt(CB_cone, V_cone);
    float eixo_cone_norma = norma(eixo_cone);
    vector<float> ex_cone_u = div_escalar(eixo_cone, eixo_cone_norma);
    // Calculamos o cos(T) com o produto escalar dos vetores normalizados
    float cosT = prod_escalar(g_cone_u, ex_cone_u);


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
                /* ============= VERIFICAÇÃO PARA CONE ============= */
                // De acordo com as equações do cone também vamos calcular o vetor v_cone = V - Po
                vector<float> v_cone = subt(V_cone, olhoPintor);

                float a_cone = prod_escalar(dr_u, d_cone)*prod_escalar(dr_u, d_cone) - prod_escalar(dr_u, dr_u)*cosT*cosT;
                float b_cone = prod_escalar(v_cone, dr_u)*cosT*cosT - prod_escalar(v_cone, d_cone)*prod_escalar(dr_u, d_cone);
                float c_cone = prod_escalar(v_cone, d_cone)*prod_escalar(v_cone, d_cone) - prod_escalar(v_cone, v_cone)*cosT*cosT;

                float tEscolhido = 100.0f;
                bool intersecao = false;

                if (a_cone == 0 && b_cone!= 0) {
                    float t_cone = -c_cone/(2.0f * b_cone);

                    if (t_cone > 0.0001f) {
                        vector<float> tdr = mult_escalar(dr_u, t_cone);
                        vector<float> Pi = soma(olhoPintor, tdr);

                        vector<float> VPi = subt(V_cone, Pi);
                        float eqCone = prod_escalar(VPi, d_cone);

                        if (eqCone >= 0 && eqCone <= h_cone) {
                            intersecao = true;
                            tEscolhido = minimo(t_cone, tEscolhido);
                        }
                    }
                }

                float delta_cone = b_cone * b_cone - a_cone * c_cone;

                if (delta_cone >= 0) {
                    float t1 = (-b_cone + (sqrt(delta_cone)))/(a_cone);
                    float t2 = (-b_cone - (sqrt(delta_cone)))/(a_cone);

                    for (float ti : {t1, t2}) {
                        if (ti > 0.0001f) {
                            vector<float> tdr = mult_escalar(dr_u, ti);
                            vector<float> Pi = soma(olhoPintor, tdr);
                            vector<float> VPi = subt(V_cone, Pi);
                            float eqCone = prod_escalar(VPi, d_cone);

                            if (eqCone >= 0 && eqCone <= h_cone) {
                                intersecao = true;
                                tEscolhido = minimo(ti, tEscolhido);
                            }
                        }
                    }
                }

                if ( prod_escalar(dr_u, d_cone) != 0 ) {
                    vector<float> PoCb = subt(olhoPintor, CB_cone);
                    float t_base = -(prod_escalar(PoCb,d_cone)/ prod_escalar(dr_u,d_cone));

                    if (t_base > 0.0001f) {
                        vector<float> tdr = mult_escalar(dr_u, t_base);
                        vector<float> Pi = soma(olhoPintor, tdr);

                        vector<float> PiCB = subt(Pi, CB_cone);

                        if (prod_escalar(PiCB, PiCB) <= rB_cone*rB_cone) {
                            intersecao = true;
                            tEscolhido = minimo(t_base, tEscolhido);
                        }
                    }
                }

                if (intersecao) {
                    vector<float> tdr = mult_escalar(dr_u, tEscolhido);
                    vector<float> Pi = soma(olhoPintor, tdr);

                    // Calcula l
                    vector<float> pfPI = subt(fonte_p_coord, Pi);
                    float pfPI_norma = norma(pfPI);
                    vector<float> l_cone = div_escalar(pfPI, pfPI_norma);

                    vector<float> I_E_cone = calcular_intensidade_olho(dr_u, l_cone,d_cone,fonte_p_int,Kd_cone, Ke_cone, m_cone, I_A);

                    // Imprimir as cores no arquivo
                    if (I_E_cone[0] > 1.0f){
                        fp.put(255);
                    } else{
                        fp.put(I_E_cone[0] * 255);
                    }
                    if (I_E_cone[1] > 1.0f){
                        fp.put(255);
                    } else{
                        fp.put(I_E_cone[1] * 255);
                    }
                    if (I_E_cone[2] > 1.0f){
                        fp.put(255);
                    } else{
                        fp.put(I_E_cone[2] * 255);
                    }
                } else {

                    /* ============= VERIFICAÇÃO PARA PLANO CHÃO ============= */
                    float ti_p_c = calcula_ti_plano(olhoPintor, P_pi_c, n_bar_c, dr_u);

                    /* ============= VERIFICAÇÃO PARA PLANO FUNDO ============= */
                    float ti_p_f = calcula_ti_plano(olhoPintor, P_pi_f, n_bar_f, dr_u);

                    if (ti_p_c <= ti_p_f && ti_p_c > 0) {
                        vector<float> l;
                        float delta_esf_c = calcula_delta_plano(&l,dr_u, ti_p_c, olhoPintor, fonte_p_coord, centroEsfera, rEsfera);
                        float delta_cone_c = calcula_delta_plano(&l, dr_u, ti_p_c, olhoPintor, fonte_p_coord, CB_cone, rB_cone);

                        if (delta_esf_c >= 0 || delta_cone_c >= 0) {
                           l = criarVetor(0.0f, 0.0f, 0.0f);
                        }

                        vector<float> I_E_c = calcular_intensidade_olho(dr_u, l,n_bar_c, fonte_p_int, Kd_c, Ke_c, m_c, I_A);

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
                        vector<float> l;

                        float delta_esf_f = calcula_delta_plano(&l,dr_u, ti_p_f, olhoPintor, fonte_p_coord, centroEsfera, rEsfera);
                        float delta_cone_f = calcula_delta_plano(&l, dr_u, ti_p_f, olhoPintor, fonte_p_coord, CB_cone, rB_cone);

                        if (delta_esf_f >= 0 || delta_cone_f >= 0) {
                            l = criarVetor(0.0f, 0.0f, 0.0f);
                        }

                        vector<float> I_E_f = calcular_intensidade_olho(dr_u, l, n_bar_f, fonte_p_int, Kd_f, Ke_f, m_f, I_A);

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
    }

    /* ============= MENSAGEM DE SUCESSO =============  */
    fp.close();
    cout << "Imagem gerada com sucesso!" << endl;
    return 0;
}
