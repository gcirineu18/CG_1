#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "cilindro.hpp"
#include <algorithm>

using namespace std;


struct EsferaStruct{
    float delta;
    float aDelta;
    float bDelta;
};

struct CilindroStruct{
    float delta;
    float aDelta;
    float bDelta;
};

struct ConeStruct{
    float delta;
    float aDelta;
    float bDelta;
    float cDelta;
};

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

float calcula_delta_plano(vector<float>*l, vector<float> dr_u, float ti_p, vector<float> olhoPintor, vector<float> fonte_p_coord, vector<float> centroEsfera, float rEsfera, Cilindro *cilindro){
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

    float cDelta = prod_escalar(s_esf, s_esf) - pow(rEsfera, 2);
    

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

Cilindro* cria_cilindro( vector<float> centroEsfera, float rEsfera){
    vector<float> centro_base = centroEsfera;
    float raio_base = rEsfera/3;
    float altura_cilindro = rEsfera * 3;
    vector<float> d_cil = criarVetor(-1/sqrt(3), 1/sqrt(3), -1/sqrt(3));
    vector<float> Kd = criarVetor( 0.2, 0.3, 0.8);
    vector<float> Ke = criarVetor( 0.2, 0.3, 0.8);
    vector<float> Ka = criarVetor( 0.2, 0.3, 0.8);

    return new Cilindro(raio_base, centro_base, altura_cilindro, d_cil, Kd, Ke, Ka);
}

CilindroStruct gera_cilindro(Cilindro* cilindro, vector<float> dr_u, vector<float> olhoPintor){
    vector<float> dCil = cilindro->getDcil();
    vector<float> baseCil = cilindro->getCentroCil();
    float raioCil = cilindro->getRaioCil();
    float delta, a ,b;

    float escalar_dr_dCil = prod_escalar(dCil, dr_u);
    vector<float> mult_dr_dCil = mult_escalar(dCil, escalar_dr_dCil);

    vector<float> w = subt(dr_u, mult_dr_dCil);


    vector<float> pb = subt(olhoPintor,baseCil);
    float escalar_pbu = prod_escalar(pb, dCil);
    vector<float> mult_dCil_pbu = mult_escalar(dCil, escalar_pbu);

    vector<float> v = subt(pb, mult_dCil_pbu);

    a = prod_escalar(w, w);
    b = prod_escalar(v, w);
    float c = prod_escalar(v, v) - pow(raioCil,2);

    delta = pow(b,2) - (a * c);
     
    return {delta, a, b};
}

float verifica_sombra_cilindro(Cilindro* cilindro, vector<float> l, vector<float> Pi, vector<float> fonte_p_coord, float t){
    vector<float> dCil = cilindro->getDcil();
    vector<float> baseCil = cilindro->getCentroCil();
    float raioCil = cilindro->getRaioCil();
    float alturaCil = cilindro->getAlturaCil();

    float a, b;

    // Calcula P(t) de acordo com a equação do Ray para o cilindro
    
    vector<float> pfPI = subt(fonte_p_coord, Pi);
    float pfPI_norma = norma(pfPI);
    l = div_escalar(pfPI, pfPI_norma);

    vector<float> lt = mult_escalar(l,t);
    vector<float> P = soma(Pi, lt);
    vector<float> PiC = subt(P, baseCil);
    float piC_dCil = prod_escalar(PiC, dCil);

    if(piC_dCil < 0 || piC_dCil > alturaCil){
        return -1.0;
    }

    float escalar_dr_dCil = prod_escalar(dCil, l);
    vector<float> mult_dr_dCil = mult_escalar(dCil, escalar_dr_dCil);

    vector<float> w = subt(l, mult_dr_dCil);


    vector<float> pb = subt(Pi,baseCil);
    float escalar_pbu = prod_escalar(pb, dCil);
    vector<float> mult_dCil_pbu = mult_escalar(dCil, escalar_pbu);

    vector<float> v = subt(pb, mult_dCil_pbu);

    a = prod_escalar(w, w);
    b = prod_escalar(v, w);
    float c = prod_escalar(v, v) - pow(raioCil,2);

    return pow(b,2) - (a * c);
}

EsferaStruct calcula_delta_esf(float x, float y, float z, vector<float> dr_u, vector<float> olhoPintor, vector<float> centroEsfera, float rEsfera){

    /* ============= VERIFICAÇÃO PARA ESFERA ============= */
    // De acordo com o processo matemático calculado, criamos o vetor w, que corresponde ao ponto olhoPintor - centroEsfera
    // P(t) = Po + t*dr -> Pj = Po + t*dr -> Pj - C = Po + t*dr - C => Po - C + t*dr => w + t*dr

    vector<float> w_esf = subt(olhoPintor, centroEsfera);

    // Com o w podemos calcular os componentes a, b e c da nossa equação do segundo grau
    float aDelta = prod_escalar(dr_u, dr_u);

    float bDelta = 2 * prod_escalar(w_esf, dr_u);

    float cDelta = prod_escalar(w_esf, w_esf) - pow(rEsfera, 2) ;

    return {(pow(bDelta,2) - (4 * aDelta * cDelta)), aDelta, bDelta};
}


ConeStruct calcula_cone(vector<float> V_cone, vector<float> dr_u, vector<float> olhoPintor, vector<float> d_cone, float cosT, float h_cone){
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

    return {delta_cone, a_cone, b_cone, c_cone};

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

    // Cone
 
    float h_con = 3.0f*rEsfera;
    vector<float> d_con = criarVetor(-1/sqrt(3), 1/sqrt(3), -1/sqrt(3));
    float d_con_norma = norma(d_con);

    vector<float> hd_con = mult_escalar(d_con, h_con);
    vector<float> CT_cone = soma(centroEsfera, hd_con);

    vector<float> CB_cone = CT_cone;
    float rB_cone = 1.5f * rEsfera;
    float h_cone = rB_cone/3.0f;
    vector<float> d_cone = div_escalar(d_con, d_con_norma);
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
    Cilindro* cilindro = cria_cilindro(centroEsfera, rEsfera);
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

            
            auto EsferaStruct = calcula_delta_esf(x, y, z, dr_u, olhoPintor, centroEsfera, rEsfera);

            auto Cstruct = gera_cilindro(cilindro, dr_u, olhoPintor);

            auto ConStruct = calcula_cone(V_cone, dr_u, olhoPintor, d_cone, cosT, h_cone);
            
            float tEsf = INFINITY;
            float tCil = INFINITY;
            float tCon = INFINITY;
            float t1, t2, t1Cil, t2Cil, t;

        /* ============= VERIFICAÇÃO PARA PLANO CHÃO ============= */
            float ti_p_c = calcula_ti_plano(olhoPintor, P_pi_c, n_bar_c, dr_u);

            /* ============= VERIFICAÇÃO PARA PLANO FUNDO ============= */
            float ti_p_f = calcula_ti_plano(olhoPintor, P_pi_f, n_bar_f, dr_u);


            // calcular menor valor de t
            if (EsferaStruct.delta >= 0) {
                t1 = (- EsferaStruct.bDelta + sqrt(EsferaStruct.delta))/(2.0f * EsferaStruct.aDelta);
                t2 = (- EsferaStruct.bDelta - sqrt(EsferaStruct.delta))/(2.0f * EsferaStruct.aDelta);
                if(t1 >= 0.0f && t2 >= 0.0f) tEsf = min(t1,t2);
                else if(t1 > 0.0f) tEsf = t1;
                else if(t2 > 0.0f) tEsf = t2;

            }

            if(Cstruct.delta >= 0){
                t1Cil = (- Cstruct.bDelta + sqrt(Cstruct.delta))/Cstruct.aDelta;
                t2Cil = (- Cstruct.bDelta - sqrt(Cstruct.delta))/Cstruct.aDelta;
                if(t1Cil >= 0.0f && t2Cil >= 0.0f) tCil = min(t1Cil,t2Cil);
                else if(t1Cil > 0.0f) tCil = t1Cil;
                else if(t2Cil > 0.0f) tCil = t2Cil;

                vector<float> dCil = cilindro->getDcil();
                vector<float> centroCil = cilindro->getCentroCil();
                float raioCil = cilindro->getRaioCil();
                vector<float> Kcil = cilindro->getKaCil();
                float alturaCil = cilindro->getAlturaCil();

                // Calcula P(t) de acordo com a equação do Ray
                vector<float> tdr = mult_escalar(dr_u, tCil);
                vector<float> Pi = soma(olhoPintor, tdr);
                vector<float> PiC = subt(Pi, centroCil);
                float piC_dCil = prod_escalar(PiC, dCil);

                if(piC_dCil< 0 || piC_dCil> alturaCil){
                    tCil = INFINITY;
                }
                
            }

            float tEscolhido = 100.0f;
                bool intersecao = false;

            if (ConStruct.aDelta == 0 && ConStruct.bDelta != 0) {
                float t_cone = -ConStruct.cDelta /(2.0f * ConStruct.bDelta);

                if (t_cone > 0.0001f) {
                    vector<float> tdr = mult_escalar(dr_u, t_cone);
                    vector<float> Pi = soma(olhoPintor, tdr);

                    vector<float> VPi = subt(V_cone, Pi);
                    float eqCone = prod_escalar(VPi, d_cone);

                    if (eqCone >= 0 && eqCone <= h_cone) {
                        intersecao = true;
                        tCon = minimo(t_cone, tEscolhido);
                    }
                }
            } else if(ConStruct.delta >= 0 ){
                float t1 = (-ConStruct.bDelta + (sqrt(ConStruct.delta)))/(ConStruct.aDelta);
                float t2 = (-ConStruct.bDelta - (sqrt(ConStruct.delta)))/(ConStruct.aDelta);

                for (float ti : {t1, t2}) {
                    if (ti > 0.0001f) {
                        vector<float> tdr = mult_escalar(dr_u, ti);
                        vector<float> Pi = soma(olhoPintor, tdr);
                        vector<float> VPi = subt(V_cone, Pi);
                        float eqCone = prod_escalar(VPi, d_cone);

                        if (eqCone >= 0 && eqCone <= h_cone) {
                            intersecao = true;
                            tCon = minimo(ti, tCon);
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
                            tCon = minimo(t_base, tCon);
                        }
                    }
                }
            }

            t = tEsf;
            if(tCil > 0.0f && (tCil < t || t < 0.0f)) t = tCil;
            if(tCon > 0.0f && (tCon < t || t < 0.0f)) t = tCon;
            if(ti_p_c > 0.0f && (ti_p_c < t || t < 0.0f)) t = ti_p_c;
            if(ti_p_f > 0.0f && (ti_p_f < t || t < 0.0f)) t = ti_p_f;   


           // Se o delta é positivo ou 0, encontramos a esfera, se não vamos testar para outros objetos utilizados na cena
            if (t == tEsf) {
               
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


            }  else if(t == tCil){

  
                vector<float> dCil = cilindro->getDcil();
                vector<float> centroCil = cilindro->getCentroCil();
                vector<float> Kcil = cilindro->getKaCil();

                /* Calcular vetores n, l, v e r */

                // Calcula P(t) de acordo com a equação do Ray
                vector<float> tdr = mult_escalar(dr_u, t);
                vector<float> Pi = soma(olhoPintor, tdr);
        
                vector<float> PiC = subt(Pi, centroCil);

                float proj = prod_escalar(PiC, dCil);
                vector<float> proj_d = mult_escalar(dCil, proj);

                // Componente perpendicular ao eixo
                vector<float> n_temp = subt(PiC, proj_d);

                // Normal do cilindro
                vector<float> n = div_escalar(n_temp, norma(n_temp));

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
                vector<float> Ifk = arroba(fonte_p_int, Kcil);
                vector<float> Id = mult_escalar(Ifk, ln_limitado);

                float rx = max(0.0f, prod_escalar(r, x));
                float rxm = pow(rx, m_esf);
                vector<float> Ie = mult_escalar(Ifk, rxm);

                vector<float> Iamb = arroba(I_A, Kcil);

                vector<float> I_E = {Ie[0] + Id[0] + Iamb[0], Ie[1] + Id[1] + Iamb[1], Ie[2] + Id[2] + Iamb[2]};

                // Imprimir as cores no arquivo
                if (I_E[0] > 1.0f){
                    fp.put(255);
                } else{
                    fp.put(I_E[0] * 255);
                }
                fp.put(I_E[1] * 255);
                fp.put(I_E[2] * 255);
                               
           } else if (t == tCon){
        
                vector<float> tdr = mult_escalar(dr_u, tCon);
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
        
            }  else {
                
                if (ti_p_c <= ti_p_f && ti_p_c > 0) {
                    vector<float> l;

                    float delta = calcula_delta_plano(&l,dr_u, ti_p_c, olhoPintor, fonte_p_coord, centroEsfera, rEsfera, cilindro);

                     // Calcula P(t) de acordo com a equação do Ray para a Esfera
                    vector<float> tdr = mult_escalar(dr_u, ti_p_c);
                    vector<float> Pi = soma(olhoPintor, tdr);

                    float deltaCil= verifica_sombra_cilindro(cilindro, l, Pi, fonte_p_coord, ti_p_c);

                    if (delta >= 0 || deltaCil >= 0) {
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

                    float delta = calcula_delta_plano(&l,dr_u, ti_p_f, olhoPintor, fonte_p_coord, centroEsfera, rEsfera, cilindro);

                    
                    // Calcula P(t) de acordo com a equação do Ray para a Esfera
                    vector<float> tdr = mult_escalar(dr_u, ti_p_f);
                    vector<float> Pi = soma(olhoPintor, tdr);

                    // Calcula P(t) de acordo com a equação do Ray para o cilindro
                    vector<float> pfPI = subt(fonte_p_coord, Pi);
                    float pfPI_norma = norma(pfPI);
                    l = div_escalar(pfPI, pfPI_norma);

                    float deltaCil= verifica_sombra_cilindro(cilindro, l, Pi, fonte_p_coord, ti_p_f);

                    if (delta >= 0 || deltaCil >= 0) {
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

    /* ============= MENSAGEM DE SUCESSO =============  */
    fp.close();
    cout << "Imagem gerada com sucesso!" << endl;
    return 0;
}
