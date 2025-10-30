#include <iostream>
#include <vector>
#include <cmath>
#include "fonteLuz.hpp"
#include "utils.hpp"

using namespace std;

class Cone {
    private:
        float raio;
        float altura;
        float m;
        float theta;
        vector<float> centro_base;
        vector<float> vertice;
        vector<float> d;
        vector<float> Kd;
        vector<float> Ke;
        vector<float> Ka;
        bool base;
    public:    
        Cone(float raio, float altura, float m,
                vector<float> centro_base, vector<float> d,
                vector<float> Kd, vector<float> Ke, vector<float> Ka,
                bool base
            ){

            vector<float> dh = mult_escalar(d, altura);

            this->raio = raio;
            this->altura = altura;
            this->vertice = soma(centro_base, dh);
            this->m = m;
            this->theta = atan(raio/altura);
            this->centro_base = centro_base;
            this->d = d;
            this->Kd = Kd;
            this->Ke = Ke;
            this->Ka = Ka;
            this->base = base;
        }
    
        vector<float> getD(){
            return d;
        } 
        vector<float> getCentroBase (){
            return centro_base;
        }  
        vector<float> getKa(){
            return Ka;
        }   
        vector<float> getKd(){
            return Kd;
        }   
        vector<float> getKe(){
            return Ke;
        }
        
        float getRaio(){
            return raio;
        } 

        float getAltura(){
            return altura;
        }

        vector<float> calculaIntersecao(vector<float> origem, vector<float> dr, Luz& fonte) {
            float cosT = cos(theta);

            vector<float> v_cone = subt(vertice, origem);
            
            float a = prod_escalar(dr, d)*prod_escalar(dr, d) - prod_escalar(dr, dr)*cosT*cosT;
            float b = prod_escalar(v_cone, dr)*cosT*cosT - prod_escalar(v_cone, d)*prod_escalar(dr, d);
            float c = prod_escalar(v_cone, d)*prod_escalar(v_cone, d) - prod_escalar(v_cone,v_cone)*cosT*cosT;

            float tEscolhido = INFINITY;
            float t = INFINITY;
            bool intersecao = false;

            if (a == 0 && b != 0) {
                t = -c/(2.0f * b);

                if (t > 0.0f) {
                    vector<float> tdr = mult_escalar(dr, t);
                    vector<float> Pi = soma(origem, tdr);

                    vector<float> VPi = subt(vertice, Pi);
                    float eqCone = prod_escalar(VPi, d);

                    if (eqCone >= 0 && eqCone <= altura) {
                        intersecao = true;
                        tEscolhido = minimo(t, tEscolhido);
                    }
                }
            }

            float delta = b * b - a * c;

            if (delta >= 0.0f) {
                float t1 = -b + (sqrt(delta)/a);
                float t2 = -b - (sqrt(delta)/a);

                for (float ti : {t1, t2}) {
                    if (ti > 0.0f) {
                        vector<float> tdr = mult_escalar(dr, ti);
                        vector<float> Pi = soma(origem, tdr);
                        vector<float> VPi = subt(vertice, Pi);
                        float eqCone = prod_escalar(VPi, d);

                        if (eqCone >= 0 && eqCone <= altura) {
                            intersecao = true;
                            tEscolhido = minimo(ti, tEscolhido);
                        }
                    }
                }
            }

            if (prod_escalar(dr, d) != 0 ) {
                vector<float> PoCb = subt(origem, centro_base);
                float t_base = -(prod_escalar(PoCb,d)/ prod_escalar(dr,d));

                if (t_base > 0.0001f) {
                    vector<float> tdr = mult_escalar(dr, t_base);
                    vector<float> Pi = soma(origem, tdr);

                    vector<float> PiCB = subt(Pi, centro_base);

                    if (prod_escalar(PiCB, PiCB) <= raio*raio) {
                        intersecao = true;
                        tEscolhido = minimo(t_base, tEscolhido);
                    }
                }
            }

            if (intersecao) {
                vector<float> IE = calculaIluminacao(origem, tEscolhido, dr, fonte);
                return IE;
            } else {
                vector<float> nulo = criarVetor(-1.0f, -1.0f, -1.0f);
                return nulo;
            }

        }

        vector<float> calculaIluminacao(vector<float> origem, float tEscolhido, vector<float> dr, Luz& fonte) {
            vector<float> fonte_coord = fonte.getCoordenadas();
            vector<float> fonte_int = fonte.getIntensidade();
            vector<float> fonte_IA = fonte.getIA();
            
            vector<float> tdr = mult_escalar(dr, tEscolhido);
            vector<float> Pi = soma(origem, tdr);
            vector<float> PiCb = subt(Pi, centro_base);
            float proj = prod_escalar(PiCb, d);

            vector<float> normal;

            if (proj < 0.0f) {
                normal = criarVetor(-d[0], -d[1], -d[2]);
            } else {
                vector<float> PiV = subt(Pi, vertice);
                float distVertice = norma(PiV);

                if (distVertice < 0.0f) {
                    normal = d;
                } else {
                    float PiVd = prod_escalar(PiV, d);
                    float cosT = cos(theta);
                    
                    vector<float> termo = mult_escalar(d, PiVd/(cosT*cosT));
                    vector<float> aux_n = subt(PiV, termo);
                    normal = div_escalar(aux_n, norma(aux_n));

                }
            }

            if (prod_escalar(normal, dr) > 0) {
                normal = criarVetor(-normal[0], -normal[1], -normal[2]);
            }
            

            // Calcula l
            vector<float> pfPI = subt(fonte_coord, Pi);
            float pfPI_norma = norma(pfPI);
            vector<float> l = div_escalar(pfPI, pfPI_norma);

            // Calcula x
            vector<float> x = criarVetor(-dr[0], -dr[1], -dr[2]);

            // Calcula r
            float ln = prod_escalar(l, normal);
            vector<float> r = {(2*ln*normal[0]) - l[0], (2*ln*normal[1]) - l[1], (2 * ln*normal[2]) - l[2]};

            // Calcular Id e Ie e Iamb
            float ln_limitado = max(0.0f, ln);
            vector<float> Ifkd = arroba(fonte_int, Kd);
            vector<float> Id = mult_escalar(Ifkd, ln_limitado);

            float rx = max(0.0f, prod_escalar(r, x));
            float rxm = pow(rx, m);
            vector<float> Ifke = arroba(fonte_int, Ke);
            vector<float> Ie = mult_escalar(Ifke, rxm);

            vector<float> Iamb = arroba(fonte_IA, Kd);

            return {Ie[0] + Id[0] + Iamb[0], Ie[1] + Id[1] + Iamb[1], Ie[2] + Id[2] + Iamb[2]};

        }
};