#ifndef LUZ_HPP
#define LUZ_HPP

#include <iostream>
#include <vector>
using namespace std;

class Luz{
    private:
        vector<float> intensidade;
        vector<float> coordenadas;
        vector<float> IA;
    public:    
        Luz(vector<float> intensidade,
                vector<float> coordenadas,
                vector<float> IAfont
            ){
            this->intensidade = intensidade;
            this->coordenadas = coordenadas;
            this->IA = IA;
        }
    
        vector<float> getIntensidade(){
            return intensidade;
        } 
        
        vector<float> getCoordenadas(){
            return coordenadas;
        } 

        vector<float> getIA(){
            return IA;
        } 
};

#endif