

#include <iostream>
#include <vector>
using namespace std;

class Cilindro{
    private:
        float raio_cilindro;
        vector<float> centro_cilindro;
        float altura_cilindro;
        vector<float> d_cil;
        vector<float> Kd;
        vector<float> Ke;
        vector<float> Ka;
    public:    
        Cilindro(float raio_cilindro, vector<float> centro_cilindro, 
                float altura_cilindro, vector<float> d_cil, vector<float> Kd, 
                vector<float> Ke, vector<float> Ka
            ){
            this->raio_cilindro = raio_cilindro;
            this->centro_cilindro = centro_cilindro;
            this->altura_cilindro = altura_cilindro;
            this->d_cil = d_cil;
            this->Kd = Kd;
            this->Ke = Ke;
            this->Ka = Ka;
        }
    
        vector<float> getDcil(){
            return d_cil;
        } 
        vector<float> getCentroCil(){
            return centro_cilindro;
        }  
        vector<float> getKaCil(){
            return Ka;
        }   
        vector<float> getKdCil(){
            return Kd;
        }   
        vector<float> getKeCil(){
            return Ke;
        }
        
        float getRaioCil(){
            return raio_cilindro;
        } 

         float getAlturaCil(){
            return altura_cilindro;
        } 
};