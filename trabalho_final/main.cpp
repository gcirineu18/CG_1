#include <iostream>
#include <fstream>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include "include/vec3.hpp"
#include "include/ray.hpp"
#include "include/camera.hpp"
#include "include/cena.hpp"
#include "include/sphere.hpp"
#include "include/plane.hpp"
#include "include/cylinder.hpp"
#include "include/cone.hpp"
#include "include/light.hpp"
#include "include/textura.hpp"
#include "include/cube.hpp"
#include "include/spotLight.hpp"
using namespace std;

#define M_PI 3.14159265358979323846

// Função para limitar os valores entre 0 e 1
float clamp(float x) {
    if (x < 0) return 0;
    if (x > 1) return 1;
    return x;
}

int main(int argc, char* argv[]) {
    
    // SDL_Init(SDL_INIT_VIDEO);
     
    // Inicializar SDL e SDL_image
    bool initialize = SDL_InitSubSystem(SDL_INIT_VIDEO);
    if(!initialize){
        cout<<SDL_GetError()<<endl;   
    }
    IMG_Init(IMG_INIT_JPG | IMG_INIT_PNG);


    const int width = 700;
    const int height = 700;

    SDL_Window* window = SDL_CreateWindow(
        "Cenário",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        width,
        height,
        SDL_WINDOW_SHOWN
    );

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    //Camera cam(Vec3(2.8, 1.6, 3.0),Vec3(3.7, 1, 2.3),Vec3(0,1,0), 1,-1.0,1.0,-1.0,1.0);
     Camera cam(
         Vec3(2.0, 1.5, 5.0),    // eye
         Vec3(2.0, 1.5, 4.0),    // at
         Vec3(0,1,0),            // up
         1,                      // d => Alterar o valor para zoom in (> 1) e zoom out (< 1)
         -1.0, 1.0, -1.0, 1.0
     );

    std::vector<Light> luzes;
    luzes.push_back(Light(Vec3(0.75, 0.35, 0.15),Vec3(0.2,0,0)));
    luzes.push_back(Light(Vec3(1.0, 2.9, 4.2),Vec3(0.2,0.2,0.2)));

    std::vector<SpotLight> spotLights;
    spotLights.push_back(SpotLight(
        Vec3(3.95, 2.75, 2.3),        
        Vec3(0.8, 0.8, 0.8),          
        Vec3(0, -1, 0),               
        M_PI / 8                  
    )); 

    Texture* texturaMadeira = new Texture();
    Texture* paredeMadeira = new Texture();
    Texture* tetoTex = new Texture();
    Texture* xGift = new Texture();
    Texture* xGift2 = new Texture();
    Texture* xGift3 = new Texture();
    Texture* xGift4 = new Texture();
    Texture* tronco = new Texture();
    Texture* fogo = new Texture();
    Texture* tijolo = new Texture();

    if (!texturaMadeira->loadImage("images/madeira.jpeg") ) {
        std::cerr << "Erro ao carregar textura madeira.jpeg!" << std::endl;
    }
    if (!paredeMadeira->loadImage("images/parede_madeira.jpeg") ) {
        std::cerr << "Erro ao carregar textura madeira.jpeg!" << std::endl;
    }
    if (!xGift->loadImage("images/embrulho_natal.jpg") ) {
        std::cerr << "Erro ao carregar embrulho_natal.jpg!" << std::endl;
    }
    if (!xGift2->loadImage("images/embrulho_natal2.jpg") ) {
        std::cerr << "Erro ao carregar embrulho_natal2.jpg!" << std::endl;
    }
    if (!xGift3->loadImage("images/embrulho_natal3.jpg") ) {
        std::cerr << "Erro ao carregar embrulho_natal3.jpg!" << std::endl;
    }
    if (!xGift4->loadImage("images/finlandia_gift.jpg") ) {
        std::cerr << "Erro ao carregar finlandia_gift.jpg!" << std::endl;
    }
    if (!tetoTex->loadImage("images/teto.jpg") ) {
        std::cerr << "Erro ao carregar teto.jpg!" << std::endl;
    }
    if (!tronco->loadImage("images/tronco.jpg") ) {
        std::cerr << "Erro ao carregar textura tronco.jpg!" << std::endl;
    }
    if (!fogo->loadImage("images/fogo.jpg") ) {
        std::cerr << "Erro ao carregar textura madeira.jpeg!" << std::endl;
    }
    if (!tijolo->loadImage("images/tijolo.jpg") ) {
        std::cerr << "Erro ao carregar textura madeira.jpeg!" << std::endl;
    }
    
    float ka = 0.2;
    
    //Material material = {cor, kd, ke, ka, m, useTexture, texture}
    Material chao = {Vec3(0.2,0.2,0.2), 0.7, 0.1, ka, 20.0, true, texturaMadeira};
    Material parede = {Vec3(0.2,0.2,0.8), 0.686, 0.933, ka, 1.0, true, paredeMadeira};
    //Material fundo = {Vec3(0.2,0.2,0.8), 0.686, 0.933, ka, 1.0, true, paredeMadeira};
    Material teto = {Vec3(0.5,0.5,0.5), 0.933, 0.933, ka, 1.0, true, tetoTex};
    Material cil = { Vec3(0.39,0.12,0), 0.824, 0.706, ka, 2.0};
    Material cadeira = { Vec3(0.22, 0.1, 0.01), 0.824, 0.706, ka, 2.0};
    Material con = {Vec3(0.22,0.73,0.37), 0.0, 1.0, ka, 2.0};
    Material esf = {Vec3(1.0, 0.73, 0.37), 0.854, 0.647, ka, 2.0};
    Material cubo_mat = {Vec3(1.0,1.,0.37), 0.0, 1.0, ka, 2.0, true, xGift};  
    Material cubo_mat2 = {Vec3(1.0,1.,0.37), 0.0, 1.0, ka, 2.0, true, xGift2};  
    Material cubo_mat3 = {Vec3(1.0,1.,0.37), 0.0, 1.0, ka, 2.0, true, xGift3}; 
    Material cubo_mat4 = {Vec3(1.0,1.,0.37), 0.0, 1.0, ka, 2.0, true, xGift4};
    Material esf2 = {Vec3(0.9, 0.95, 1.0), 0.8, 0.6, ka, 128.0};
    Material suporte_esf= {Vec3(0.1, 0.2, 0.9), 0.5, 0.8, ka, 32.0};
    Material prateleira = {Vec3(0.55, 0.27, 0.07), 0.3, 0.7, ka, 10.0, false, nullptr}; 
    Material decoracao1 = {Vec3(0.8, 0.1, 0.1), 0.5, 0.8, ka, 32.0}; 
    Material decoracao2 = {Vec3(0.2, 0.5, 0.9), 0.2, 0.8, ka, 32.0}; 
    Material filme1 = {Vec3(0.6, 0.2, 0.8), 0.2, 0.8, ka, 64.0};
    Material filme2 = {Vec3(0.0, 0.8, 0.8), 0.2, 0.8, ka, 64.0};
    Material filme3 = {Vec3(0.1, 0.8, 0.1), 0.2, 0.8, ka, 64.0};
    Material lustre = {Vec3(0.9, 0.7, 0.0), 0.3, 0.7, ka, 10.0};
    Material lenha = {Vec3(0.2,0.2,0.8), 0.686, 0.933, ka, 1.0, true, tronco};
    Material lareira = {Vec3(0.2,0.2,0.8), 0.686, 0.933, ka, 1.0, true, fogo};
    Material chamine = {Vec3(0.2,0.2,0.8), 0.686, 0.933, ka, 1.0, true, tijolo};

    
    Cube* cube1 = new Cube(Vec3(2.1, 0.2, 2.6), 0.5, cubo_mat, "Gift1");
    cube1->rotateY(M_PI / 4);
    cube1->translate(0.2,0,0);

    Cube* cube2 = new Cube(Vec3(1.68, 0.2, 2.4), 0.5, cubo_mat3, "Gift2");
    cube2->scaleTransform(1.2, 0.7, 0.7);
    cube2->rotateY(M_PI / 4);

    Cube* cube3 = new Cube(Vec3(2.8, 0.2, 2.8), 0.5, cubo_mat2, "Gift3");
    cube3->scaleTransform(1.5, 1.0, 1.0);

    Cube* cube4 = new Cube(Vec3(1.9, 0.2, 2.844), 0.5, cubo_mat4, "Gift4");  
    cube4->scaleTransform(0.5, 0.8, 2.0);
    cube4->rotateY(M_PI/6 );


    Cube* prateleira1 = new Cube(Vec3(3.85, 1.0, 2.3), 0.3, prateleira, "prateleira1");
    prateleira1->scaleTransform(1.23, 0.3, 4.0);  

    Cube* prateleira2 = new Cube(Vec3(3.85, 1.3, 2.3), 0.3, prateleira, "prateleira2");
    prateleira2->scaleTransform(1.23, 0.3, 4.0); 

    Cube* movie1 = new Cube(Vec3(3.85, 1.14, 2.35), 0.3, filme1, "movie1");
    movie1->scaleTransform(1.0, 0.4, 1.0); 
    movie1->rotateX(M_PI/2);
    movie1->shearTransform(0.3, 0, 0, 0, 0, 0); 

    Cube* movie2 = new Cube(Vec3(3.85, 1.14, 2.42), 0.3, filme2, "movie2");
    movie2->scaleTransform(1.0, 0.4, 1.0); 
    movie2->rotateX(M_PI/2);
    movie2->shearTransform(0.3, 0, 0, 0, 0, 0); 
    
    Cube* movie3 = new Cube(Vec3(3.85, 1.14, 2.28), 0.3, filme3, "movie3");
    movie3->scaleTransform(1.0, 0.4, 1.0); 
    movie3->rotateX(M_PI/2);
    movie3->shearTransform(0.3, 0, 0, 0, 0, 0); 

    Cone* coneDec1 = new Cone(Vec3(0.0, 0.0, 0.0), Vec3(0,1,0), 0.5, 0.8, decoracao1, true, "ampulheta1");
    coneDec1->scaleTransform(0.3, 0.3, 0.3);
    coneDec1->translate(3.85, 1.33, 2.42);
    
    Cone* coneDec2 = new Cone(Vec3(0.0, 0.0, 0.0), Vec3(0,1,0), 0.4, 0.7, decoracao1, true, "ampulheta2");  
    coneDec2->scaleTransform(0.3, 0.3, 0.3);
    coneDec2->rotateX(M_PI);
    coneDec2->translate(3.85, 1.68, 2.42);

    Sphere* cristal = new Sphere(Vec3(3.82, 1.46, 2.22), 0.06, esf2, "cristal");

    Cube* suporte = new Cube(Vec3(3.82, 1.369, 2.22), 0.3, suporte_esf, "suporte cristal");
    suporte->scaleTransform(1.0, 0.2, 1.0); 

    Cone* lustre1 = new Cone(Vec3(0.0, 0.0, 0.0), Vec3(0,1,0), 0.4, 0.7, lustre, false, "lustre");
    lustre1->scaleTransform(0.6, 0.6, 0.9); 
    lustre1->translate(4.0,2.5,2.3);

    Cube* lenha1 = new Cube(Vec3(0.0,0.0,0.0),0.5,lenha,"lenha1");
    lenha1->scaleTransform(2.0,0.5,0.5);
    lenha1->rotateY(-M_PI/6);
    lenha1->translate(0.75, 0.0, 0.0);

    Cube* lenha2 = new Cube(Vec3(0.0,0.0,0.0),0.5,lenha,"lenha2");
    lenha2->scaleTransform(2.0,0.5,0.5);
    lenha2->rotateY(M_PI/6);
    lenha2->translate(0.75, 0.0, 0.0);

    Cube* lareira1 = new Cube(Vec3(0.0,0.0,0.0),0.5, lareira, "lareira1");
    lareira1->scaleTransform(2.0,2.0,0.1);
    lareira1->translate(0.75, 0.35, 0.05);

    Cube* tijolo1 = new Cube(Vec3(0.0,0.0,0.0), 0.5, chamine, "tijolo1");
    tijolo1->scaleTransform(0.5, 2.0, 0.5);
    tijolo1->translate(0.4,0.35,0.0);
    
    Cube* tijolo2 = new Cube(Vec3(0.0,0.0,0.0), 0.5, chamine, "tijolo2");
    tijolo2->scaleTransform(0.5, 2.0, 0.5);
    tijolo2->translate(1.1,0.35,0.0);

    Cube* tijolo3 = new Cube(Vec3(0.0,0.0,0.0), 0.5, chamine, "tijolo3");
    tijolo3->scaleTransform(2.6, 0.5, 0.5);
    tijolo3->translate(0.75,0.766,0.0);

    Cube* tijolo4 = new Cube(Vec3(0.0,0.0,0.0), 0.5, chamine, "tijolo4");
    tijolo4->scaleTransform(2.0, 6.5, 0.5);
    tijolo4->translate(0.75,1.93,0.0);

    Cylinder* peCadeira1 = new Cylinder(Vec3(1.16, 0.0, 2.85),Vec3(0,1,0),0.02, 0.4, cadeira, true, true, "pe1");
    Cylinder* peCadeira2 = new Cylinder(Vec3(0.96, 0.0, 3.12),Vec3(0,1,0),0.02, 0.4, cadeira, true, true, "pe2");
    Cylinder* peCadeira3 = new Cylinder(Vec3(0.78, 0.0, 2.72),Vec3(0,1,0),0.02, 0.4, cadeira, true, true, "pe3");
    Cylinder* peCadeira4 = new Cylinder(Vec3(0.63, 0.0, 2.94),Vec3(0,1,0),0.02, 0.4, cadeira, true, true, "pe4");
    Cube* assentoCadeira = new Cube(Vec3(0.9, 0.4, 2.9), 0.3, cadeira, "assento");
    assentoCadeira->scaleTransform(1.8, 0.2, 2.2); 
    assentoCadeira->rotateY(M_PI/3); 

    Cube* encostoCadeira = new Cube(Vec3(0.72, 0.65, 2.82), 0.3, cadeira, "encosto");
    encostoCadeira->scaleTransform(1.8, 0.2, 2.5); 
    encostoCadeira->rotateX(M_PI/2);
    encostoCadeira->rotateY(M_PI/3);

    Cena cena;

    cena.adicionar(movie1);
    cena.adicionar(movie2);
    cena.adicionar(movie3);
        
    cena.adicionar(cristal);
    cena.adicionar(suporte);

    cena.adicionar(prateleira1);
    cena.adicionar(prateleira2);
    cena.adicionar(coneDec1);
    cena.adicionar(coneDec2);

    cena.adicionar(cube1);
    cena.adicionar(cube3);
    cena.adicionar(cube2);
    cena.adicionar(cube4);

   cena.adicionar(lustre1);

    cena.adicionar(lenha1);
    cena.adicionar(lenha2);
    cena.adicionar(lareira1);
    cena.adicionar(tijolo1);
    cena.adicionar(tijolo2);
    cena.adicionar(tijolo3);
    cena.adicionar(tijolo4);

    cena.adicionar(peCadeira1);
    cena.adicionar(peCadeira2);
    cena.adicionar(peCadeira3);
    cena.adicionar(peCadeira4);
    cena.adicionar(assentoCadeira);
    cena.adicionar(encostoCadeira);

    cena.adicionar(new Plane(Vec3(2.0, 0.0, 4.0),Vec3(0,1,0),chao, "chao"));
    cena.adicionar(new Plane(Vec3(4.0, 0.0, 4.0),Vec3(-1,0,0),parede, "parede"));
    cena.adicionar(new Plane(Vec3(0.0, 0.0, 0.0),Vec3(0,0,1),parede, "fundo"));
    cena.adicionar(new Plane(Vec3(0.0, 0.0, 4.0),Vec3(1,0,0),parede, "parede"));
    cena.adicionar(new Plane(Vec3(Vec3(2.0, 3.0, 4.0)),Vec3(0,-1,0),teto, "teto"));

    cena.adicionar(new Cylinder(Vec3(2.2, 0.0, 2.0),Vec3(0,1,0),0.05, 0.9, cil, true, true, "tronco"));
    cena.adicionar(new Cone(Vec3(2.2, 0.9, 2.0), Vec3(0,1,0),0.9,1.5,con, true, "arvore"));
    cena.adicionar(new Sphere(Vec3(2.2, 2.45, 2.0),0.05,esf, "topo arvore"));

    std::vector<uint8_t> framebuffer(width * height * 3);

    SDL_Texture* screenTexture = SDL_CreateTexture(
        renderer,
        SDL_PIXELFORMAT_RGB24,
        SDL_TEXTUREACCESS_STREAMING,
        width,
        height
    );

    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            float u = (float)x /(width - 1);
            float v = (float)y /(height - 1);

            Ray r = cam.gerarRaio(u, v);
            HitRecord hit;
            Vec3 corPixel;

            if (cena.trace(r, hit)) {               
                
                Vec3 corLuzes = Light::calcula_luzes(luzes, hit, cena, cam);
                Vec3 corSpot = SpotLight::calcula_luzSpot(spotLights, hit, cena, cam);
                corPixel = corSpot + corLuzes;

            } else {
                corPixel = cena.bg_color;
            }

            int idx = ((height - 1 - y) * width + x) * 3;

            framebuffer[idx + 0] = (uint8_t)(clamp(corPixel.x) * 255);
            framebuffer[idx + 1] = (uint8_t)(clamp(corPixel.y) * 255);
            framebuffer[idx + 2] = (uint8_t)(clamp(corPixel.z) * 255);
        }
    }


    std::cout << "Imagem gerada com sucesso!" << std::endl;

    SDL_UpdateTexture(
        screenTexture,
        nullptr,
        framebuffer.data(),
        width * 3
    );

    SDL_RenderClear(renderer);
    SDL_RenderCopy(renderer, screenTexture, nullptr, nullptr);
    SDL_RenderPresent(renderer);
    
    delete texturaMadeira;
    IMG_Quit();
    bool running = true;
    SDL_Event e;

    while (running) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT)
                running = false;
        }
        SDL_Delay(16);
    }
    SDL_Quit();
    
    return 0;
}