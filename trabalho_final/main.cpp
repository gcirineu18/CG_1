#include <iostream>
#include <fstream>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include "vec3.hpp"
#include "ray.hpp"
#include "camera.hpp"
#include "cena.hpp"
#include "sphere.hpp"
#include "plane.hpp"
#include "cylinder.hpp"
#include "cone.hpp"
#include "light.hpp"
#include "textura.hpp"
#include "cube.hpp"
using namespace std;

#define M_PI 3.14159265358979323846

// Função para limitar os valores entre 0 e 1
double clamp(double x) {
    if (x < 0) return 0;
    if (x > 1) return 1;
    return x;
}

int main() {
    
    // SDL_Init(SDL_INIT_VIDEO);
     
    // Inicializar SDL e SDL_image
    bool initialize = SDL_InitSubSystem(SDL_INIT_VIDEO);
    if(!initialize){
        cout<<SDL_GetError()<<endl;   
    }
    IMG_Init(IMG_INIT_JPG | IMG_INIT_PNG);


    const int width = 800;
    const int height = 800;

    SDL_Window* window = SDL_CreateWindow(
        "Cenário",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        width,
        height,
        SDL_WINDOW_SHOWN
    );

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);


    Camera cam(Vec3(0, 0, 1),Vec3(0,0,1),Vec3(0,1,0), 1,-1.0,1.0,-1.0,1.0);

    std::vector<Light> luzes;
    luzes.push_back(Light(Vec3(-1,1.4,0.2),Vec3(0.7,0.7,0.7)));
    
    Texture* texturaMadeira = new Texture();
    Texture* paredeMadeira = new Texture();
    Texture* tetoTex = new Texture();
    Texture* xGift = new Texture();
    Texture* xGift2 = new Texture();
    Texture* xGift3 = new Texture();
    Texture* xGift4 = new Texture();
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
    
    //Material material = {cor, ka, kd, ke, m, useTexture, texture}
    Material chao = {Vec3(0.2,0.2,0.2), 0.7, 0.1, 0.2, 20.0, true, texturaMadeira};
    Material parede = {Vec3(0.2,0.2,0.8), 0.686, 0.933, 0.933, 1.0, true, paredeMadeira};
    Material fundo = {Vec3(0.2,0.2,0.8), 0.686, 0.933, 0.933, 1.0, true, paredeMadeira};
    Material teto = {Vec3(0.5,0.5,0.5), 0.933, 0.933, 0.933, 1.0, true, tetoTex};
    Material cil = {Vec3(0.39,0.12,0), 0.824, 0.706, 0.549, 2.0};
    Material con = {Vec3(0.22,0.73,0.37), 0.0, 1.0, 0.498, 2.0};
    Material esf = {Vec3(1.0, 0.73, 0.37), 0.854, 0.647, 0.125, 2.0};
    Material cubo_mat = {Vec3(1.0,1.,0.37), 0.0, 1.0, 0.498, 2.0, true, xGift};  
    Material cubo_mat2 = {Vec3(1.0,1.,0.37), 0.0, 1.0, 0.498, 2.0, true, xGift2};  
    Material cubo_mat3 = {Vec3(1.0,1.,0.37), 0.0, 1.0, 0.498, 2.0, true, xGift3}; 
    Material cubo_mat4 = {Vec3(1.0,1.,0.37), 0.0, 1.0, 0.498, 2.0, true, xGift4}; 
   
    
    Cube* cube1 = new Cube(Vec3(-0.1, -1.3,-1.4), 0.5, cubo_mat, 9);
    cube1->rotateY(M_PI / 4);
    cube1->translate(0.2,0,0);


    Cube* cube2 = new Cube(Vec3(-0.9, -1.3,-1.289), 0.5, cubo_mat3, 11);
    cube2->scaleTransform(1.5, 1.0, 1.0);
    cube2->rotateY(M_PI / 4);

    Cube* cube3 = new Cube(Vec3(0.6, -1.3,-1.2), 0.5, cubo_mat2, 10);
    cube3->scaleTransform(1.5, 1.0, 1.0);

    Cube* cube4 = new Cube(Vec3(-0.3, -1.3,-1.156), 0.5, cubo_mat4, 16);
    
    cube4->scaleTransform(0.5, 0.8, 2.0);
    cube4->rotateY(M_PI/6 );



    Cone* cone = new Cone(Vec3(-0.9, -1.3,-0.9), Vec3(0,1,0), 0.9, 1.5, con, 7);

    
    cone->rotateX(M_PI); 

    cone->scaleTransform(0.1, 0.1, 0.1);


   // cone->shearTransform(0.3, 0, 0, 0, 0, 0);  

    Material prateleira = {Vec3(0.55, 0.27, 0.07), 0.3, 0.7, 0.2, 10.0, false, nullptr}; // Cor marrom madeira
    Material decoracao1 = {Vec3(0.8, 0.1, 0.1), 0.2, 0.8, 0.5, 32.0};  // Vermelho brilhante
    Material decoracao2 = {Vec3(0.2, 0.5, 0.9), 0.2, 0.8, 0.5, 32.0};  // Azul
    Material decoracao3 = {Vec3(0.9, 0.8, 0.1), 0.2, 0.8, 0.6, 64.0};  // Dourado

    Cube* prateleira1 = new Cube(Vec3(1.7, -0.5,-1.7), 0.3, prateleira, 12);
    prateleira1->scaleTransform(1.23, 0.3, 4.0);  

    Cube* prateleira2 = new Cube(Vec3(1.7, -0.2,-1.7), 0.3, prateleira, 13);
    prateleira2->scaleTransform(1.23, 0.3, 4.0); 
    

    Cone* coneDec1 = new Cone(Vec3(1.65, -0.15, -1.7), Vec3(0,1,0), 0.5, 0.8, decoracao1, 14);

    Cone* coneDec2 = new Cone(Vec3(1.65, -0.15, -1.7), Vec3(0,1,0), 0.4, 0.7, decoracao1, 15);

    coneDec1->scaleTransform(0.1, 0.1, 0.1);
    coneDec1->translate(0.37, -0.036, 0.4);
    
    coneDec2->scaleTransform(0.07, 0.07, 0.07);
    coneDec2->rotateX(M_PI);
    coneDec2->translate(0.185, 0.036, 0.45);

    Cena cena;
        

    cena.adicionar(prateleira1);
    cena.adicionar(prateleira2);
    cena.adicionar(coneDec1);
    cena.adicionar(coneDec2);

    cena.adicionar(cube1);
    cena.adicionar(cube3);
    cena.adicionar(cube2);
    cena.adicionar(cube4);

    cena.adicionar(new Plane(Vec3(0,-1.5,0),Vec3(0,1,0),chao,1));
    cena.adicionar(new Plane(Vec3(2,-1.5,0),Vec3(-1,0,0),parede,2));
    cena.adicionar(new Plane(Vec3(2,-1.5,-4),Vec3(0,0,1),fundo,3));
    cena.adicionar(new Plane(Vec3(-2,-1.5,0),Vec3(1,0,0),parede,4));
    cena.adicionar(new Plane(Vec3(0,1.5,0),Vec3(0,-1,0),teto,5));

    cena.adicionar(new Cylinder(Vec3(0, -1.5,-2.0),Vec3(0,1,0),0.05, 0.9, cil,6));
    cena.adicionar(new Cone(Vec3(0, -0.6, -2), Vec3(0,1,0),0.9,1.5,con,7));
    cena.adicionar(new Sphere(Vec3(0, 0.95,-2.0),0.05,esf,8));

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
            double u = (double)x /(width - 1);
            double v = (double)y /(height - 1);

            Ray r = cam.gerarRaio(u, v);
            HitRecord hit;
            Vec3 corPixel;

            if (cena.trace(r, hit)) {
                Vec3 baseColor = hit.mat.color;
                
                // Aplicar textura se disponível
                if (hit.mat.useTexture && hit.mat.texture != nullptr) {
                    baseColor = hit.mat.texture->sample(hit.u, hit.v);
                }
                
                Vec3 corFinal = baseColor * hit.mat.ka;

                for (const auto& luz : luzes){
                    Vec3 L = (luz.position - hit.p).normalize();
                    double distL = (luz.position - hit.p).length();

                    Ray shadowRay(hit.p + hit.normal * 0.01, L);
                    HitRecord shadowHit;

                    bool emSombra =  false;
                    if (cena.trace(shadowRay, shadowHit)){
                        if (shadowHit.objectID != hit.objectID && shadowHit.t < distL) emSombra = true;
                    }

                    if(!emSombra){
                        double dotNL = std::max(0.0, hit.normal.dot(L));
                        corFinal = corFinal + (baseColor * luz.color)* (hit.mat.kd * dotNL);

                        Vec3 V_d = (cam.eye - hit.p).normalize();
                        Vec3 R_d = (hit.normal * (2.0 * hit.normal.dot(L))) - L;
                        double dotRV = std::max(0.0, R_d.dot(V_d));

                        corFinal = corFinal + (luz.color) * (hit.mat.ks * std::pow(dotRV, hit.mat.shininess));
                    }
                }
                corPixel = corFinal;

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