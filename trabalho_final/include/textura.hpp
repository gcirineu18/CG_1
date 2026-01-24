#ifndef TEXTURE_HPP
#define TEXTURE_HPP

#include <vector>
#include <string>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include "vec3.hpp"

class Texture {
public:
    std::vector<Vec3> pixels;
    int width, height;

    Texture() : width(0), height(0) {}

    ~Texture() {}

    // Carregar imagem JPEG/PNG usando SDL2_image
    bool loadImage(const std::string& filename) {
        SDL_Surface* surface = IMG_Load(filename.c_str());
        if (!surface) {
            return false;
        }

        width = surface->w;
        height = surface->h;
        
        // Converter para formato RGB se necessário
        SDL_Surface* rgbSurface = SDL_ConvertSurfaceFormat(surface, SDL_PIXELFORMAT_RGB24, 0);
        SDL_FreeSurface(surface);
        
        if (!rgbSurface) {
            return false;
        }

        pixels.resize(width * height);
        
        Uint8* pixelData = (Uint8*)rgbSurface->pixels;
        int pitch = rgbSurface->pitch;
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int offset = y * pitch + x * 3;
                pixels[y * width + x] = Vec3(
                    pixelData[offset] / 255.0,
                    pixelData[offset + 1] / 255.0,
                    pixelData[offset + 2] / 255.0
                );
            }
        }
        
        SDL_FreeSurface(rgbSurface);
        return true;
    }

    // retornar a cor de um pixel da textura com coordenadas UV
    Vec3 sample(float u, float v) const {
        if (pixels.empty()) return Vec3(1, 0, 1); // Magenta se erro
        
        // Repetir textura (wrap)
        u = u - floor(u);
        v = v - floor(v);
        
        int i = (int)(u * (width - 1));
        int j = (int)(v * (height - 1));
        
        i = std::max(0, std::min(width - 1, i));
        j = std::max(0, std::min(height - 1, j));
        
        // Acessa o vetor 1D usando fórmula 2D
        return pixels[j * width + i];
    }
};

#endif