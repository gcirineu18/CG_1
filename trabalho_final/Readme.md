## Pré-Requisitos (Windows)
1. MinGW-w64: Certifique-se de ter o GCC instalado (via MSYS ou similar)
2. SDL2 & SDL2_image: *As bibliotecas devem estar instaladas em C:\SDL2.
- Certifique-se de que a estrutura seja C:\SDLS2\bin e C:\SDLS2\lib

## Estrutura
A pasta images contém as texturas utilizadas na cena.
A pasta include contém os arquivos de cabeçalho .hpp
- Mantemos uma cópia dos cabeçalhos do SDL no repositório na pastas include\SDL

## Como rodar

1. Antes de rodar, copie as DLLs da sua pasta C:\SLD2\bin para a raiz do projeto.

2. Entre na pasta trabalho_final
```bash
cd trabalho_final
```

3. Digite o comando para compilar e executar:
```bash
g++ main.cpp -o programa.exe -I.\include -L"C:\SDL2\lib" - lmingw32 - lSDL2main -lSDL2 -lSDL2_image && programa.exe
```

ou

```bash
g++ -o programa main.cpp  -lSDL2 -lSDL2_image && ./programa
```


