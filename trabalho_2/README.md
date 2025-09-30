## Tarefa 02 Esfera iluminada por uma fonte de luz pontual

Modificar o método da Tarefa 01 para que, caso haja interseção de um raio com a esfera, a cor retornada seja dada pela energia luminosa que vem do ponto de interseção PI em direção ao olho do observador.  Essa energia luminosa é o resultado da interação entre a energia luminosa emitida pela fonte pontual e o material da esfera no ponto de interseção. 
Ela é composta de duas parcelas: a reflexão DIFUSA (I_d) e a reflexão Especular. (I_e), onde

I_d =( I_F@K)* (l . n)

I_e = (I_F@K)*(v . r)^m.

Use os seguintes atributos da fonte de luz pontual:

I_F = (0.7, 0.7, 0.7)  // Intensidade da fonte pontual

P_F = (0, 5, 0)   // Posição da fonte pontual situada a 5 metros acima do olho do observador.