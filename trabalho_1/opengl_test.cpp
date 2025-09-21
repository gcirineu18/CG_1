#include <GL/glut.h>
#include <iostream>


void init(){
    glClearColor(100.0f, 100.0f, 100.0f, 1.0f);

    glEnable(GL_DEPTH_TEST);
}

void display(){
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_QUADS);
        
        glColor3f(255, 0, 0);  
        
        // Define vértices do triângulo usando arrays
        float vertexA[] = {-0.0f, -0.8f, 1.0f};
        float vertex2[] = {0.0f, -0.8f, 0.0f};
        float vertex3[] = {0.0f, -0.2f, 0.0f};
        float vertex4[] = {-0.0f, -0.2f, 0.0f};
        
        // Usa glVertex3fv para definir coordenadas usando ponteiros
        glVertex3fv(vertexA);
        glVertex3fv(vertex2);
        glVertex3fv(vertex3);
        glVertex3fv(vertex4);
    glEnd();
    
    glDisable(GL_BLEND);
    glutSwapBuffers();

}


void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // Configura projeção ortogonal
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}


int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(900, 700);
    glutCreateWindow("OpenGL - Demo");

    init();

    glutDisplayFunc(display);
    //glutReshapeFunc(reshape);

    glutMainLoop();
}
