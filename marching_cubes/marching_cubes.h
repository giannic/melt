#ifndef MARCHING_CUBES
#define MARCHING_CUBES

#include <stdio.h>
#include <math.h>
//This program requires the OpenGL and GLUT libraries
// You can obtain them for free from http://www.opengl.org
#include "../common/GL/glut.h"

struct GLvector
{
        GLfloat fX;
        GLfloat fY;
        GLfloat fZ;     
};

void vIdle();
void vDrawScene(); 
void vResize(GLsizei, GLsizei);
void vKeyboard(unsigned char cKey, int iX, int iY);
void vSpecial(int iKey, int iX, int iY);

GLvoid vPrintHelp();
GLvoid vSetTime(GLfloat fTime);
GLfloat fSample1(GLfloat fX, GLfloat fY, GLfloat fZ);
GLfloat fSample2(GLfloat fX, GLfloat fY, GLfloat fZ);
GLfloat fSample3(GLfloat fX, GLfloat fY, GLfloat fZ);

GLvoid vMarchingCubes();
GLvoid vMarchCube1(GLfloat fX, GLfloat fY, GLfloat fZ, GLfloat fScale);
GLvoid vMarchCube2(GLfloat fX, GLfloat fY, GLfloat fZ, GLfloat fScale);

#endif