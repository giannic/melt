#include "visualizer.h"
const float Visualizer::vertices [V_SIZE]= {
        // front
        -hx, hy, hz,
        -hx, -hy, hz,

        -hx, -hy, hz,
        hx, -hy, hz,

        hx, -hy, hz,
        hx, hy, hz,

        hx, hy, hz,
        -hx, hy, hz,

        //back
        -hx, hy, 0,
        -hx, -hy, 0,

        -hx, -hy, 0,
        hx, -hy, 0,

        hx, -hy, 0,
        hx, hy, 0,

        hx, hy, 0,
        -hx, hy, 0,

        //sides
        -hx, hy, hz,
        -hx, hy, 0,

        -hx, -hy, hz,
        -hx, -hy, 0,

        hx, -hy, hz,
        hx, -hy, 0,

        hx, hy, hz,
        hx, hy, 0
};

const float Visualizer::colors [C_SIZE]= {
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,

    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,

    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,

    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,

    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,

    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f
};


Visualizer::Visualizer()
{
}

void Visualizer::draw_bound()
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, vertices);
    glColorPointer(3, GL_FLOAT, 0, colors);
    glDrawArrays(GL_LINES, 0, V_SIZE/3);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
}