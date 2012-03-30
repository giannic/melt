#ifndef VISUALIZER
#define VISUALIZER

#include "gl_helper.h"

#define hx 20.0 //half size
#define hy 20.0
#define hz 40.0

class Visualizer {

public:
    Visualizer();

    static void draw_bound();
private:
    static const int V_SIZE = 72;
    static const int C_SIZE = 72;
    static const float vertices[V_SIZE];
    static const float colors[C_SIZE];
};
#endif