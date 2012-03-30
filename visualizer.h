#ifndef VISUALIZER
#define VISUALIZER

#include "gl_helper.h"

#define hx 30.0 //half size
#define hy 30.0
#define hz 30.0

class Visualizer {

public:
    Visualizer();

    static void draw_bound();
private:
    static const int SIZE = 72;
    static const float vertices[SIZE];
    static const float colors[SIZE];
};
#endif