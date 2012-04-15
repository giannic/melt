#ifndef MYDEFS
#define MYDEFS
#include <math.h>

// State of particle
enum Status { SOLID, LIQUID};

static const float AMBIENT_T = 283;
static const float DIFF_T = 0.5;
static const float MIN_T = 253;
static const float MAX_T = 373;
static const float ICE_T = 273;
static const float K_WATER = 0.5;//71.97;
static const float K_ICE = 2.0;//75.64;
static const float MASS_H2O =  0.0008; //2.99;// * pow(10.0, -23);
static const float VISC_WATER = 0.3;
static const float EFFECTIVE_RADIUS = 0.0055; //0.01
static const float INT_STIFF = 1.0;  // 0.5
static const float EXT_STIFF = 5000;

static const float THERMAL_CONDUCTIVITY_ICE = 2.18;// in watts per meter kelvin
static const float THERMAL_CONDUCTIVITY_WATER = 0.58;// in watts per meter kelvin
static const float THERMAL_CONDUCTIVITY = 0.00267; //IUDN
static const float HEAT_CAPACITY_ICE = 2.11; // units: kJ/kg-K
static const float HEAT_CAPACITY_WATER = 4.181; // units: kJ/kg-K

// Init fluid system
static const float VOLMIN_X = -5;
static const float VOLMIN_Y = -5;
static const float VOLMIN_Z = 0;

static const float VOLMAX_X = 33;//20;
static const float VOLMAX_Y = 20;//20;
static const float VOLMAX_Z = 30;//40;

// Hacking for now...Need to find Good mapping
static const float INITMIN_X = 0;//-30;
static const float INITMIN_Y = 0;//-30;
static const float INITMIN_Z = 1;

static const float INITMAX_X = 30;
static const float INITMAX_Y = 20;
static const float INITMAX_Z = 20;

// rendering radius for voxel space
static const float RENDER_GRID_DIV = 1;
static const float RENDER_RADIUS = 1;

// Init fluid system with OBJ
#define ADJUST_SCALE 1.521587
#define ADJUST_OFFSET_X 0//-30
#define ADJUST_OFFSET_Y 0//-10
#define ADJUST_OFFSET_Z 1//-30

#endif MYDEFS