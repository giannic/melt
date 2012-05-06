#ifndef MYDEFS
#define MYDEFS
#include <math.h>

// State of particle
enum Status { SOLID, LIQUID};

static const float AMBIENT_T = 373; //283;
static const float C_ICE = 0.5;
static const float C_WATER = 0.1;
static const float MIN_T = 253;
static const float MAX_T = 373;
static const float ICE_T = 273;
static const float K_WATER = 0.01 ;//71.97;
static const float K_ICE = 0.5; //* 10000;//75.64;rr
static const float MASS_H2O =  0.0008; //2.99;// * pow(10.0, -23);
static const float VISC_WATER = 0.2;
static const float EFFECTIVE_RADIUS = 0.0043; //0.0053 //1.1/2;  // 0.01
static const float INT_STIFF_ICE = 0.3;  // 0.5
static const float INT_STIFF_WATER = 0.5; 
static const float EXT_STIFF = 15000;
static const float P_PRADIUS = 1.1; // 1.1;

static const float THERMAL_CONDUCTIVITY_ICE = 0.00267;//2.18;// in watts per meter kelvin
static const float THERMAL_CONDUCTIVITY_WATER = 0.00267;//0.58;// in watts perr meter kelvin
static const float THERMAL_CONDUCTIVITY = 0.00267; //IUDN
static const float HEAT_CAPACITY_ICE = 2.11; // units: kJ/kg-K
static const float HEAT_CAPACITY_WATER = 4.181; // units: kJ/kg-K

// Init fluid system
static const float VOLMIN_X = -2;
static const float VOLMIN_Y = -2;
static const float VOLMIN_Z = 0;

static const float VOLMAX_X = 40;//20;
static const float VOLMAX_Y = 20;//20;
static const float VOLMAX_Z = 40;//40;

// Hacking for now...Need to find Good mapping
static const float INITMIN_X = 2;//-30;
static const float INITMIN_Y = 2;//-30;
static const float INITMIN_Z = 5;

static const float INITMAX_X = 35;
static const float INITMAX_Y = 25;
static const float INITMAX_Z = 30;

// rendering radius for voxel space
static const float RENDER_GRID_DIV = 1;
static const float RENDER_RADIUS = 1;

// Init fluid system with OBJ
#define ADJUST_SCALE 1.521587
#define ADJUST_OFFSET_X 2//-30
#define ADJUST_OFFSET_Y 2//-10
#define ADJUST_OFFSET_Z 5//-30

// Marching cube
static const double MARCH_THRESHOLD = 0.001;
static const double MARCH_RESO = 800;

#endif MYDEFS