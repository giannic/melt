#ifndef MYDEFS
#define MYDEFS

// State of particle
enum Status { SOLID, LIQUID};

static const float AMBIENT_T = 296;
static const float DIFF_T = 0.5;
static const float MIN_T = 273.15;
static const float MAX_T = 373.15;
static const float ICE_T = 273;
static const float K_W = 0.1;
static const float K_ICE = 0.2;
static const float THERMAL_CONDUCTIVITY_ICE = 2.18;// in watts per meter kelvin
static const float THERMAL_CONDUCTIVITY_WATER = 0.58;// in watts per meter kelvin
static const float HEAT_CAP = 15;

// Init fluid system
static const float VOLMIN_X = -20;
static const float VOLMIN_Y = -20;
static const float VOLMIN_Z = 0;

static const float VOLMAX_X = 90;//20;
static const float VOLMAX_Y = 90;//20;
static const float VOLMAX_Z = 90;//40;

// Hacking for now...Need to find Good mapping
static const float INITMIN_X = 0;//-30;
static const float INITMIN_Y = 0;//-30;
static const float INITMIN_Z = 0;

static const float INITMAX_X = 60;//30;
static const float INITMAX_Y = 60;//30;
static const float INITMAX_Z = 60; //60;

// Init fluid system with OBJ
#define ADJUST_SCALE  1.6
#define ADJUST_OFFSET_X -15//-30
#define ADJUST_OFFSET_Y -10//-10
#define ADJUST_OFFSET_Z 1//-30

#endif MYDEFS