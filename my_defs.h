#ifndef MYDEFS
#define MYDEFS

#define AMBIENT_T 1.0
#define MIN_T 0.0
#define MAX_T 1.0

// State of particle
enum Status { SOLID, LIQUID};
static const float diff_T = 0.5;
static const float HEAT_CONDUCT = 0.001;
static const float HEAT_CAP = 15;

// Init fluid system
static const float VOLMIN_X = -20;
static const float VOLMIN_Y = -20;
static const float VOLMIN_Z = 0;

static const float VOLMAX_X = 20;
static const float VOLMAX_Y = 20;
static const float VOLMAX_Z = 40;

static const float INITMIN_X = -30;
static const float INITMIN_Y = -30;
static const float INITMIN_Z = 1;

static const float INITMAX_X = 30;
static const float INITMAX_Y = 30;
static const float INITMAX_Z = 60;

// Init fluid system with OBJ
#define ADJUST_SCALE  1.6
#define ADJUST_OFFSET_X -15//-30
#define ADJUST_OFFSET_Y 0//-10
#define ADJUST_OFFSET_Z 5//-30

#endif MYDEFS