#ifndef MYDEFS
#define MYDEFS

#define AMBIENT_T 1.0
#define MIN_T 0.0
#define MAX_T 1.0

// State of particle
enum Status { SOLID, LIQUID};
static const float diff_T = 0.5;
static const float heat_conduct = 0.001;
static const float heat_cap = 15;

#endif MYDEFS