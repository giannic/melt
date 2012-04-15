#ifndef RENDER_GRID_H
#define RENDER_GRID_H

#include "../my_defs.h"
#include "vector.h"

class RenderGrid {

public:
	// Constructor
	RenderGrid();
	~RenderGrid();

	int theDim[3];
	float voxelSize[3];
	bool ***data; // 1 and 0 for now
};

#endif
