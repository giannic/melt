#ifndef VOXEL_GRID_H
#define VOXEL_GRID_H

#include "vector.h"

class VoxelGrid {

public:
	// Constructor
	VoxelGrid();
	VoxelGrid(const char* filename) {
		loadGrid(filename);
	}

	// Function
	void loadGrid(const char* filename);
    
    // If it is not in the voxelgrid then return vector of -1 
	Vector3DF inVoxelGrid(double x, double y, double z);

	Vector3DF getCellCenter(int i, int j, int k);
	
	//Dimensions of the Grid (Resolution)
	int theDim[3];
	//The Size of a voxel along each axis.
	float voxelSize[3];
	float offset[3];
	bool ***data;
    short ***adj; //adjacency list
};

#endif
