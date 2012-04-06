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
	bool inVoxelGrid(double x, double y, double z);
	Vector3DF getCellCenter(int i, int j, int k);
	
	//Dimensions of the Grid (Resolution)
	int theDim[3];
	//The Size of a voxel along each axis.
	float voxelSize[3];
	float offset[3];

protected:
	//unsigned int data
	bool ***data;

		
};

#endif
