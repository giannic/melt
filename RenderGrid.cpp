#include "RenderGrid.h"

RenderGrid::RenderGrid() {
    theDim[0] = VOLMAX_X - VOLMIN_X;
    theDim[1] = VOLMAX_Y - VOLMIN_Y;
    theDim[2] = VOLMAX_Z - VOLMIN_Z;
	data = new bool**[theDim[0]];
	for (int i = 0; i < theDim[0]; ++i) {
		data[i] = new bool*[theDim[1]];
		for (int j = 0; j < theDim[1]; ++j) {
			data[i][j] = new bool[theDim[2]];
			for (int k = 0; k < theDim[2]; ++k) {
				data[i][j][k] = false;
			}
		}
	}
}

RenderGrid::~RenderGrid() {
	for(int i=0;i<theDim[0];i++){
		for(int j=0;j<theDim[1];j++){
			delete[] data[i][j];
		}
		delete[] data[i];
	}
	delete[] data;
}