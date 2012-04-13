/**************************************************************************
	ORIGINAL AUTHOR: 
		Emud Mokhberi (emud@ucla.edu)
	MODIFIED BY:
	
	CONTRIBUTORS:
	

-----------------------------------------------
	
 ***************************************************************
 ******General License Agreement and Lack of Warranty ***********
 ****************************************************************

 This software is distributed for noncommercial use in the hope that it will 
 be useful but WITHOUT ANY WARRANTY. The author(s) do not accept responsibility
 to anyone for the consequences of using it or for whether it serves any 
 particular purpose or works at all. No guarantee is made about the software 
 or its performance.

 You are allowed to modify the source code, add your name to the
 appropriate list above and distribute the code as long as 
 this license agreement is distributed with the code and is included at 
 the top of all header (.h) files.

 Commercial use is strictly prohibited.
***************************************************************************/


#ifndef MARCHCUBES_H
#define MARCHCUBES_H

#include "impsurface.h"

class MarchCube;
class CubeEdge;
class CubeVtx;

extern int edgeTable[256];
extern int triTable[256][16];

/*
 * bit flags for the eight vertices of a cube. Set to 1 << {0...7}
 * so a single char can contain all the information for one cube. See
 * the comments above int edgeTable[] in marchcubes.cpp for a diagram of the
 * geometry. v2 is the lower-left-back corner. v2->v3 points in the +z direction,
 * v2->v6 in the +x, and v2->v1 in the +y.
 */
const int VULF = 1 << 0;
const int VULB = 1 << 1;
const int VLLB = 1 << 2;
const int VLLF = 1 << 3;
const int VURF = 1 << 4;
const int VURB = 1 << 5;
const int VLRB = 1 << 6;
const int VLRF = 1 << 7;

/*
 * bit flags for the 12 edges of a cube. Set to 1 << {0...11} so that all
 * possible configurations are stored in the range {0...4095}
 */
const int TOPLEFT		= 1 << 0;
const int LEFTBACK		= 1 << 1;
const int BOTLEFT		= 1 << 2;
const int LEFTFRONT		= 1 << 3;
const int TOPRIGHT		= 1 << 4;
const int RIGHTBACK		= 1 << 5;
const int BOTRIGHT		= 1 << 6;
const int RIGHTFRONT	= 1 << 7;
const int TOPFRONT		= 1 << 8;
const int TOPBACK		= 1 << 9;
const int BOTBACK		= 1 << 10;
const int BOTFRONT		= 1 << 11;

class MarchCube
{
public:
	MarchCube	();
	~MarchCube	();

	void	march		(IsoSurface& surface);

	void	setThreshold	(Double threshold_);
	void	setSize			(Double x, Double y, Double z);
	void	setRes			(int x, int y, int z);
	void	setCenter		(const Point3d& center_);
	void	setCenter		(Double x, Double y, Double z);

private:
	void	clearGrids		();
	void	initGrids		();
	void	initVertices	(int level, 
							 CubeVtx** toSet,
							 ImpSurface* function);
	void	setCubeFlags	();

	int		getVertex		(int cubeX, int cubeY, int edgeNum);

	/*
	 * vtxGrid is a 2x(width + 1)x(height + 1) array. It corresponds to all
	 * the vertices from all the cubes in a one-cube-deep slice of the volume
	 * being marched. vtxGrid[0] contains those in the plane behind those in
	 * vtxGrid[1] (where "forward" is measured in the +z direction). Looking in
	 * the -z direction, vtxGrid[1][0][0] corresponds to the front-lower-left corner
	 * of the lower-left cube in the slice. vtxGrid[0][0][0] the rear-lower-left from
	 * the same cube, and vtxGrid[1][width][0] the front-lower-right of the lower-right
	 * cube.
	 */
	CubeVtx***	vtxGrid;

	/*
	 * The first three dimensions of edgeGrid correspond exactly with those from vtxGrid.
	 * edgeGrid[0][0][0] points to an int* corresponding to the rear-lower-left-corner
	 * of the lower-left cube and so forth. edgeGrid[0][0][0][0] points in the +y direction
	 * from this vertex, edgeGrid[0][0][0][1] points in the +x direction, and 
	 * edgeGrid[0][0][0][2] in the +z. The values in edgeGrid are indices into the vertex
	 * array which correspond to the vertices on the aforementioned edges on the isosurface.
	 */
	int****		edgeGrid;

	/*
	 * edgeFlags is an array of values from edgeTable. It's a width x height array where 
	 * each value corresponds to a set of twelve bit-flags, each flag corresponding to
	 * one edge of the cube at index [i][j]. If a flag is '1' the edge under consideration 
	 * crosses the isosurface, otherwise the flag will be zero.
	 */
	int**		edgeFlags;

	/*
	 * Each int in vtxFlags contains a set of 8 one-bit flags, each corresponding to a vertex
	 * 
	 */
	int**		vtxFlags;

	Double		threshold;
	Double		sizex;
	Double		sizey;
	Double		sizez;
	int			resx;
	int			resy;
	int			resz;
	Point3d		center;
};

class CubeVtx
{
public:
	CubeVtx		();
	CubeVtx		(Point3d& pos_);
	CubeVtx		(Double x, Double y, Double z);
	CubeVtx		(Point3d& pos_, Double value_);
	CubeVtx		(Double x, Double y, Double z, Double value_);
	~CubeVtx	();

	Point3d				findSurface	(CubeVtx& endPoint,
									 Double threshold);

	const Point3d&	getPos		() 
							{ return pos; }
	void			setPos		(Double x, Double y, Double z);

	void			setVal		(Double value_)
							{ value = value_; }
	bool			isInside	(Double threshold) 
							{ return (value <= threshold); }
private:
	Point3d		pos;
	Double		value;
};

#endif // MARCHCUBES_H