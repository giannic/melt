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


#ifndef IMPSURFACE_H
#define IMPSURFACE_H

#include "geometry.h"

class ImpSurface;
class IsoSurface;
class Camera;

const Double GRAD_EPSILON = 1e-5;
const Double GRAD_EPS_INV = 1.0 / GRAD_EPSILON;

class ImpSurface
{
public:
	ImpSurface			();
	virtual	~ImpSurface	();

	virtual Double	eval	(Double x, Double y, Double z) 
								{ return eval(Point3d(x, y, z)); }
	virtual Double	eval	(Double x, Double y, Double z, Double t) 
								{ return eval(Point3d(x, y, z), t); }
	virtual Double	eval	(const Point3d& location) = 0;
	virtual Double	eval	(const Point3d& location, Double t);

	Vector3d		grad	(const Point3d& location);
	Vector3d		normal	(const Point3d& location);
private:
	
};

class IsoSurface
{
public:
	IsoSurface			();
	IsoSurface			(ImpSurface* function_);
	~IsoSurface			();

	int		addVertex	(Point3d& toAdd);
	void	addFace		(int v1, int v2, int v3);
	void	addFace		(MeshTriangle& toAdd);

	void	glDraw		();

	void	calcVNorms	();

	void	clear		();

	ImpSurface*	getFunction	();

	friend ostream& operator <<		(ostream& out, const IsoSurface& s);
private:
	ImpSurface*				function;
	vector<Point3d>			vertices;
	vector<MeshTriangle>	faces;
	vector<Vector3d>		vNormals;
};

#endif // IMPSURFACE_H