#include "impsurface.h"
//#include "MACGrid.h"
//#include "camera.h"
#include "GL/glut.h"

/*******************************************************************
 * Class ImpSurface
 * represents an implicit function on R^3
 ******************************************************************/
ImpSurface::ImpSurface ()
{
}

ImpSurface::~ImpSurface ()
{
}

Vector3d ImpSurface::grad (const Point3d& location)
{
	Double dx = this->eval(location[0] + GRAD_EPSILON,
						   location[1],
						   location[2]) -
				this->eval(location[0] - GRAD_EPSILON,
						   location[1],
						   location[2]);
	Double dy = this->eval(location[0],
						   location[1] + GRAD_EPSILON,
						   location[2]) -
				this->eval(location[0],
						   location[1] - GRAD_EPSILON,
						   location[2]);
	Double dz = this->eval(location[0],
						   location[1],
						   location[2] + GRAD_EPSILON) -
				this->eval(location[0],
						   location[1],
						   location[2] - GRAD_EPSILON);

	Vector3d g(dx, dy, dz);
	g *= GRAD_EPS_INV * 0.5;
	g.normalize();
	return g;
}

Vector3d ImpSurface::normal (const Point3d& location)
{
	return grad(location).normalize();
}

Double ImpSurface::eval (const Point3d& location, Double t)
{
	return eval(location);
}

/*******************************************************************
 * Class IsoSurface
 * holds a tesselation of a specified implicit function representing
 * an isosurface of that function.
 ******************************************************************/

IsoSurface::IsoSurface ()
:function(NULL)
{
}

IsoSurface::IsoSurface (ImpSurface* function_)
:function(function_)
{
}

IsoSurface::~IsoSurface ()
{
}

int IsoSurface::addVertex (Point3d& toAdd)
{
	int index = (int) vertices.size();
	try 
	{ 
		vertices.push_back(toAdd); 
	}
	catch(...)
	{
		cerr << "couldn't push vertex #" << index << endl;
	}
	return index;
}

void IsoSurface::addFace (int v1, int v2, int v3)
{
	try
	{
		faces.push_back(MeshTriangle(v1, v2, v3));
	}
	catch(...)
	{
		cerr << "couldn't push face #" << (int) faces.size() << endl;
	}
}

void IsoSurface::addFace (MeshTriangle& toAdd)
{
	try
	{
		faces.push_back(toAdd);
	}
	catch(...)
	{
		cerr << "couldn't push face #" << (int) faces.size() << endl;
	}
}

void IsoSurface::glDraw ()
{
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glVertexPointer(3, GL_DOUBLE, sizeof(Vector3d), &vertices[0][0]);
	glNormalPointer(GL_DOUBLE, sizeof(Vector3d), &vNormals[0][0]);
	glDrawElements(GL_TRIANGLES, (GLsizei) faces.size() * 3, GL_UNSIGNED_INT, &faces[0][0]);
	glDisableClientState(GL_VERTEX_ARRAY); 
    glDisableClientState(GL_NORMAL_ARRAY); 
}

void IsoSurface::calcVNorms ()
{
	int numV = (int)vertices.size();
	vNormals.clear();

	for (int i = 0; i < numV; ++i)
	{
		try
		{
			vNormals.push_back(Vector3d(0.0,0.0,0.0));
		}
		catch (...)
		{
			cerr << "couldn't push normal #" << i << endl;
			break;
		}
	}
}

void IsoSurface::clear ()
{
	vertices.clear();
	faces.clear();
	vNormals.clear();
}

ImpSurface* IsoSurface::getFunction ()
{
	return function;
}

/* 
 * Prints the Mesh if Wavefront OBJ Format
 */ 
ostream& operator << (ostream& out, const IsoSurface& s)
{
	out << "#Mesh Animation OBJ Exporter" << endl;
	out << "#Vertices" << endl;

	int numVertices = (int)s.vertices.size();
	
	for (int i = 0; i < numVertices; ++i)
	{
		out << "v " << s.vertices[i][0] << " " << s.vertices[i][1] << " " << s.vertices[i][2] << endl;
	}

	//out << "#Vertex Normals" << endl;
	//for (int i = 0; i < numVertices; ++i)
	//{
	//	out << "vn " << s.vNormals[i][0] << " " << s.vNormals[i][1] << " " << s.vNormals[i][2] << endl;
	//	out << endl;
	//}
	//out << endl;

	//out << "#Faces" << endl;
	int numFaces = (int) s.faces.size();
	for (int i = 0; i < numFaces; ++i)
	{
		out << "f " << s.faces[i] << endl;
	}

	return out;
}