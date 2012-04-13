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


#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "main.h"

class Tuple3d;
class Point3d;
class Vector3d;
class Triangle;

class Tuple3d
{
 public:
  Tuple3d				();
  Tuple3d				(const Tuple3d& toCopy);
  virtual ~Tuple3d		();
  Tuple3d& operator =	(const Tuple3d& toCopy);

  Double const * ptr	() const
				{ return data; }
  Double&	operator []	(int index) 
				{ return data[index]; }
  Double	operator []	(int index) const 
				{ return data[index]; }
  Double&	get			(int index)
				{ return data[index]; }
  Double	get			(int index) const
				{ return data[index]; }

  virtual void glLoad	() = 0;
  
  friend ostream& operator << (ostream& out, const Tuple3d& toPrint);
 protected:
  Double data[3];
};

class Point3d : public Tuple3d
{
 public:
  Point3d			();
  Point3d			(Double x, Double y, Double z);
  ~Point3d	();
  Point3d			(const Point3d& toCopy);
  Point3d& operator =	(const Point3d& toCopy);

  bool   operator ==	(const Point3d& toCompare);
  bool   operator !=	(const Point3d& toCompare);

  Point3d  operator +	(const Vector3d& toAdd) const;
  Point3d& operator +=	(const Vector3d& toAdd);

  Point3d  operator *	(const Point3d& toMult) const;
  Point3d& operator *=	(const Point3d& toMult);

  Point3d  operator /	(const Point3d& toDiv) const;
  Point3d& operator /=	(const Point3d& toDiv);

  Point3d  operator *	(Double scaleFactor) const;
  Point3d& operator *=	(Double scaleFactor);
  
  Point3d  operator -	(const Vector3d& toSub);
  Point3d& operator -=	(const Vector3d& toSub);

  friend Point3d operator - (const Vector3d& v, const Point3d& p);

  void inline glLoad	();
};


class Vector3d : public Tuple3d
{
 public:
  Vector3d		();
  Vector3d		(Double x, Double y, Double z);
			 
  virtual ~Vector3d	();
  Vector3d		(const Point3d& dest);
  Vector3d		(const Point3d& start, const Point3d& end);
  Vector3d		(const Vector3d& toCopy);

  Vector3d& operator =	(const Point3d& dest); 
  Vector3d& operator =	(const Vector3d& toCopy); 
  
  bool    operator ==	(const Vector3d& toCompare) const;
  bool    operator !=	(const Vector3d& toCompare) const;

  Vector3d  operator +	(const Vector3d& toAdd) const;
  Vector3d& operator +=	(const Vector3d& toAdd);
  
  Vector3d  operator -    () const;
  Vector3d  operator -	(const Vector3d& toSub) const;
  Vector3d& operator -=	(const Vector3d& toSub);

  Vector3d  operator *	(const Vector3d& toMult) const;
  Vector3d& operator *=	(const Vector3d& toMult);

  Vector3d  operator /	(const Vector3d& toDiv) const;
  Vector3d& operator /=	(const Vector3d& toDiv);

  Vector3d  operator *	(Double scaleFactor) const;
  Vector3d& operator *=	(Double scaleFactor);

  friend Vector3d operator * (Double scaleFactor, const Vector3d& v);

  Vector3d  operator /	(Double scaleFactor) const;
  Vector3d& operator /=	(Double scaleFactor);

  Double dot			(const Vector3d& toDot);
  Vector3d cross          (const Vector3d& toCross);

  Vector3d& normalize 	();
  Vector3d  getUnit 	();

  Double length			();
  void   inline glLoad 	();
};


class Color3d : public Tuple3d 
{
 public:
  Color3d			();
  Color3d			(Double r, Double g, Double b);
  virtual ~Color3d	();
  Color3d 		(const Color3d& toCopy);  
  Color3d& operator = 	(const Color3d& toCopy);

  Color3d  operator +	(const Color3d& toAdd) const;
  Color3d& operator +=	(const Color3d& toAdd);

  Color3d  operator -	(const Color3d& toSub) const;
  Color3d& operator -=	(const Color3d& toSub);

  Color3d  operator *	(const Color3d& toMult) const;
  Color3d& operator *=	(const Color3d& toMult);

  Color3d  operator *	(Double scaleFactor) const;
  Color3d& operator *=	(Double scaleFactor);

  friend Color3d operator * (Double scaleFactor, const Color3d& c);

  Color3d  operator /	(Double scaleFactor) const;
  Color3d& operator /=	(Double scaleFactor);

  bool operator == 	(const Color3d& toCompare);
  bool operator != 	(const Color3d& toCompare);
  
  Double intensity		() const;
  void inline glLoad	();
  void clampTo			(Double min, Double max);
};

class MeshTriangle
{
public:
	MeshTriangle	();
	MeshTriangle	(unsigned int v0, unsigned int v1, unsigned int v2);
	~MeshTriangle	();

	unsigned int&	operator	[]	(int index);
private:
	unsigned int	indices[3];

	friend ostream& operator << (ostream& out, const MeshTriangle& toPrint);

};
#endif // GEOMETRY_H
