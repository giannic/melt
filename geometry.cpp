#include "GL/glut.h"
//#include "common.h"
#include "geometry.h"

ostream& operator << (ostream& out, const Tuple3d& toPrint)
{
	out << toPrint[0] << ", " << toPrint[1] << ", " << toPrint[2];
	return out;
}

Tuple3d::Tuple3d ()
{
}

Tuple3d::Tuple3d(const Tuple3d& toCopy)
{
	data[0] = toCopy[0];
	data[1] = toCopy[1];
	data[2] = toCopy[2];
}

Tuple3d::~Tuple3d ()
{
}

Tuple3d& Tuple3d::operator = (const Tuple3d& toCopy)
{
	data[0] = toCopy[0];
	data[1] = toCopy[1];
	data[2] = toCopy[2];

  return *this;
}

Point3d::Point3d ()
{
}

Point3d::Point3d (Double x, Double y, Double z)
{
    data[0] = x;
    data[1] = y;
    data[2] = z;
}

Point3d::~Point3d ()
{
}

Point3d::Point3d (const Point3d& toCopy)
{
	data[0] = toCopy[0];
	data[1] = toCopy[1];
	data[2] = toCopy[2];
}


Point3d& Point3d::operator = (const Point3d& toCopy)
{
	data[0] = toCopy[0];
	data[1] = toCopy[1];
	data[2] = toCopy[2];

	return *this;
}

bool Point3d::operator == (const Point3d& compare)
{
  for (int i = 0; i < 3; ++i)
    if (data[i] != compare[i])
      return false;
  
  return true;
}

bool Point3d::operator != (const Point3d& compare)
{
  for (int i = 0; i < 3; ++i)
    if (data[i] != compare[i])
      return true;

  return false;
}

Point3d Point3d::operator + (const Vector3d& toAdd) const
{
  Point3d sum = *this;
  sum += toAdd;
  return sum;
}

Point3d& Point3d::operator += (const Vector3d& toAdd)
{
  for (int i = 0; i < 3; ++i)
    data[i] += toAdd[i];
  return *this;
}

Point3d Point3d::operator * (const Point3d& toMult) const
{
  Point3d product = *this;
  product *= toMult;
  return product;
}

Point3d& Point3d::operator *= (const Point3d& toMult)
{  
  for (int i = 0; i < 3; ++i)
    data[i] *= toMult[i];
  return *this;
}

Point3d Point3d::operator / (const Point3d& toDiv) const
{
  Point3d quotient = *this;
  quotient /= toDiv;
  return quotient;
}

Point3d& Point3d::operator /= (const Point3d& toDiv)
{
  for (int i = 0; i < 3; ++i)
    data[i] /= toDiv[i];
  return *this;
}

Point3d Point3d::operator * (Double scaleFactor) const
{
  Point3d product = *this;
  product *= scaleFactor;
  return product;
}

Point3d& Point3d::operator *= (Double scaleFactor)
{
  for (int i = 0; i < 3; ++i)
    data[i] *= scaleFactor;
  return *this;
}

Point3d Point3d::operator - (const Vector3d& toSub)
{
  Point3d diff = *this;
  diff -= toSub;
  return diff;
}

Point3d& Point3d::operator -= (const Vector3d& toSub)
{
  for (int i = 0; i < 3; ++i)
    data[i] -= toSub[i];
  return *this;
}

Point3d operator - (const Vector3d& v, const Point3d& p)
{
  Point3d diff;
  for (int i = 0; i < 3; ++i)
    diff[i] = v[i] - p[i];
  return diff;
}

void Point3d::glLoad()
{
	glVertex3dv(data);
}

Vector3d::Vector3d ()
{
}

Vector3d::Vector3d (Double x, Double y, Double z)
{
	data[0] = x;
	data[1] = y;
	data[2] = z;
}

Vector3d::~Vector3d ()
{
}

Vector3d::Vector3d (const Vector3d& toCopy)
{
	data[0] = toCopy[0];
	data[1] = toCopy[1];
	data[2] = toCopy[2];
}

Vector3d::Vector3d(const Point3d& dest)
{
	data[0] = dest[0];
	data[1] = dest[1];
	data[2] = dest[2];
}

Vector3d::Vector3d(const Point3d& start, const Point3d& end)
{
	data[0] = end[0] - start[0];
	data[1] = end[1] - start[1];
	data[2] = end[2] - start[2];
}

Vector3d& Vector3d::operator = (const Point3d& dest)
{
  Tuple3d::operator= (dest);

  return *this;
}

Vector3d& Vector3d::operator = (const Vector3d& toCopy)
{
  Tuple3d::operator= (toCopy);

  return *this;
}

Vector3d Vector3d::operator + (const Vector3d& toAdd) const
{
  Vector3d sum = *this;
  sum += toAdd;
  return sum;
}

Vector3d& Vector3d::operator += (const Vector3d& toAdd)
{
  for (int i = 0; i < 3; ++i)
    data[i] += toAdd[i];
  return *this;
}

Vector3d Vector3d::operator - () const
{
  Vector3d neg = *this;
  for (int i = 0; i < 3; ++i)
    neg[i] *= -1;
  return neg;
}

Vector3d Vector3d::operator - (const Vector3d& toSub) const
{
  Vector3d diff = *this;
  diff -= toSub;
  return diff;
}

Vector3d& Vector3d::operator -= (const Vector3d& toSub)
{
  for (int i = 0; i < 3; ++i)
    data[i] -= toSub[i];
  return *this;
}

Vector3d Vector3d::operator * (const Vector3d& toMult) const
{
  Vector3d product = *this;
  product *= toMult;
  return product;
}

Vector3d& Vector3d::operator *= (const Vector3d& toMult)
{  
  for (int i = 0; i < 3; ++i)
    data[i] *= toMult[i];
  return *this;
}

Vector3d Vector3d::operator / (const Vector3d& toDiv) const
{
  Vector3d quotient = *this;
  quotient /= toDiv;
  return quotient;
}

Vector3d& Vector3d::operator /= (const Vector3d& toDiv)
{
  for (int i = 0; i < 3; ++i)
    data[i] /= toDiv[i];
  return *this;
}

Vector3d Vector3d::operator * (Double scaleFactor) const
{
  Vector3d product = *this;
  product *= scaleFactor;
  return product;
}

Vector3d& Vector3d::operator *= (Double scaleFactor)
{
  for (int i = 0; i < 3; ++i)
    data[i] *= scaleFactor;
  return *this;
}

Vector3d operator * (Double scaleFactor, const Vector3d& v)
{
  return v * scaleFactor;
}

Vector3d Vector3d::operator / (Double scaleFactor) const
{
  Vector3d quotient = *this;
  quotient /= scaleFactor;
  return quotient;
}

Vector3d& Vector3d::operator /= (Double scaleFactor)
{
  if (scaleFactor == 0)
    {
      cerr << "cannot divide by zero " << endl;
      return *this;
    }
  for (int i = 0; i < 3; ++i)
    data[i] /= scaleFactor;
  return *this;
}

Double Vector3d::dot (const Vector3d& toDot)
{
  Double dotProduct = 0;
  for (int i = 0; i < 3; ++i)
    dotProduct += data[i] * toDot[i];
  return dotProduct;
}

/*
 * Take the cross product of two 3-dimensional vectors. Note that if called 
 * on anything but 3-dimensional vectors a compiler error will result.
 */
Vector3d Vector3d::cross (const Vector3d& toCross)
{
  Vector3d crossProduct;
  crossProduct[0] = data[1] * toCross[2] - data[2] * toCross[1];
  crossProduct[1] = data[2] * toCross[0] - data[0] * toCross[2];
  crossProduct[2] = data[0] * toCross[1] - data[1] * toCross[0];
  return crossProduct;
}

Vector3d& Vector3d::normalize ()
{
  Double len = this->length();

  if (len == 0)
    {
      cerr << "Cannot normalize 0-length vector " << endl;
      return *this;
    }

  for (int i = 0; i < 3; ++i)
    data[i] /= len;
  
  return *this;
}

Vector3d Vector3d::getUnit ()
{
  return Vector3d(this->normalize());;
}

Double Vector3d::length()
{
  Double sumSquares = 0;
  for (int i = 0; i < 3; ++i)
    sumSquares += data[i] * data[i];

  return sqrt(sumSquares);
}

void Vector3d::glLoad ()
{
    glNormal3dv(data);
}

Color3d::Color3d ()
: Tuple3d()
{
}

Color3d::Color3d (Double r, Double g, Double b)
: Tuple3d()
{
	data[0] = r;
	data[1] = g;
	data[2] = b;
}

Color3d::~Color3d ()
{
}

Color3d::Color3d (const Color3d& toCopy)
: Tuple3d(toCopy)
{
}

Color3d& Color3d::operator = (const Color3d& toCopy)
{
  if (this == &toCopy)
    return *this;

  Tuple3d::operator= (toCopy);

  return *this;
}

Color3d Color3d::operator + (const Color3d& toAdd) const
{
  Color3d sum = *this;
  sum += toAdd;
  return sum;
}

Color3d& Color3d::operator += (const Color3d& toAdd)
{  
  for (int i = 0; i < 3; ++i)
    data[i] += toAdd[i];
  return *this;
}

Color3d Color3d::operator - (const Color3d& toSub) const
{
  Color3d diff = *this;
  diff -= toSub;
  return diff;
}

Color3d& Color3d::operator -= (const Color3d& toSub)
{  
  for (int i = 0; i < 3; ++i)
    data[i] -= toSub[i];
  return *this;
}

Color3d Color3d::operator * (const Color3d& toMult) const
{
  Color3d product = *this;
  product *= toMult;
  return product;
}

Color3d& Color3d::operator *= (const Color3d& toMult)
{  
  for (int i = 0; i < 3; ++i)
    data[i] *= toMult[i];
  return *this;
}

Color3d Color3d::operator / (Double scaleFactor) const
{
  Color3d product = *this;
  product /= scaleFactor;
  return product;
}

Color3d& Color3d::operator /= (Double scaleFactor)
{
  if (scaleFactor == 0)
    {
      cerr << "cannot divide by 0" << endl;
      return *this;
    }
  for (int i = 0; i < 3; ++i)
    data[i] /= scaleFactor;
  return *this;
}


Color3d Color3d::operator * (Double scaleFactor) const
{
  Color3d product = *this;
  product *= scaleFactor;
  return product;
}

Color3d& Color3d::operator *= (Double scaleFactor)
{
  for (int i = 0; i < 3; ++i)
    data[i] *= scaleFactor;
  return *this;
}

Color3d operator * (Double scaleFactor, const Color3d& c)
{
  return c * scaleFactor;
}

bool Color3d::operator == (const Color3d& compare)
{
  for (int i = 0; i < 3; ++i)
    if (data[i] != compare[i])
      return false;
  
  return true;
}

bool Color3d::operator != (const Color3d& compare)
{
  for (int i = 0; i < 3; ++i)
    if (data[i] != compare[i])
      return true;

  return false;
}

Double Color3d::intensity () const
{
  Double sum = 0;
  for (int i = 0; i < 3; ++i)
    sum += data[i] * data[i];

  return (sum / 3);
}

void Color3d::glLoad ()
{
    glColor3dv(data);
}

void Color3d::clampTo(Double min, Double max)
{
  for (int i = 0; i < 3; ++i)
    {
      if (data[i] < min)
	data[i] = min;
      else if (data[i] > max)
	data[i] = max;
    }  
}

MeshTriangle::MeshTriangle ()
{
}

MeshTriangle::MeshTriangle (unsigned int v0, unsigned int v1, unsigned int v2)
{
	indices[0] = v0;
	indices[1] = v1;
	indices[2] = v2;
}

MeshTriangle::~MeshTriangle ()
{
}

unsigned int& MeshTriangle::operator [] (int index)
{
	return indices[index];
}

ostream& operator << (ostream& out, const MeshTriangle& toPrint)
{
	int x = toPrint.indices[0] + 1;
	int y = toPrint.indices[1] + 1;
	int z = toPrint.indices[2] + 1;

	out << x << " " 
		<< y << " " 
		<< z ;
	return out;
}