/*
-----------------------------------------------------------------------------
This source file is part of OGRE
    (Object-oriented Graphics Rendering Engine)
For the latest info, see http://www.ogre3d.org/

Copyright (c) 2000-2006 Torus Knot Software Ltd
Also see acknowledgements in Readme.html

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place - Suite 330, Boston, MA 02111-1307, USA, or go to
http://www.gnu.org/copyleft/lesser.txt.

You may alternatively use this source under the terms of a specific version of
the OGRE Unrestricted License provided you have obtained such a license from
Torus Knot Software Ltd.
-----------------------------------------------------------------------------
*/
#ifndef __Vector2_H__
#define __Vector2_H__


#include "OgreCommon.h"

namespace Ogre
{

    /** Standard 2-dimensional vector.
        @remarks
            A direction in 2D space represented as distances along the 2
            orthoganal axes (x, y). Note that positions, directions and
            scaling factors can be represented by a vector, depending on how
            you interpret the values.
    */
    class _OgreExport Vector2
    {
    public:
        float x, y;

    public:
        inline Vector2()
        {
        }

        inline Vector2(const float fX, const float fY )
            : x( fX ), y( fY )
        {
        }

        inline explicit Vector2( const float scaler )
            : x( scaler), y( scaler )
        {
        }

        inline explicit Vector2( const float afCoordinate[2] )
            : x( afCoordinate[0] ),
              y( afCoordinate[1] )
        {
        }

        inline explicit Vector2( const int afCoordinate[2] )
        {
            x = (float)afCoordinate[0];
            y = (float)afCoordinate[1];
        }

        inline explicit Vector2( float* const r )
            : x( r[0] ), y( r[1] )
        {
        }

        inline Vector2( const Vector2& rkVector )
            : x( rkVector.x ), y( rkVector.y )
        {
        }

		inline float operator [] ( const size_t i ) const
        {
            assert( i < 2 );

            return *(&x+i);
        }

		inline float& operator [] ( const size_t i )
        {
            assert( i < 2 );

            return *(&x+i);
        }

		/// Pointer accessor for direct copying
		inline float* ptr()
		{
			return &x;
		}
		/// Pointer accessor for direct copying
		inline const float* ptr() const
		{
			return &x;
		}

        /** Assigns the value of the other vector.
            @param
                rkVector The other vector
        */
        inline Vector2& operator = ( const Vector2& rkVector )
        {
            x = rkVector.x;
            y = rkVector.y;

            return *this;
        }

		inline Vector2& operator = ( const float fScalar)
		{
			x = fScalar;
			y = fScalar;

			return *this;
		}

        inline bool operator == ( const Vector2& rkVector ) const
        {
            return ( x == rkVector.x && y == rkVector.y );
        }

        inline bool operator != ( const Vector2& rkVector ) const
        {
            return ( x != rkVector.x || y != rkVector.y  );
        }

        // arithmetic operations
        inline Vector2 operator + ( const Vector2& rkVector ) const
        {
            return Vector2(
                x + rkVector.x,
                y + rkVector.y);
        }

        inline Vector2 operator - ( const Vector2& rkVector ) const
        {
            return Vector2(
                x - rkVector.x,
                y - rkVector.y);
        }

        inline Vector2 operator * ( const float fScalar ) const
        {
            return Vector2(
                x * fScalar,
                y * fScalar);
        }

        inline Vector2 operator * ( const Vector2& rhs) const
        {
            return Vector2(
                x * rhs.x,
                y * rhs.y);
        }

        inline Vector2 operator / ( const float fScalar ) const
        {
            assert( fScalar != 0.0f );

            float fInv = 1.0f / fScalar;

            return Vector2(
                x * fInv,
                y * fInv);
        }

        inline Vector2 operator / ( const Vector2& rhs) const
        {
            return Vector2(
                x / rhs.x,
                y / rhs.y);
        }

        inline const Vector2& operator + () const
        {
            return *this;
        }

        inline Vector2 operator - () const
        {
            return Vector2(-x, -y);
        }

        // overloaded operators to help Vector2
        inline friend Vector2 operator * ( const float fScalar, const Vector2& rkVector )
        {
            return Vector2(
                fScalar * rkVector.x,
                fScalar * rkVector.y);
        }

        inline friend Vector2 operator / ( const float fScalar, const Vector2& rkVector )
        {
            return Vector2(
                fScalar / rkVector.x,
                fScalar / rkVector.y);
        }

        inline friend Vector2 operator + (const Vector2& lhs, const float rhs)
        {
            return Vector2(
                lhs.x + rhs,
                lhs.y + rhs);
        }

        inline friend Vector2 operator + (const float lhs, const Vector2& rhs)
        {
            return Vector2(
                lhs + rhs.x,
                lhs + rhs.y);
        }

        inline friend Vector2 operator - (const Vector2& lhs, const float rhs)
        {
            return Vector2(
                lhs.x - rhs,
                lhs.y - rhs);
        }

        inline friend Vector2 operator - (const float lhs, const Vector2& rhs)
        {
            return Vector2(
                lhs - rhs.x,
                lhs - rhs.y);
        }
        // arithmetic updates
        inline Vector2& operator += ( const Vector2& rkVector )
        {
            x += rkVector.x;
            y += rkVector.y;

            return *this;
        }

        inline Vector2& operator += ( const float fScaler )
        {
            x += fScaler;
            y += fScaler;

            return *this;
        }

        inline Vector2& operator -= ( const Vector2& rkVector )
        {
            x -= rkVector.x;
            y -= rkVector.y;

            return *this;
        }

        inline Vector2& operator -= ( const float fScaler )
        {
            x -= fScaler;
            y -= fScaler;

            return *this;
        }

        inline Vector2& operator *= ( const float fScalar )
        {
            x *= fScalar;
            y *= fScalar;

            return *this;
        }

        inline Vector2& operator *= ( const Vector2& rkVector )
        {
            x *= rkVector.x;
            y *= rkVector.y;

            return *this;
        }

        inline Vector2& operator /= ( const float fScalar )
        {
            assert( fScalar != 0.0f );

            float fInv = 1.0f / fScalar;

            x *= fInv;
            y *= fInv;

            return *this;
        }

        inline Vector2& operator /= ( const Vector2& rkVector )
        {
            x /= rkVector.x;
            y /= rkVector.y;

            return *this;
        }

        /** Returns the length (magnitude) of the vector.
            @warning
                This operation requires a square root and is expensive in
                terms of CPU operations. If you don't need to know the exact
                length (e.g. for just comparing lengths) use squaredLength()
                instead.
        */
        inline float length () const
        {
            return sqrtf( x * x + y * y );
        }

        /** Returns the square of the length(magnitude) of the vector.
            @remarks
                This  method is for efficiency - calculating the actual
                length of a vector requires a square root, which is expensive
                in terms of the operations required. This method returns the
                square of the length of the vector, i.e. the same as the
                length but before the square root is taken. Use this if you
                want to find the longest / shortest vector without incurring
                the square root.
        */
        inline float squaredLength () const
        {
            return x * x + y * y;
        }

        /** Calculates the dot (scalar) product of this vector with another.
            @remarks
                The dot product can be used to calculate the angle between 2
                vectors. If both are unit vectors, the dot product is the
                cosine of the angle; otherwise the dot product must be
                divided by the product of the lengths of both vectors to get
                the cosine of the angle. This result can further be used to
                calculate the distance of a point from a plane.
            @param
                vec Vector with which to calculate the dot product (together
                with this one).
            @returns
                A float representing the dot product value.
        */
        inline float dotProduct(const Vector2& vec) const
        {
            return x * vec.x + y * vec.y;
        }

        /** Normalises the vector.
            @remarks
                This method normalises the vector such that it's
                length / magnitude is 1. The result is called a unit vector.
            @note
                This function will not crash for zero-sized vectors, but there
                will be no changes made to their components.
            @returns The previous length of the vector.
        */
        inline float normalise()
        {
            float fLength = sqrtf( x * x + y * y);

            // Will also work for zero-sized vectors, but will change nothing
            if ( fLength > 1e-08 )
            {
                float fInvLength = 1.0f / fLength;
                x *= fInvLength;
                y *= fInvLength;
            }

            return fLength;
        }



        /** Returns a vector at a point half way between this and the passed
            in vector.
        */
        inline Vector2 midPoint( const Vector2& vec ) const
        {
            return Vector2(
                ( x + vec.x ) * 0.5f,
                ( y + vec.y ) * 0.5f );
        }

        /** Returns true if the vector's scalar components are all greater
            that the ones of the vector it is compared against.
        */
        inline bool operator < ( const Vector2& rhs ) const
        {
            if( x < rhs.x && y < rhs.y )
                return true;
            return false;
        }

        /** Returns true if the vector's scalar components are all smaller
            that the ones of the vector it is compared against.
        */
        inline bool operator > ( const Vector2& rhs ) const
        {
            if( x > rhs.x && y > rhs.y )
                return true;
            return false;
        }

        /** Sets this vector's components to the minimum of its own and the
            ones of the passed in vector.
            @remarks
                'Minimum' in this case means the combination of the lowest
                value of x, y and z from both vectors. Lowest is taken just
                numerically, not magnitude, so -1 < 0.
        */
        inline void makeFloor( const Vector2& cmp )
        {
            if( cmp.x < x ) x = cmp.x;
            if( cmp.y < y ) y = cmp.y;
        }

        /** Sets this vector's components to the maximum of its own and the
            ones of the passed in vector.
            @remarks
                'Maximum' in this case means the combination of the highest
                value of x, y and z from both vectors. Highest is taken just
                numerically, not magnitude, so 1 > -3.
        */
        inline void makeCeil( const Vector2& cmp )
        {
            if( cmp.x > x ) x = cmp.x;
            if( cmp.y > y ) y = cmp.y;
        }

        /** Generates a vector perpendicular to this vector (eg an 'up' vector).
            @remarks
                This method will return a vector which is perpendicular to this
                vector. There are an infinite number of possibilities but this
                method will guarantee to generate one of them. If you need more
                control you should use the Quaternion class.
        */
        inline Vector2 perpendicular(void) const
        {
            return Vector2 (-y, x);
        }
        /** Calculates the 2 dimensional cross-product of 2 vectors, which results
			in a single floating point value which is 2 times the area of the triangle.
        */
        inline float crossProduct( const Vector2& rkVector ) const
        {
            return x * rkVector.y - y * rkVector.x;
        }

        /** Returns true if this vector is zero length. */
        inline bool isZeroLength(void) const
        {
            float sqlen = (x * x) + (y * y);
            return (sqlen < (1e-06 * 1e-06));
        }

        /** As normalise, except that this vector is unaffected and the
            normalised vector is returned as a copy. */
        inline Vector2 normalisedCopy(void) const
        {
            Vector2 ret = *this;
            ret.normalise();
            return ret;
        }

        /** Calculates a reflection vector to the plane with the given normal .
        @remarks NB assumes 'this' is pointing AWAY FROM the plane, invert if it is not.
        */
        inline Vector2 reflect(const Vector2& normal) const
        {
            return Vector2( *this - ( 2 * this->dotProduct(normal) * normal ) );
        }

        // special points
        static const Vector2 ZERO;
        static const Vector2 UNIT_X;
        static const Vector2 UNIT_Y;
        static const Vector2 NEGATIVE_UNIT_X;
        static const Vector2 NEGATIVE_UNIT_Y;
        static const Vector2 UNIT_SCALE;

        /** Function for writing to a stream.
        */
        inline _OgreExport friend std::ostream& operator <<
            ( std::ostream& o, const Vector2& v )
        {
            o << "Vector2(" << v.x << ", " << v.y <<  ")";
            return o;
        }

    };

}
#endif
