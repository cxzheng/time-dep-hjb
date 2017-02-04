/*
 * =====================================================================================
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  vector2.hpp
 *
 *        Version:  1.0
 *        Created:  6/12/13 10:44:29
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
 
#ifndef VECMATH_VECTOR2_H
#   define VECMATH_VECTOR2_H

#include <assert.h>
#include <iostream>

/*!
 * Class for two dimensional vector.
 */
template <typename T> 
class vector2 
{
    public:
        static const vector2<T> ZERO;

        typedef T   element;

        union 
        {
            T x;    // first element of vector, alias for X-coordinate.
            T s;    // alias for S-coordinate (texture notation)
            T u;    // alias u, for fluid simulator
        };
 
        union 
        {
            T y;    // second element of vector, alias for Y-coordinate.
            T t;    // alias for T-coordinate(texture notation)
            T v;    // third alias v, for fluid simulator
        };

        vector2():x(0),y(0) { }
        vector2(T nx, T ny):x(nx), y(ny) { }
        /*! Copy constructor */
        vector2(const vector2<T>& src):x(src.x), y(src.y) { }
        /*! Copy casting constructor. */
        template <typename FromT> 
        vector2(const vector2<FromT>& src):
                x(static_cast<T>(src.x)), 
                y(static_cast<T>(src.y)) 
        { }

        /*! Directly set the fields */
        template <typename FromT>
        void set(FromT x, FromT y)
        {
            this->x = static_cast<T>(x);
            this->y = static_cast<T>(y);
        }

        void zero()
        {
            x = static_cast<T>(0);
            y = static_cast<T>(0);
        }

        /*!
         * provide the similar method with j3d.vecmath.Tuple2d
         * this = s*t1 + t2
         */
        void scaleAdd(T s, const vector2<T>& t1, const vector2<T>& t2)
        {
            x = s*t1.x + t2.x;
            y = s*t1.y + t2.y;
        }

        /*!
         * this = this + s*t1
         */
        void scaleAdd(T s, const vector2<T>& t1)
        {
            x += s*t1.x;
            y += s*t1.y;
        }

        /*!
         * this = s*this + t1
         */
        void selfScaleAdd(T s, const vector2<T>& t1)
        {
            x = x*s + t1.x;
            y = y*s + t1.y;
        }

        void clamp(const vector2<T>& minV, const vector2<T>& maxV)
        {
            x = x <= minV.x ? minV.x : (x > maxV.x ? maxV.x : x);
            y = y <= minV.y ? minV.y : (y > maxV.y ? maxV.y : y);
        }
         
        //================= operators ==================
        /*! Copy casting operator */
        template <typename FromT>
        vector2<T>& operator=(const vector2<FromT>& rhs)
        {
            x = static_cast<T>(rhs.x);
            y = static_cast<T>(rhs.y);
            return *this;
        }

        /*! Copy operator */
        vector2<T>& operator=(const vector2<T>& rhs)
        {
            x = rhs.x;
            y = rhs.y;
            return * this;
        }

        /*!
         * Array access operator
         * @param n Array index
         * @return For n = 0, reference to x coordinate, else reference to y 
         * y coordinate.
         */
        T& operator[](int n)
        {
            assert(n >= 0 && n <= 1);
            //return n == 0 ? x : y;
            return ((T*)this)[n];
        }
 
        /*! Addition operator */
        vector2<T> operator + (const vector2<T>& rhs) const
        {
            return vector2<T> (x + rhs.x, y + rhs.y);
        }

        /*! Substraction operator */
        vector2<T> operator - (const vector2<T>& rhs) const
        {
            return vector2<T>(x - rhs.x, y - rhs.y);
        }

        vector2<T> operator - () const
        {
            return vector2<T>(-x, -y);
        }

        /*! Multiplication operator */
        vector2<T> operator * (const vector2<T>& rhs) const
        {
            return vector2<T> (x * rhs.x, y * rhs.y);
        }

        /*! Division operator */
        vector2<T> operator / (const vector2<T>& rhs) const 
        {
            return vector2<T> (x / rhs.x, y / rhs.y);
        }

        /*! Addition operator */
        vector2<T>& operator += (const vector2<T>& rhs)
        {
            x += rhs.x;
            y += rhs.y;
            return * this;
        }

        /*! Substraction operator */
        vector2<T>& operator -= (const vector2<T>& rhs)
        {
            x -= rhs.x;
            y -= rhs.y;
            return *this;
        }

        template <typename FromT>
        vector2<T>& operator -= (const vector2<FromT>& rhs)
        {
            x -= static_cast<T>(rhs.x);
            y -= static_cast<T>(rhs.y);
            return  *this;
        }
        
        /*! Multiplication operator */
        vector2<T>& operator *= (const vector2<T>& rhs)
        {
            x *= rhs.x;
            y *= rhs.y;
            return * this;
        }
        
        /*! Division operator */
        vector2<T>& operator /= (const vector2<T>& rhs)
        {
            x /= rhs.x;
            y /= rhs.y;
            return * this;
        }

        //============== scalar vector operator ===============
        /*! Addition operator */
        vector2<T> operator + (T rhs) const
        {
            return vector2<T> (x + rhs, y + rhs);
        }
        
        /*! Substraction operator */
        vector2<T> operator - (T rhs) const
        {
            return vector2<T> (x - rhs, y - rhs);
        }
        
        /*! Multiplication operator */
        vector2<T> operator * (T rhs) const
        {
            return vector2<T>(x * rhs, y * rhs);
        }
        
        friend vector2<T> operator * (T lhs, const vector2<T>& rhs)
        {
            return vector2<T>(lhs*rhs.x, lhs*rhs.y);
        }

        /*! Division operator */
        vector2<T> operator / (T rhs) const
        {
            return vector2<T> (x / rhs, y / rhs);
        }

        /*! Addition operator */
        vector2<T>& operator += (T rhs)
        {
            x += rhs;
            y += rhs;
            return * this;
        }

        /*! Substraction operator */
        vector2<T>& operator -= (T rhs) 
        {
            x -= rhs;
            y -= rhs;
            return * this;
        }

        /*! Multiplication operator */
        vector2<T>& operator *= (T rhs)
        {
            x *= rhs;
            y *= rhs;
            return * this;
        }

        /*! Division operator */
        vector2<T>& operator /= (T rhs) 
        {
            x /= rhs;
            y /= rhs;
            return * this;
        }

        ///*!
        // * Equality test operator
        // * @param rhs Right hand side argument of binary operator.
        // * @note Test of equality is based of threshold _epsilon value. To be two
        // * values equal, must satisfy this condition | lws.x - rhs.y | < _epsilon,
        // * same for y-coordinate.
        // */
        //bool operator == (const vector2<T>& rhs) const 
        //{
        //    return std::abs(x - rhs.x) < _epsilon && std::abs(y - rhs.y) < _epsilon;
        //}
        ///*! Inequality test operator */
        //bool operator!=(const vector2<T>& rhs) const 
        //{
        //    return std::abs(x - rhs.x) >= _epsilon || std::abs(y - rhs.y) >= _epsilon;
        //}


        bool equals(const vector2<T>& rhs) const
        {
            return x == rhs.x && y == rhs.y;
        }
        
        /*! Dot product of two vectors */
        T dotProduct(const vector2<T>& rhs) const
        {
            return x*rhs.x + y*rhs.y;
        }

        /*! Cross product opertor */    
        T crossProduct(const vector2<T>& rhs) const
        {
            return x*rhs.y - rhs.x*y;
        }

        /*! Get lenght of vector. */
        T length() const
        {
            return (T)sqrt(x * x + y * y);
        }

        /**
         * Normalize vector
         */
        void normalize()
        {
            T s = length();
            x /= s;
            y /= s;
        }

        /**
         * Normalize vector and return its original length
         */
        T normalize2()
        {
            T s = length();
            x /= s;
            y /= s;
            return s;
        }

        /*! Return square of length. */
        T lengthSqr() const
        {
            return x * x + y * y;
        }

        /*! distance between two vectors */
        T distance(const vector2<T>& v) const
        {
            T dx = x - v.x;
            T dy = y - v.y;
            return (T)sqrt(dx*dx + dy*dy);
        }
        
        T distanceSqr(const vector2<T>& v) const
        {
            T dx = x - v.x;
            T dy = y - v.y;
            return dx*dx + dy*dy;
        }

        /*!
         * Linear interpolation of two vectors
         * @param fact Factor of interpolation. For translation from positon
         * of this vector to vector r, values of factor goes from 0.0 to 1.0.
         * @param r Second vector for interpolation 
         * @note Hovewer values of fact parameter are reasonable only in interval
         * [0.0 , 1.0], you can pass also values outside of this interval and you
         * can get result (extrapolation?)
         */
        vector2<T> lerp(T fact, const vector2<T>& r) const
        {
            return (*this) + (r - (*this)) * fact;    
        }
        
        
        /*!
         * Conversion to pointer operator
         * @return Pointer to internaly stored (in managment of class vector2<T>)
         *         used for passing vector2<T> values to gl*2[fd] functions.
         */
        operator T*() { return (T*) this; }
        /**
         * Conversion to pointer operator
         * @return Constant Pointer to internaly stored (in managment of class vector2<T>)
         * used for passing vector2<T> values to gl*2[fd] functions.
         */
        operator const T*() const { return (const T*) this; }

        /*! Output to stream operator */
        friend std::ostream& operator<<(std::ostream& lhs, const vector2<T>& rhs) 
        {
            lhs << "[" << rhs.x << "," << rhs.y << "]";
            return lhs;
        }

        /*! input from a stream */
        friend std::istream& operator>>(std::istream& in, vector2<T>& obj)
        {
            in >> obj.x >> obj.y;
            return in;
        }
};

template <typename T>
const vector2<T> vector2<T>::ZERO(0, 0, 0);

typedef class vector2 <float>   vector2f;
typedef class vector2 <double>  vector2d;
typedef class vector2 <int>     vector2i;

#endif


