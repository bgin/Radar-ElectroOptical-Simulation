#ifndef __GMS_COORDINATE_POINTS_H__
#define __GMS_COORDINATE_POINTS_H__

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <cmath>
#include <initializer_list>
#include <immintrin.h>


namespace file_version {

    const unsigned int GMS_COORDINATE_POINTS_MAJOR = 1U;
    const unsigned int GMS_COORDINATE_POINTS_MINOR = 0U;
    const unsigned int GMS_COORDINATE_POINTS_MICRO = 0U;
    const unsigned int GMS_COORDINATE_POINTS_FULLVER =
      1000U*GMS_COORDINATE_POINTS_MAJOR+
      100U*GMS_COORDINATE_POINTS_MINOR+
      10U*GMS_COORDINATE_POINTS_MICRO;
    const char * const GMS_COORDINATE_POINTS_CREATION_DATE = "05-03-2016 09:39 AM +00200";
    const char * const GMS_COORDINATE_POINTS_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COORDINATE_POINTS_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COORDINATE_POINTS_DESCRIPTION   = "Coordinate points templated class."

}




namespace  gms {

          namespace math {


	/**********************************************
	   Cartesian 3D Corrdinate Points class
	   Primary template declaration - unimplemented.
	 **********************************************/
	template<typename T>  class CoordPoints;

	/**********************************************
	    float specialization for __m128 data type.

	***********************************************/

	template<>   class CoordPoints<__m128> {


	public:

		/********************************************
		      Constructors and Destructor
		*********************************************/

		/*
		@brief    Default Ctor initializes m_F32Coord3D to zero.
		*/
		CoordPoints()noexcept(true);

		/*
		@brief    Ctor initializes m_F32Coord3D object from user passed 
		          scalar float values: x, y, z
		*/
		CoordPoints(_In_ const float, _In_ const float, _In_ const float)noexcept(true);

		/*
		@brief     Ctor initializes m_F32Coord3D object from user passed 
		           scalar float value: x (implies that objects coordinates are equal)
		*/
		CoordPoints(_In_ const float)noexcept(true);

		/*
		@brief    explicit Ctor initializes m_F32Coord3D object from user passed 
		           scalar double precision values: x, y, z
				   Warning: Implied precision loss conversion.
		*/
		explicit CoordPoints(_In_ const double, _In_ const double, _In_ const double)noexcept(true);

		/*
		@brief     Ctor initializes m_F32Coord3D object from the user passed
		           reference to static array.Last element should be set NAN

		*/
		CoordPoints(_In_ float(&)[4])noexcept(true);

		/*
		@brief     Ctor initializes m_F32Coord3D object from the user
		           std::initializer_list<float>
		*/
		CoordPoints(_In_ const std::initializer_list<float> &);

		/*
		@brief     Ctor initializes m_F32Coord3D object from the user
		           passed __m128 union.
		*/
		CoordPoints(_In_ const __m128 &)noexcept(true);

		/*
		@brief     Copy-Ctor
		*/
		CoordPoints(_In_ const CoordPoints &);

		/*
		@brief     Move-Ctor.
		*/
		CoordPoints(_In_ CoordPoints &&);

		/*
		@brief     Dtor = default
		*/
		~CoordPoints()noexcept(true) = default;

		/****************************************************************
		        Member and friend operators
		*****************************************************************/

		/*
		@brief       *this = rhs (copy operation).
		*/
		inline    auto      operator=(_In_ const CoordPoints &)->CoordPoints &;

		/*
		@brief       *this = rhs (move operation).
		*/
		inline    auto      operator=(_In_  const CoordPoints &&)->CoordPoints &;

		/*
		@brief       *this += rhs
		*/
		inline    auto      operator+=(_In_ const CoordPoints &)->CoordPoints &;

		/*
		@brief       *this += rhs, where rhs is of type float.
		*/
		inline    auto      operator+=(_In_ const float)->CoordPoints &;

		/*
		@brief       *this -= rhs
		*/
		inline    auto      operator-=(_In_ const CoordPoints &)->CoordPoints &;

		/*
		@brief       *this -= rhs, where rhs is of type float.
		*/
		inline    auto      operator-=(_In_ const float)->CoordPoints &;

		/*
		@brief       *this *= rhs
		*/
		inline    auto      operator*=(_In_ const CoordPoints &)->CoordPoints &;

		/*
		@brief        *this *= rhs, where rhs is of type float.
		*/
		inline    auto      operator*=(_In_ const float)->CoordPoints &;

		/*
		@brief         *this /= rhs 
		*/
		inline   auto       operator/=(_In_ const CoordPoints &)->CoordPoints &;

		/*
		@brief         *this /= rhs, where rhs is of type float.
		*/
		inline   auto       operator/=(_In_ const float)->CoordPoints &;

		/*
		@brief          *this == rhs
		*/
		inline   auto       operator==(_In_ const CoordPoints &)const-> __m128;

		/*
		@brief         *this == rhs,  where rhs is of type float.
		*/
		inline   auto       operator==(_In_ const float)const-> __m128;

		/*
		@brief         *this != rhs
		*/
		inline   auto       operator!=(_In_ const CoordPoints &)const-> __m128;

		/*
		@brief         *this != rhs, where rhs is of type float
		*/
		inline   auto       operator!=(_In_ const float)const-> __m128;

		/*
		@brief         conversion operator from __m128 to m_F32Coord3D (mutable)
		*/
		operator  __m128 ();

		/*
		@brief        conversion operator from __m128 to m_F32Coord3D (immutable)
		*/
		operator  __m128 const () const;

		/*
		@brief        operator<< - writes to the ostream.
		*/
		friend auto inline  operator<<(_In_ std::ostream &, _In_ const CoordPoints &)->std::ostream &;

		/******************************************************************
		                     Accessor methods
		*******************************************************************/

		/*
		@brief       Returns m_F32Coord3D member.
		*/
		inline    auto     getF32Coord3D()const-> __m128;

		/*
		@brief       Returns X coordinate scalar value
		*/
		inline    auto     X()const-> float;

		/*
		@brief       Returns Y coordinate scalar value.
		*/
		inline    auto     Y()const-> float;

		/*
		@brief       Returns Z coordinate scalar value.
		*/
		inline    auto     Z()const-> float;

		/*
		@brief       Converts m_F32Coord3D data is mutable.
		*/
		inline    auto     CArray()->double *;

		/*
		@brief       Converts m_F32Coord3D data is immutable.
		*/
		inline    auto     CArray()const->double const *;


		/*
		@brief       Sets X corrdinate while keeping two others unchanged.
		*/
		inline    auto     setX(_In_ const float)->void;

		/*
		@brief       Sets Y coordinate while keeping two others unchanged.
		*/
		inline    auto     setY(_In_ const float)->void;

		/*
		@brief       Sets Z coordinate while keeping two others unchanged.
		*/
		inline    auto     setZ(_In_ const float)->void;

		/*****************************************************************
		                   Class member functions
		******************************************************************/



		/*
		@brief        Returns coordinate points magnitude.
		*/
		inline    auto      magnitude()const->float;

		/*
		@brief        Returns coordinate points perpendicular.
		*/
		inline    auto      perpendicular()const->float;

		/*
		@brief        Returns coordinate points conversion to Rho angle.
		*/
		inline   auto       rho()const->float;

		/*
		@brief       Returns coordinate points conversion to theta angle.
		*/
		inline   auto       theta()const->float;

		/*
		@brief       Returns coordinate points conversion to phi angle.
		*/
		inline   auto        phi()const->float;


	private:

		/*
		@brief    class member variable m_F32Coord3D
		*/
		__m128   m_F32Coord3D;
	};




#include "GMS_coordinate_points.inl"
    }
   
}
#endif /*__GMS_COORDINATE_POINTS_H__*/
