#ifndef __GMS_CARTESIAN_SYSTEM_H__
#define __GMS_CARTESIAN_SYSTEM_H__



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
#include <array>
#include <immintrin.h>
#include <string>
#include "GMS_coord_axis.h"

namespace file_version {

    const unsigned int GMS_CARTESIAN_SYSTEM_MAJOR = 1U;
    const unsigned int GMS_CARTESIAN_SYSTEM_MINOR = 0U;
    const unsigned int GMS_CARTESIAN_SYSTEM_MICRO = 0U;
    const unsigned int GMS_CARTESIAN_SYSTEM_FULLVER =
      1000U*GMS_CARTESIAN_SYSTEM_MAJOR+
      100U*GMS_CARTESIAN_SYSTEM_MINOR+
      10U*GMS_CARTESIAN_SYSTEM_MICRO;
    const char * const GMS_CARTESIAN_SYSTEM_CREATION_DATE = "25-02-2016 09:39 AM +00200";
    const char * const GMS_CARTESIAN_SYSTEM_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CARTESIAN_SYSTEM_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CARTESIAN_SYSTEM_DESCRIPTION   = "Cartesian system templated class."

}



namespace    gms {

        namespace math {


	/*
	
	                 Cartesian coordinate system templated on the Axis parameter.
	 
	*/

	template<class Coords>  class Cartesian3D {

	public:

		/*
		@brief     Default Ctor = delete
		*/
		Cartesian3D() = delete;

		/*
		@brief     "Main" class Ctor constructs Cartesian 3D Coordinate System.
		*/
		Cartesian3D(_In_ const Coords &, _In_ const Coords &,
			_In_ const Coords &, const std::string &, const bool);

		/*
		@brief     explicit Ctor initializes *this from the different Axis Ctor.
		           OtherCoords must have a members of type similar to
				   types which Axis has: Position,Orientation.
		*/
		template<class OtherCoords> explicit Cartesian3D(_In_ const OtherCoords &, _In_ const OtherCoords &,
			_In_ const OtherCoords &, const std::string &, const bool);

		/*
		@brief     Class Copy-Ctor.
		*/
		Cartesian3D(_In_ const Cartesian3D &);

		/*
		@brief     Class Move-Ctor.
		*/
		Cartesian3D(_In_ Cartesian3D &&);

		/*
		@brief      Dtor = default
		*/
		~Cartesian3D() = default;

		/*

		      ***   Class member operators   ***
		*/

		/*
		@brief      *this = rhs (copy operation).
		*/
		inline     auto      operator=(_In_ const Cartesian3D &)->Cartesian3D &;

		/*
		@brief      *this = const OtherCoords &rhs
		             OtherCoords must implement Axis objects.
		*/
	    template<class OtherCoords>  inline  auto operator=(_In_ const OtherCoords &)->Cartesian3D &;

		/*
		@brief       *this = rhs (move operation).
		*/
		inline     auto      operator=(_In_ Cartesian3D &&)->Cartesian3D &;

		/*
		@brief       *this == rhs 
		@params       implicit *this
		@params       const Cartesian3D &
		@returns      std::array<__m256d,6U> which contains __m256d data types initilized 
		              by the call to _mm256_cmp_pd intrinsic.
		*/
		inline     auto      operator==(_In_ Cartesian3D &)->std::array<__m256d,6U>;

		/*
		@brief       *this != rhs
		@params       implicit *this
		@params       const Cartesian3D &
		@returns      std::array<__m256d,6U> which contains __m256d data types initilized
		              by the call to _mm256_cmp_pd intrinsic.
		*/
		inline     auto       operator!=(_In_ Cartesian3D &)->std::array<__m256d, 6U>;

		

		/*
		 
		                ***      Accessors     ***
		*/

		/*
		"brief     Returns mX_Axis class member object
		*/
		inline      auto      getX_Axis() const->Coords;

		/*
		@brief     Returns mY_Axis class member object.
		*/
		inline      auto      getY_Axis() const->Coords;

		/*
		@brief     Returns mZ_Axis class member object.
		*/
		inline      auto      getZ_Axis() const->Coords;

		/*
		@brief     Returns m_Type class member object.
		*/
		inline      auto      getType() const->std::string;

		/*
		@brief     Returns m_isFixed class member object.
		*/
		inline      auto      isFixed() const->bool;

		/*
		@brief     Returns m_OrientationQuaternions class member object.
		*/
		inline      auto      getOrientQuaternions() const->std::array<mathlib::Quaternion, 3U>;

		/*
		@brief     Returns  mX_Origin class member object
		*/
		inline     auto       getXOrigin()const-> double;

		/*
		@brief     Returns  mY_Origin class member object.
		*/
		inline     auto       getYOrigin()const-> double;

		/*
		@brief     Returns  mZ_Origin  class member object.
		*/
		inline     auto       getZOrigin()const-> double;


		/***************************************************************
		               Class member functions
		***************************************************************/

		/*
		@brief     Displays the context of the Cartesian3D object.
		*/
		inline      auto      displayCartesian3D() const->void;

		

		/*
		@brief     Computes Cartesian3D cross product from class members type: Coords.
		*/
		inline      auto      cross_prod()->Cartesian3D &;


	private:

		/*
		@brief    class member mX_Axis which represents X-axis of Cartesian system.
		*/
		Coords  mX_Axis;

		/*
		@brief    class member mY_Axis which represents Y-axis of Cartesian system.
		*/
		Coords  mY_Axis;
			
			

		/*
		@brief    class member mZ_Axis which represents Z-Axis of Cartesian system.
		*/
		Coords  mZ_Axis;

		/*
		@brief    name of coordinate system i.e fixed earth frame, wind frame, missile local frame
		*/
		std::string  m_Type;

		/*
		@brief    boolean class member which denotes if specific coordinate system is fixed.
		*/
		bool       m_isFixed;

		/*
		@brief     Cartesian3D Origin x-coordinate. 
		
		*/
		double     mX_Origin;

		/*
		@brief     Cartesian3D Origin y-coordinate.
		*/
		double     mY_Origin;

		/*
		@brief     Cartesian3D Origin z-coordinate.
		*/
		double     mZ_Origin;
		/*
		@brief    class member m_OrientationQuaternions.		
		*/
		std::array<gms::math::Quaternion, 3U>  m_OrientationQuaternions;
	};


	/*
	
	        ***    Free standing functions operating on CartesianSystem class   ***
			TO DO: Implement special class for inter-coordinate system distance calculation.
	*/

	

	/*
	@brief       Cartesian3D magnitude.
	*/
	//template<class Coords, typename T> auto  Cartesian3DMag(_In_ const mathlib::Cartesian3D<Coords> &)->T;

	/*
	@brief       Cartesian3D perpendicular.
	*/
	//template<class Coords, typename T> auto Cartesian3DPer(_In_ const mathlib::Cartesian3D<Coords> &)->T;

	/*
	@brief       Cartesian3D rho value.
	*/
	//template<class Coords, typename T> auto Cartesian3DRho(_In_ const mathlib::Cartesian3D<Coords> &)->T;

	/*
	@brief       Cartesian3D theta value.
	*/
	//template<class Coords, typename T> auto Cartesian3DTheta(_In_ const mathlib::Cartesian3D<Coords> &)->T;

	/*
	@brief        Cartesian3D phi value.
	*/
	//template<class Coords, typename T> auto Cartesian3DPhi(_In_ const mathlib::Cartesian3D<Coords> &)->T;



#include "GMS_cartesian_system.inl"
   }
   
}
#endif /*__GMS_CARTESIAN_SYSTEM_H__*/
