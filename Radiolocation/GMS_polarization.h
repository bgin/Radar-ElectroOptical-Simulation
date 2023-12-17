#ifndef __GMS_POLARIZATION_H__
#define __GMS_POLARIZATION_H__

/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
 file Polarization.h
 class JonesVector
 class JonesMatrix
@aulthor: Bernard Gingold
@version:  1.0  19/10/2015

*/
//#include "RadLibDefs.h"
#include <complex>

namespace gms {
       namespace radiolocation {



	class JonesVector
	{

	public:

		/*
		@Brief: Explicitly default Ctor - null vector.
		*/
		 JonesVector() = default;

		 JonesVector(  std::complex<double> const,   std::complex<double> const);

		 JonesVector(  struct JonesVectorParams const&);

		 JonesVector( JonesVector const&);

		JonesVector( JonesVector &&);

		~JonesVector() = default;

		JonesVector&     operator=( JonesVector const&);

		JonesVector&     operator=( JonesVector &&);

		friend    std::ostream&    operator<<( std::ostream&,  JonesVector const&);

		JonesVector&     operator*=( JonesVector const&);

		JonesVector&     operator*=( std::complex<double> const&);

		JonesVector&     operator/=( JonesVector const&);

		JonesVector&     operator/=( std::complex<double> const&);

		JonesVector&     operator-=( JonesVector const&);

		JonesVector&     operator-=( std::complex<double> const&);

		JonesVector&     operator+=( JonesVector const&);

		JonesVector&     operator+=( std::complex<double> const&);

		JonesVector      operator-() const;

		

		std::complex<double>  operator*( JonesVector const&);

		friend	JonesVector   operator*( JonesVector const&,  JonesVector const&);

		friend  JonesVector   operator*( JonesVector const&,  std::complex<double> const&);

		friend  JonesVector   operator/( JonesVector const&,  JonesVector const&);

		friend  JonesVector   operator/( JonesVector const&,  std::complex<double> const&);

		friend  JonesVector   operator-( JonesVector const&,  JonesVector const&);

		friend  JonesVector   operator-( JonesVector const&,  std::complex<double> const&);

		friend  JonesVector   operator+( JonesVector const&,  JonesVector const&);

		friend  JonesVector   operator+( JonesVector const&,  std::complex<double> const&);
		

		__forceinline 	std::complex<double>   h();

		__forceinline   std::complex<double>   v();

		__forceinline  const  std::complex<double>  h() const;

		__forceinline  const  std::complex<double>  v() const;

		double     field_intensity() const;

		double     field_atan() const;

		double     field_phase_diff() const;

		double     degree_polarization() const;

	private:

		/*
		@Brief: horizontal component of electrical field
		*/
		std::complex<double> m_h;

		/*
		@Brief: vertical component of electrical field.
		*/
		std::complex<double> m_v;

	};

	class JonesMatrix
	{


	public:

		JonesMatrix() = default;

		JonesMatrix( std::complex<double> const,  std::complex<double> const,  std::complex<double> const,
			 std::complex<double> const);

		JonesMatrix( struct JonesMatrixParams const&);

		JonesMatrix( JonesMatrix const&);

		JonesMatrix( JonesMatrix &&);

		~JonesMatrix() = default;

		inline   std::complex<double>    s0();

		inline   const   std::complex<double>  s0() const;

		inline   std::complex<double>    s1();

		inline   const   std::complex<double>  s1() const;

		inline   std::complex<double>    s2();

		inline   const  std::complex<double>  s2() const;

	        inline   std::complex<double>    s3();

		inline   const  std::complex<double>  s3() const;

		inline   std::complex<double>  *matrix();

		inline   const  std::complex<double>  *matrix() const;

		   JonesMatrix&    operator=( JonesMatrix const&);

		   

		   JonesMatrix&    operator=( JonesMatrix &&);

		   JonesMatrix&    operator+=(  JonesMatrix const&);

		   JonesMatrix&    operator+=( std::complex<double> const&);

		   JonesMatrix&    operator+=(  JonesMatrix &&);

		   JonesMatrix&    operator-=(  JonesMatrix const&);

		   JonesMatrix&    operator-=(  std::complex<double> const&);

		   JonesMatrix&    operator-=(  JonesMatrix &&);

		   JonesMatrix&    operator*=(  JonesMatrix const&);

		   JonesMatrix&    operator*=(  std::complex<double> const&);

		   JonesMatrix&    operator*=(  JonesMatrix &&);

		   JonesMatrix&    operator/=(  JonesMatrix const&);

		   JonesMatrix&    operator/=(  std::complex<double> const&);

		   JonesMatrix&    operator/=(  JonesMatrix &&);

		   std::complex<double>     operator[]( const int);
		   
		   const std::complex<double>   operator[]( const int) const;

		   JonesMatrix     operator-() const;

		   

		   friend  JonesMatrix   operator+( JonesMatrix const&,  JonesMatrix const&);

		   friend  JonesMatrix   operator+( JonesMatrix const& , std::complex<double> const&);

		   friend  JonesMatrix   operator-( JonesMatrix const&,  JonesMatrix const&);

		   friend JonesMatrix   operator-( JonesMatrix const&,  std::complex<double> const&);

		   friend  JonesMatrix   operator*( JonesMatrix const&,  JonesMatrix const&);

		   friend  JonesMatrix   operator*( JonesMatrix const&,  std::complex<double> const&);

		   friend  JonesMatrix   operator/( JonesMatrix const&,  JonesMatrix const&);

		  

		   friend  JonesMatrix   operator/( JonesMatrix const&,  std::complex<double> const&);

		   friend  std::ostream&     operator<<( std::ostream&,  JonesMatrix const&);

		   JonesMatrix   matrix_transpose() const;

		   JonesMatrix   matrix_hermitian() const;
	private:

		std::complex<double> m_matrix[4];
	};


	class StokesVector
	{

	public:

		
		StokesVector() = default;

		StokesVector( const double,  const double,  const double,  const double);

		StokesVector( struct StokesVectorParams const&);

		StokesVector( __m256d const&);

		StokesVector( StokesVector const&);

		StokesVector( StokesVector &&);

		~StokesVector() = default;

		/*
		@Brief 1st element of Stokes vector.
		*/
	        inline 	double            s0() const;
	/*
	@Brief: 2nd element of Stokes vector.
	*/
	        inline   double            s1() const;
	/*
	@Brief: 3rd element of Stokes vector.
	*/
	        inline   double            s2() const;
	/*
	@Brief: 4th element of Stokes vector.
	*/
	        inline   double            s3() const;
	/*
	@Brief: Constant pointer to this->m_stokes.
	*/
	        inline  const  double           *stokes() const;

	/*
	@Brief: Normalized Stokes vector elements.
	*/

	       inline  double  normalized_s1() const;

	       inline  double  normalized_s2() const;

	       inline  double  normalized_s3() const;

	/*
	@Brief: Stokes vector dot product.
	*/
	       inline  double            dot( StokesVector const&) const;

	/*
	@Brief: Principal angle of polarization
	*/
	double        poa() const;

	/*
	@Brief: Circular polarization degree.
	*/
	double        dcp() const;

	/*
	@Brief: Polarization ellipticity.
	*/
	double        polarization_ellipticity() const;


	/*
	@Brief: Elliptical plarization eccentrity
	*/
	double        polarization_eccentrity() const;

	/*
	@Brief: Zero Stokes vector.
	*/
	static   StokesVector    zero_stokes_vector();

	/*
	@Brief: Non-polarized Stokes vector.
	*/
	static   StokesVector    nonpolarized_stokes_vector();

		StokesVector&     operator=( StokesVector const&);

		StokesVector&     operator=( StokesVector &&);

		friend	std::ostream&     operator<<( std::ostream&,  StokesVector const&);

		double            operator[]( const int);

		const  double     operator[]( const int) const;

		StokesVector&     operator+=( StokesVector const&);

		StokesVector&     operator-=( StokesVector const&);

		StokesVector&     operator*=( StokesVector const&);

		StokesVector&     operator/=( StokesVector const&);


		StokesVector&     operator+=( const double);

		StokesVector&     operator-=( const double);

		StokesVector&     operator*=( const double);

		StokesVector&     operator/=( const double);

		

		StokesVector      operator-() const;

		friend   StokesVector   operator+( StokesVector const&,  StokesVector const&);

		

		friend   StokesVector   operator-( StokesVector const&,  StokesVector const&);

		friend   StokesVector   operator*( StokesVector const&,  StokesVector const&);

		friend   StokesVector   operator/( StokesVector const&,  StokesVector const&);

		friend   StokesVector   operator+( StokesVector const&,  const double);

		friend   StokesVector   operator-( StokesVector const&,  const double);

		friend   StokesVector   operator*( StokesVector const&,  const double);

		friend   StokesVector   operator/( StokesVector const&,  const double);

	private:

		double m_stokes[4];
	};

	

	/*
	"Brief:  Aggregate structure of JonesMatrix parameters.
	*/
	
	struct JonesMatrixParams
	{
		std::complex<double> ps0;
		std::complex<double> ps1;
		std::complex<double> ps2;
		std::complex<double> ps3;
	};

	/*
	@Brief: Aggregated structure of JonesVector parametrs.
	*/

	struct JonesVectorParams
	{
		std::complex<double> ps0;
		std::complex<double> ps1;
	};

	/*
	@Brief: Aggregated structure of StokesVectore parameters.
	*/

	struct StokesVectorParams
	{
		double ps0;
		double ps1;
		double ps2;
		double ps3;
	};

#include "GMS_polarization.inl"
   }
}
#endif /*__GMS_POLARIZATION_H__*/
