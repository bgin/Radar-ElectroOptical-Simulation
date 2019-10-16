
#ifndef __GMS_QUATERNION_H__
#define __GMS_QUATERNION_H__



namespace file_info {
#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif

    const unsigned int gGMS_QUATERNION_MAJOR = gms::common::gVersionInfo.m_VersionMajor;
    const unsigned int gGMS_QUATERNION_MINOR = gms::common::gVersionInfo.m_VersionMinor;
    const unsigned int gGMS_QUATERNION_MICRO = gms::common::gVersionInfo.m_VersionMicro;
    const unsigned int gGMS_QUATERNION_FULLVER =
      1000U*gGMS_QUATERNION_MAJOR+100U*gGMS_QUATERNION_MINOR+10U*gGMS_QUATERNION_MICRO;
    const char * const pgGMS_QUATERNION_CREATE_DATE = "06-12-2017 08:13 +00200 (WED 06 DEC 2017 GMT+2)";
    const char * const pgGMS_QUATERNION_BUILD_DATE  = "00-00-0000 00:00";
    const char * const pgGMS_QUATERNION_AUTHOR      = "Programmer: Bernard Gingold contact: beniekg@gmail.com";
    const char * const pgGMS_QUATERNION_SYNOPSIS    = "Implementation of Quaternion class.";
}

#include <iosfwd>
#include <iomanip>
#include <cstdint>
#include <complex>
  
namespace gms {
	namespace math {

	using C64 = std::complex<double>;
#if !defined (CHECK_QUATERNION_LEPS)
#define CHECK_QUATERNION_LEPS(value) \
	do {							 \
	     if ((value) < std::numeric_limits<double>::epsilon()) \
		  (value) += std::numeric_limits<double>::epsilon();  \
	} while(0); 
#endif

#if !defined (GMS_QUATERNION_ACCESSORS)
#define GMS_QUATERNION_ACCESSORS    \
									\
	inline	double get_x() const {	\
		return (this->m_x);			\
	}								\
									\
	inline double get_y() const {	\
		return (this->m_y);			\
	}								\
									\
	inline double get_z() const {	\
		return (this->m_z);			\
	}								\
									\
	inline double get_w() const {	\
		return (this->m_w);			\
	}								\
									\
	inline double real()  const {	\
		return (this->get_x());		\
	}								\
									\
	inline Quaternion unreal() const { \
		return (Quaternion{ 0.0,		\
		this->m_y, this->m_z, this->m_w }); \
	}								\
									\
	inline C64 get_c1() const {		\
		return (C64{ m_x, m_y });	\
	}								\
									\
	inline C64 get_c2() const {		\
		return (C64{ m_z, m_w });	\
	}



#endif

#if !defined (GMS_QUATERNION_CTORS_DECLARATIONS)
#define GMS_QUATERNION_CTORS_DECLARATIONS    \
											 \
	Quaternion();							 \
											 \
											 \
	explicit Quaternion(_In_ const double);  \
											 \
											 \
	Quaternion(const double,			 \
		   const double,					 \
		   const double);					 \
											 \
	Quaternion(const double,			 \
		const double,					 \
		const double,					 \
		const double);					 \
											 \
											 \
	Quaternion(const C64 &,			 \
		   const C64 &);			 \
											 \
										     \
	Quaternion(const Quaternion &);     \
											
#endif

		class Quaternion {

			
			public:

			GMS_QUATERNION_CTORS_DECLARATIONS
			//
			// Destructor - default
			//
			~Quaternion() = default;

			//
			// Member operators
			//
			Quaternion & operator=(const Quaternion &);

			Quaternion & operator=(const double );

			Quaternion & operator=(const C64 &);

			friend std::ostream & operator<<(std::ostream &,
							 const Quaternion &);

			const double operator[](const uint32_t) const;

			//
			// Getters
			//
			GMS_QUATERNION_ACCESSORS

			private:

			// Quaternion 4-tuple components.

			double m_x;

			double m_y;

			double m_z;

			double m_w;
		};

#if !defined (Q_x)
#define Q_x(obj) obj.get_x()
#endif

#if !defined (Q_y)
#define Q_y(obj) obj.get_y()
#endif

#if !defined (Q_z)
#define Q_z(obj) obj.get_z()
#endif

#if !defined (Q_w)
#define Q_w(obj) obj.get_w()
#endif

#if !defined (Q_SQR)
#define Q_SQR(x) ((x) * (x))
#endif



		//
		//	Global operators
		//

		static inline Quaternion operator+(const Quaternion &,
						   const Quaternion &);

		static inline Quaternion operator+(const Quaternion &,
						   const std::complex<double> &);

		static inline Quaternion operator+(const Quaternion &,
						   const double);

		static inline Quaternion operator+(const std::complex<double> &,
						   const Quaternion &);

		static inline Quaternion operator+(const double,
						   const Quaternion &);

		static inline Quaternion operator+=(Quaternion &,
						    const Quaternion &);

		static inline Quaternion operator+=(Quaternion &,
						    const std::complex<double> &);

		static inline Quaternion operator+=(const C64 &,
						    Quaternion &);

		static inline Quaternion operator+=(Quaternion &,
						    const double);

		static inline Quaternion operator+=(const double,
						    Quaternion &);

		static inline Quaternion operator-(const Quaternion &,
						   const Quaternion &);

		static inline Quaternion operator-(const Quaternion &,
						   const std::complex<double> &);

		static inline Quaternion operator-(const Quaternion &,
						   const double);

		static inline Quaternion operator-(const std::complex<double> &,
						   const Quaternion &);

		static inline Quaternion operator-(const double,
						   const Quaternion &);

		static inline Quaternion operator-=(Quaternion &,
						   const Quaternion &);

		static inline Quaternion operator-=(Quaternion &,
						     const std::complex<double> &);

		static inline Quaternion operator-=(const C64 &,
						     Quaternion &);

		static inline Quaternion operator-=(Quaternion &,
						    const double);

		static inline Quaternion operator-=( const double,
						      Quaternion &);

		static inline Quaternion operator*( const Quaternion &,
						     const Quaternion &);

		static inline Quaternion operator*( const Quaternion &,
						         const std::complex<double> &);

		static inline Quaternion operator*( const Quaternion &,
						      const double);
		
		static inline Quaternion operator*(const std::complex<double> &,
						    const Quaternion &);

		static inline Quaternion operator*(const double,
						   const Quaternion &);

		static inline Quaternion operator*=(Quaternion &,
						     const Quaternion &);

		static inline Quaternion operator*=(Quaternion &,
						    const std::complex<double> &);

		static inline Quaternion operator*=(const C64 &,
						    Quaternion &);

		static inline Quaternion operator*=(Quaternion &,
						     const double );

		static inline Quaternion operator*=(const double,
					             Quaternion &);

		static inline Quaternion operator/( const Quaternion &,
						    const Quaternion &);

		static inline Quaternion operator/( const Quaternion &,
						    const std::complex<double> &);

		static inline Quaternion operator/(const Quaternion &,
						    const double);

		static inline Quaternion operator/(const std::complex<double> &,
						    const Quaternion &);

		static inline Quaternion operator/(const double,
						   const Quaternion &);

		static inline Quaternion operator/=(Quaternion &,
						    const Quaternion &);

		static inline Quaternion operator/=(Quaternion &,
						    const std::complex<double> &);

		static inline Quaternion operator/=(const C64 &,
						    Quaternion &);

		static inline Quaternion operator/=(Quaternion &,
					            const double);

		static inline Quaternion operator/=(const double,
						     Quaternion &);

		static inline bool       operator==(const Quaternion &,
						     const Quaternion &);

		static inline bool       operator==(const Quaternion &,
						     const std::complex<double> &);

		static inline bool       operator==( const Quaternion &,
						     const double);

		static inline bool       operator==(const std::complex<double> &,
						    const Quaternion &);

		static inline bool       operator==(const double,
						    const Quaternion &);

		static inline bool       operator!=(const Quaternion &,
						    const Quaternion &);

		static inline bool       operator!=(const Quaternion &,
						     const std::complex<double> &);

		static inline bool       operator!=( const Quaternion &,
						      const double);

		static inline bool       operator!=( const std::complex<double> &,
						     const Quaternion &);

		static inline bool       operator!=( const double,
						     const Quaternion &);

		static inline bool       operator>(const Quaternion &,
					          const Quaternion &);

		static inline bool       operator>=(const Quaternion &,
						     const Quaternion &);

		static inline bool       operator<(const Quaternion &,
					           const Quaternion &);

		static inline bool       operator<=(const Quaternion &,
						    const Quaternion &);

		//
		// Free-standing functions operating on
		// various Quaternion-related characteristics like:
		// Norm, Vector part norm, conjugate,distance ...etc
		//

		static inline Quaternion conjugate(const Quaternion &);

		static inline double norm(const Quaternion &);

		static inline double vnorm(const Quaternion &);

		static inline double distance( const Quaternion &,
						 const Quaternion &);

		static inline Quaternion unit(const Quaternion &);

		static inline Quaternion polar_decomp( const Quaternion &);

		static inline Quaternion reciprocal(const Quaternion &);

		static inline void mat4x4(const Quaternion &,
					    double (&m4x4)[4][4],
					    const int32_t );

		static inline Quaternion exp(const Quaternion &);

		static inline Quaternion ln(const Quaternion &);

		//
		// Based on  boost::quaternion
		//
		static inline Quaternion spherical(const double,
						   const double,
						   const double,
						   const double);

		static inline Quaternion semipolar(const double,
						   const double,
						    const double,
						    const double);

		static inline Quaternion multipolar(const double,
						    const double,
						     const double,
						     const double);

		static inline Quaternion cylindrospherical(const double,
							   const double,
							   const double,
							   const double);

		static inline Quaternion cylindrical(const double,
						     const double,
						     const double,
						     const double);

		

#include "GMS_quaternion.inl"
	}
}

#endif /*__GMS_QUATERNION_H__*/
