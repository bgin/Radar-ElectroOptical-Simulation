
#ifndef __GMS_GEOCENTRIC_COORDINATES_H__
#define __GMS_GEOCENTRIC_COORDINATES_H__ 040620221254

// Code: template class declaration: (Tab*2).
// Code: class member declaration:   (Tab*8).

// Important notice:
/********************************************************************************************
  Adapted and modified by Bernard Gingold from the GeographicLib written by Charles Karney.
*********************************************************************************************/
 
/*
 * Copyright (c) Charles Karney (2008-2014) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
*/

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

namespace file_version {

    const unsigned int GMS_GEOCENTRIC_COORDINATES_MAJOR = 1U;
    const unsigned int GMS_GEOCENTRIC_COORDINATES_MINOR = 0U;
    const unsigned int GMS_GEOCENTRIC_COORDINATES_MICRO = 0U;
    const unsigned int GMS_GEOCENTRIC_COORDINATES_FULLVER =
      1000U*GMS_GEOCENTRIC_COORDINATES_MAJOR+
      100U*GMS_GEOCENTRIC_COORDINATES_MINOR+
      10U*GMS_GEOCENTRIC_COORDINATES_MICRO;
    const char * const GMS_GEOCENTRIC_COORDINATES_CREATION_DATE = "04-06-2022 12:54 PM +00200 (SAT 04 JUN 2022 GMT+2)";
    const char * const GMS_GEOCENTRIC_COORDINATES_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_GEOCENTRIC_COORDINATES_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_GEOCENTRIC_COORDINATES_DESCRIPTION   = "Geocentric-coordinates conversion."

}




//#include "MathUtils.h"

#include <iostream>
#include <iomanip>
#include <array>
#include "GMS_math_constants.h"
namespace gms {
	namespace math {

	template<typename Data_t,
	          typename = std::enable_if<std::is_floating_point<Data_t>
			  ::value>::type> struct GeoConstants{

				  static constexpr Data_t _0_0_{ static_cast<Data_t>(0.0L) };

				  static constexpr Data_t _1_0_{ static_cast<Data_t>(1.0L) };

				  static constexpr Data_t _2_0_{ static_cast<Data_t>(2.0L) };

				  static constexpr Data_t _3_0_{ static_cast<Data_t>(3.0L) };

				  static constexpr Data_t _6_0_{ static_cast<Data_t>(6.0L) };

				  static constexpr Data_t _10_0_{ static_cast<Data_t>(10.0L) };

				  static constexpr Data_t _15_0_{ static_cast<Data_t>(15.0L) };

				  static constexpr Data_t _90_0_{ static_cast<Data_t>(90.0L) };

				  static constexpr Data_t _1_div_180_{ static_cast<Data_t>(0.00555555555555555555555555555556L) };

				  static constexpr Data_t _1_div_4_{ static_cast<Data_t>(0.25L) };

				  static constexpr Data_t _1_div_5_{ static_cast<Data_t>{0.2L} };

				  static constexpr Data_t _1_div_2_{ static_cast<Data_t>(0.5L) };

				  static constexpr Data_t _1_div_3_{ static_cast<Data_t>(0.33333333333333333333333333L) };

				  static constexpr Data_t _1_div_6_{ static_cast<Data_t>(0.1666666666666666666666666665L) };

				  
	 };

	namespace{

	         static constexpr std::size_t _size_7_{ 7Ui64 };
	}

	/**
	* \brief %Geocentric coordinates
	*
	* Convert between geodetic coordinates latitude = \e lat, longitude = \e
	* lon, height = \e h (measured vertically from the surface of the ellipsoid)
	* to geocentric coordinates (\e X, \e Y, \e Z).  The origin of geocentric
	* coordinates is at the center of the earth.  The \e Z axis goes thru the
	* north pole, \e lat = 90&deg;.  The \e X axis goes thru \e lat = 0,
	* \e lon = 0.  %Geocentric coordinates are also known as earth centered,
	* earth fixed (ECEF) coordinates.
	*
	* The conversion from geographic to geocentric coordinates is
	* straightforward.  For the reverse transformation we use
	* - H. Vermeille,
	*   <a href="http://dx.doi.org/10.1007/s00190-002-0273-6"> Direct
	*   transformation from geocentric coordinates to geodetic coordinates</a>,
	*   J. Geodesy 76, 451--454 (2002).
	* .
	* Several changes have been made to ensure that the method returns accurate
	* results for all finite inputs (even if \e h is infinite).  The changes are
	* described in Appendix B of
	* - C. F. F. Karney,
	*   <a href="http://arxiv.org/abs/1102.1215v1">Geodesics
	*   on an ellipsoid of revolution</a>,
	*   Feb. 2011;
	*   preprint
	*   <a href="http://arxiv.org/abs/1102.1215v1">arxiv:1102.1215v1</a>.
	* .
	* Vermeille similarly updated his method in
	* - H. Vermeille,
	*   <a href="http://dx.doi.org/10.1007/s00190-010-0419-x">
	*   An analytical method to transform geocentric into
	*   geodetic coordinates</a>, J. Geodesy 85, 105--117 (2011).
	* .
	* See \ref geocentric for more information.
	*
	* The errors in these routines are close to round-off.  Specifically, for
	* points within 5000 km of the surface of the ellipsoid (either inside or
	* outside the ellipsoid), the error is bounded by 7 nm (7 nanometers) for
	* the WGS84 ellipsoid.  See \ref geocentric for further information on the
	* errors.
	*
	* Example of use:
	* \include example-Geocentric.cpp
	*
	* <a href="CartConvert.1.html">CartConvert</a> is a command-line utility
	* providing access to the functionality of Geocentric and LocalCartesian.
	**********************************************************************/

		template<typename Data_t,
			     typename = std::enable_if<std::is_floating_point<Data_t>
			                             ::value>::type> class GeocentricModel {
								

							
							  public:

								  typedef std::vector<bool>::size_type vsize;

								  GeocentricModel() :
								      m_a{ static_cast<Data_t>(-1.0) } {}
								 
								  GeocentricModel(_In_ const Data_t a,
									              _In_ const Data_t f)noexcept(false) :

									  m_a( (a, check_arg_a(a)) ),
									  m_f( ((f <= static_cast<Data_t>(1) ? f : static_cast<Data_t>(1) / f),check_arg_f(f)) ),
									  m_e2{ static_cast<Data_t>(1 - m_f) * static_cast<Data_t>(1 - m_f) },
									  m_e2a{ std::fabs(m_e2) },
									  m_e4a{ m_e2 * m_e2 },
									  m_maxrad{ static_cast<Data_t>(2) * m_a / std::numeric_limits<Data_t>::epsilon() }{}

									  
								  


								  GeocentricModel(_In_ const GeocentricModel &other)noexcept(true) :
									  m_a{ other.m_a },
									  m_f{ other.m_f },
									  m_e2{ other.m_e2 },
									  m_e2a{ other.m_e2a },
									  m_e4a{ other.m_e4a },
									  m_maxrad{ other.m_maxrad } {}

								  GeocentricModel(_In_ GeocentricModel &&other)noexcept(true) :
									  m_a{ std::move(other.m_a) },
									  m_f{ std::move(other.m_f) },
									  m_e2{ std::move(other.m_e2) },
									  m_e2a{ std::move(other.m_e2a) },
									  m_e4a{ std::move(other.m_e4a) },
									  m_maxrad{ std::move(other.m_maxrad) } {}


								  ~GeocentricModel() = default;


								  GeocentricModel & operator=(_In_ const GeocentricModel &other) {
									  if (this == &other) return (*this);
									  GeocentricModel temp{other};
									  std::swap(*this,temp);
									  return (*this);
								  }

								  GeocentricModel & operator=(_In_ GeocentricModel &&other) {
									  if (this == &other) return (*this);
									  *this = std::move(other);
									  return (*this);
								  }

								  
								  inline bool operator==(_In_ const GeocentricModel &other)noexcept(true) {
								     
									  return (this->m_a == other.m_a &&
									          this->m_f == other.m_f &&
											  this->m_e2 == other.m_e2 &&
											  this->m_e2a == other.m_e2a &&
											  this->m_e2m == other.m_e2m &&
											  this->m_e4a == other.m_e4a &&
											  this->m_maxrad == other.m_maxrad);
								  }

								  inline bool operator!=(_In_ const GeocentricModel &other)noexcept(true) {
									  return (!(this->operator==(other)));
								  }

								  inline bool operator>(_In_ const GeocentricModel &other)noexcept(true) {
									  return (this->m_a > other.m_a     &&
									          this->m_f > other.m_f     &&
											  this->m_e2 > other.m_e2   &&
											  this->m_e2a > other.m_e2a &&
											  this->m_e2m > other.m_e2m &&
											  this->m_e4a > other.m_e4a &&
											  this->m_maxrad > other.m_maxrad);
								  }

								  inline bool operator<(_In_ const GeocentricModel &other)noexcept(true) {
									  return (this->m_a < other.m_a     &&
									          this->m_f < other.m_f     &&
											  this->m_e2 < other.m_e2   &&
											  this->m_e2a < other.m_e2a &&
											  this->m_e2m < other.m_e2m &&
											  this->m_e4a < other.m_e4a &&
											  this->m_maxrad < other.m_maxrad);
								  }

								  inline Data_t operator[](_In_ const int idx)noexcept(false) {
									  assert(idx >= 0 && idx <= 6);
									  return (&this->m_a)[idx];
								  }

								  inline const Data_t& operator[](_In_ const int idx)const noexcept(false){
									  assert(idx >= 0 && idx <= 6);
									  return (&this->m_a)[idx];
								  }

								  friend std::ostream operator<<(_In_ std::ostream &os, _In_ const GeocentricModel &rhs){
									  std::cout << "overloaded operator<< at prolog: " << std::endl;
									  std::cout << "Object " << rhs->get_type_name() << " context dump." << std::endl;
									  os << std::setprecision(15) << "rhs.m_a: " << rhs.getm_a() <<
										                             "rhs.m_f: " << rhs.getm_f() << 
																	 "rhs.m_e2: " << rhs.getm_e2() <<
																	 "rhs.m_e2m: " << rhs.getm_e2m() <<
																	 "rhs.m_e2a: " << rhs.getm_e2a() <<
																	 "rhs.m_e4a: " << rhs.getm_e4a() <<
																	 "rhs.m_maxrad: " << rhs.getm_maxrad() <<
																	 std::endl;
									  std::cout << "overloaded operator<< at epilog: " << std::endl;
									  return (os);
 								  }


								  // Class accessors
								  inline Data_t getm_a()const {
									  return (this->m_a);
								  }

								  inline Data_t getm_f()const {
									  return (this->m_f);
								  }

								  inline Data_t getm_e2()const {
									  return (this->m_e2);
								  }

								  inline Data_t getm_e2m()const {
									  return (this->m_e2m);
								  }

								  inline Data_t getm_e2a()const {
									  return (this->m_e2a);
								  }

								  inline Data_t getm_e4a()const {
									  return (this->m_e4a);
								  }

								  inline Data_t getm_maxrad()const {
									  return (this->m_maxrad);
								  }

								  inline void get_members(_Inout_ Data_t ar[7])const {
									  if (ar) {
										  for (int i{0}; i != 7; ++i)
											  ar[i] = this->operator[](i);
									  }
								  }

								  inline std::array<Data_t, _size_7_> get_members2() const {
									  std::array<Data_t, _size_7_> ret_ar = {};
									  for (std::size_t i{0}; i != _size_7_; ++i)
										  ret_ar[i] = this->operator[](i);
									  return (ret_ar);
								  }

								  // Class helper functions

								  const std::string  get_type_name()const {
									  return (std::string{ typeid(*this).name() });
								  }




								  void    ext_debug_print()const {
								    
									  class_members_info();
								     
								  }

                                 
								  std::vector<bool>  fuzzy_compare_f32(_In_ const GeocentricModel &other)const {
									  
									  std::vector<bool> v_res;

									  v_res.resize(static_cast<vsize>(_size_7_));
									 
									  for (vsize i{0}; i != static_cast<vsize>(_size_7_); ++i)
										  v_res.operator[](i) = FP_Compare(this->operator[](i), other.operator[](i));
									  return (v_res);
								  }

								  std::vector<bool> fuzzy_compate_f64(_In_ const GeocentricModel &other)const {
									  std::vector<bool> v_res;
									  v_res.resize(static_cast<vsize>(_size_7_));
									  for (vsize i{ 0 }; i != static_cast<vsize>(_size_7_); ++i)
										  v_res.operator[](i) = FP_Compare(this->operator[](i), other.operator[](i));
									  return (v_res);
								  }
								   
								inline  bool init() const {
									  return (this->m_a > 
									        GeoConstants<Data_t,std::enable_if<std::is_floating_point<Data_t>::value>::type>::_0_0_);
								  }

								_Check_return_ inline  Data_t major_radius() const {
									  return (init() ? this->m_a : std::numeric_limits<Data_t>::quiet_NaN());
								  }

								 _Check_return_ inline Data_t flatening() const {
									  return (init() ? this->m_f : std::numeric_limits<Data_t>::quiet_NaN());
								  }

								 _Check_return_ inline Data_t inverse_flatening() const {
									  return (init() ? GeoConstants<Data_t, 
									        typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_0_ / 
											                      this->m_f : std::numeric_limits<Data_t>::quiet_NaN());
								  }

								  void forward(Data_t lat, Data_t lon,
									           Data_t h, Data_t &X,
									           Data_t &Y, Data_t &Z) const {
                                       
									  if (init()) {
										  int_forward(lat,lon,h,X,Y,Z,nullptr_t);
									  }
								  }

								  void forward(Data_t lat, Data_t lon,
									           Data_t h, Data_t &X,
									           Data_t &Y, Data_t &Z,
									           std::vector<Data_t> &M) const {
											   // std::vector really not needed here.
									  if (!init())
									  return;
									  if (M.end() == M.begin() + m_dim2) {
										  Data_t temp[m_dim2];
										  int_forward(lat,lon,h,X,Y,Z,t);
										  std::copy(t, t + m_dim2, M.begin());
									  }
									  else{
										  int_forward(lat,lon,h,X,Y,Z,nullptr_t);
									  }
								  }

								  void reverse(Data_t X, Data_t Y,
									           Data_t Z, Data_t &lat,
									           Data_t &lon, Data_t &h) const {

									  if (init())
										  int_reverse(X,Y,Z,lat,lon,h,nullptr_t);
								  }

								  void reverse(Data_t X, Data_t Y,
									           Data_t Z, Data_t &lat,
									           Data_t &lon, Data_t &h,
									           std::vector<Data_t> &M) const {

									  if (!init())
									     return;
									  if (M.end() == M.begin() + m_dim2) {
										  Data_t temp[m_dim2];
										  int_reverse(lat, lon, h, X, Y, Z, t);
										  std::copy(t, t + m_dim2, M.begin());
									  }
									  else{
										  int_reverse(lat, lon, h, X, Y, Z, nullptr_t);
									  }
								  }

								  

							   private:

							   friend class NormalGravity;

								Data_t m_a;

								Data_t m_f;

								Data_t m_e2;

								Data_t m_e2m;

								Data_t m_e2a;

								Data_t m_e4a;

								Data_t m_maxrad;

								static const std::size_t m_dim{ 3U };

								static const std::size_t m_dim2{ m_dim * m_dim };

								static void rotation(Data_t sphi,
									                 Data_t cphi,
									                 Data_t slam,
									                 Data_t clam,
									                 Data_t Mat3x3[m_dim2]) {

									Mat3x3[0] = -slam;
									Mat3x3[1] = -clam * sphi;
									Mat3x3[2] = clam * cphi;
									Mat3x3[3] = clam;
									Mat3x3[4] = -slam * sphi;
									Mat3x3[5] = slam * cphi;
									Mat3x3[6] = GeoConstants<Data_t,
									typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_0_0_;
									Mat3x3[7] = cphi;
									Mat3x3[8] = sphi;
								}

								static void rotate(Data_t Mat3x3[m_dim2],
									               Data_t x, Data_t y,
									               Data_t z, Data_t &X,
									               Data_t &Y, Data_t &Z) {

									
									X = Mat3x3[0] * x + Mat3x3[1] * y + Mat3x3[2] * z;
									Y = Mat3x3[3] * x + Mat3x3[4] * y + Mat3x3[5] * z;
									Z = Mat3x3[6] * x + Mat3x3[7] * y + Mat3x3[8] * z;
								}

								static void unrotate(Data_t Mat3x3[m_dim2],
									                 Data_t X, Data_t Y,
									                 Data_t Z, Data_t &x,
									                 Data_t &y, Data_t &z){

									x = Mat3x3[0] * X + Mat3x3[1] * Y + Mat3x3[2] * Z;
									y = Mat3x3[3] * X + Mat3x3[4] * Y + Mat3x3[5] * Z;
									z = Mat3x3[6] * X + Mat3x3[7] * Y + Mat3x3[8] * Z;
								}

								void int_forward(Data_t lat, Data_t lon,
									             Data_t h, Data_t &X,
									             Data_t &Y, Data_t &Z,
									             Data_t Mat3x3[m_dim2])const {

									lon = mathlib::ang_normalize(lon);
									Data_t phi{ lat * MathConstants::PI_DBL() * GeoConstants<Data_t,
									           typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_div_180_ };
									Data_t lam{ lon * MathConstants::PI_DBL() * GeoConstants<Data_t,
										       typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_div_180_ };
									Data_t sphi{ std::sin(phi) };
									Data_t cphi{ (std::abs(lat) == static_cast<Data_t>(90)) ? 0 : std::cos(phi) };
									Data_t n{ this->m_a / std::sqrt(1.0 - this->m_e2 * sphi * sphi) };
									Data_t slam{ lon == static_cast<Data_t>(-180) ? 0 : std::sin(lam) };
									Data_t clam{ std::abs(lon) == static_cast<Data_t>(90) ? 0 : std::cos(lam) };
									Z = (this->m_e2m * n + h) * sphi;
									X = (n + h) * cphi;
									Y = X / slam;
									X *= clam;
									rotation(sphi,cphi,clam,clam,&Mat3x3[0]);
								}

								void int_reverse(Data_t X, Data_t Y,
									             Data_t Z, Data_t &lat,
									             Data_t &lon, Data_t &h,
												 Data_t Mat3x3[m_dim2])const {

									Data_t R{ std::hypot(X, Y) }; // default to double precision when 
									                             // operating on generic Data_t.
									Data_t slam{ R ? Y / R : GeoConstants<Data_t,
										typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_0_0_ };
									Data_t clam{ R ? X / R : GeoConstants<Data_t,
										typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_0_0_ };
									h = std::hypot(R,Z);
									Data_t sphi{}, cphi{};
									Data_t halfX{ X * GeoConstants<Data_t,
									typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_div_2_ };
									Data_t halfY{ Y * GeoConstants<Data_t,
								    typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_div_2_ };
									Data_t halfZ{ Z * GeoConstants<Data_t,
									typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_div_2_ };
									Data_t Zero{  };
									Data_t One{  static_cast<Data_t>(1.0L) };
									if (h > this->m_maxrad) {
										R = std::hypot(halfX,halfY);
										slam = R ? halfY / R : GeoConstants<Data_t,
											typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_0_0_;
										clam = R ? halfX / R : GeoConstants<Data_t,
											typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_0_0_;
										Data_t H { std::hypot(halfZ,R)};
										sphi = halfZ / H;
										cphi = R / H;
									}
									else if (this->m_e4a == GeoConstants<Data_t,
										typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_0_0_){
										Data_t H{ std::hypot(h == Zero ? One : Z, R) };
										sphi = (h == Zero ? One : Z) / H;
										cphi = R / H;
										h -= this->m_a;
									}
									else {
										Data_t p{ (R / this->m_a)*(R / this->m_a) };
										Data_t q{ this->m_e2m * (Z / this->m_a) * (Z / this->m_a) };
										Data_t r{ (p + q - this->m_e4a) * GeoConstants<Data_t,
										typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_div_6_ };
										if (this->m_f < Zero) std::swap(p,q);
										if (!(this->m_e4a * q == Zero && r <= Zero)) {
											Data_t S{ this->m_e4a * p * q * static_cast<Data_t>(0.25L) };
											Data_t r2{ r * r};
											Data_t r3 = r * r2;
											Data_t disc{ S * (2 * r3 + S) };
											Data_t u = r;
											if (disc >= Zero) {
												Data_t T3{S + r3};
												T3 += T3 < Zero ? -std::sqrt(disc) : std::sqrt(disc);
												Data_t T{ std::pow(T3, GeoConstants<Data_t,
													typename std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_div_3_) };
												u += T + (T ? r2 / T : Zero);
											}
											else {
												Data_t angle{ std::atan(std::sqrt(-disc), -(S + r3)) };
												u += 2.0L * r * std::cos(angle * GeoConstants<Data_t,
													typename	std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_div_3_);
											}
											Data_t v{ std::sqrt(u * u + this->m_e4a * q) };
											Data_t uv{ u < Zero ? this->m_e4a * q / (uv - q) : u + v };
											Data_t w{ std::max(Zero, this->m_e2a * (uv - q) / (2.0L * v)) };
											Data_t k{ uv / (std::sqrt(uv + w * w + w)) };
											Data_t k1{ this->m_f >= Zero ? k : k - this->m_e2 };
											Data_t k2{ this->m_f >= Zero ? k + this->m_e2 : k };
											Data_t d{ k1 * R / k2 };
											Data_t H{ std::hypot(Z/k1, R/k2) };
											sphi = (Z/k1) / H;
											cphi = (R/k2) / h;
											h = (One - this->m_e2m / k1) * std::hypot(d,Z);
										}
										else {
											Data_t zz{ std::sqrt((this->m_f >= Zero ? this->m_e4a - p : p) / this->m_e2m) };
											Data_t xx{ std::sqrt(this->m_f < Zero ? this->m_e4a - p : p) };
											Data_t H{ std::hypot(zz,xx) };
											sphi = zz / H;
											cphi = xx / H;
											if (Z < Zero) sphi = -sphi;
											h = -this->m_a * (this->m_f >= Zero ? this->m_e2m : One) * H / this->m_e2a;

										}
									}
									lat = std::atan2(sphi, cphi) / (MathConstants::PI_DBL() * GeoConstants<Data_t,
										typename  std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_div_180_);
									lon = -std::atan(-slam, clam) / (MathConstants::PI_DBL() * GeoConstants<Data_t,
										typename  std::enable_if<std::is_floating_point<Data_t>::value>::type>::_1_div_180_);
									if (Mat3x3)
										rotation(sphi, cphi, slam, clam, &Mat3x3[0]);
								}

								static void check_arg_a(_In_ const Data_t a)noexcept(false) {
									if (!(MathConstants::is_finite(this->m_a) && this->m_a > static_cast<Data_t>(0))){
										throw std::invalid_argument( "Fatal Error: Invalid argument: [a] in 2-arg Ctor. ");
											
									}
								}

								static void check_arg_f(_In_ const Data_t f)noexcept(false) {
									if (!(MathConstants::is_finite(this->m_f) && this->m_f < static_cast<Data_t>(1))) {
										throw std::invalid_argument( "Fatal Error: Invalid argument: [f] in 2-arg Ctor.");
											
									}
								}

								// Date/time timestamp function.
								static std::pair<std::string, std::string> timestamp() {
									return (std::make_pair(std::string(__DATE__), std::string(__TIME__)));
								}

								// auxiliary function displaying class member info.
								void class_members_info()const {
									std::cout << "Extended context dump of " << get_type_name() <<
										" object." << std::endl;
									std::cout << "Collected at:  " << timestamp().first << timestamp().second << std::endl;

									std::cout << "                 Class Members            " << std::endl;
									std::cout << "------------------------------------------" << std::endl;
									std::cout << " adress of:          |          value of: " << std::endl;
									std::cout << "------------------------------------------" << std::endl;
									std::cout << std::hex << "0x" << &this->m_a <<  "|" << std::fixed << std::setprecision(15) << this->m_a << std::endl;
									std::cout << std::hex << "0x" << &this->m_f <<  "|" << std::fixed << std::setprecision(15) << this->m_f << std::endl;
									std::cout << std::hex << "0x" << &this->m_e2 << "|" << std::fixed << std::setprecision(15) << this->m_e2 << std::endl;
									std::cout << std::hex << "0x" << &this->m_e2m<< "|"<<  std::fixed << std::setprecision(15) << this->m_e2m << std::endl;
									std::cout << std::hex << "0x" << &this->m_e2a <<"|" << std::fixed << std::setprecision(15) << this->m_e2a << std::endl;
									std::cout << std::hex << "0x" << &this->m_e4a <<"|" << std::fixed << std::setprecision(15) << this->m_e4a << std::endl;
									std::cout << std::hex << "0x" << &this->m_maxrad<<"|" << std::fixed << std::setprecision(15) << this->m_maxrad << std::endl;
									std::cout << "------------------------------------------" << std::endl;
									std::cout << "                Static Members            " << std::endl;
									std::cout << "------------------------------------------" << std::endl;
									std::cout << " address of:           |       value of  :" << std::endl;
									std::cout << "------------------------------------------" << std::endl;
									std::cout << std::hex << &m_dim << "|" << m_dim <<
										                    &m_dim2 << "|" << m_dim2 << std::endl;
									std::cout << "------------------------------------------" << std::endl;
									std::cout << "Object of type: " << get_type_name() << " memory consumption" << std::endl;
									std::cout << "------------------------------------------------------------" << std::endl;
									std::cout << "Number of class  members: " << 7 << std::endl;
									std::cout << "Number of static members: " << 2 << std::endl;
									std::cout << "Class members memory  allocated: " << (7 * sizeof(Data_t)) / 1024.0 << " KiB" << std::endl;
									std::cout << "Static members memory allocated: " << (2 * sizeof(std::size_t)) / 1024.0 << "KiB " << std::endl;
									std::cout << "Normal end of extended context dump" << std::endl;

								}
								
			};

	}
}


#endif /*__GMS_GEOCENTRIC_CORDINATES_H__*/
