
#ifndef __GMS_MATRIX2D_H__
#define __GMS_MATRIX2D_H__ 110520231053

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

    const unsigned int GMS_MATRIX2D_MAJOR = 1U;
    const unsigned int GMS_MATRIX2D_MINOR = 0U;
    const unsigned int GMS_MATRIX2D_MICRO = 0U;
    const unsigned int GMS_MATRIX2D_FULLVER =
      1000U*GMS_MATRIX2D_MAJOR+
      100U*GMS_MATRIX2D_MINOR+
      10U*GMS_MATRIX2D_MICRO;
    const char * const GMS_MATRIX2D_CREATION_DATE = "11-05-2023 10:53 PM +00200 (THR 11 MAY 2023 GMT+2)";
    const char * const GMS_MATRIX2D_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_MATRIX2D_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_MATRIX2D_DESCRIPTION   = "Simple Matrix 2D.";

}

#include <stdexcept>
#include <memory>
#include <type_traits>
#include <limits>
#include <iostream>
#include "GMS_config.h"


namespace gms {

       template<typename T, size_t N> class Matrix2D {
       
	public:


		//using directives
	using DELETER = void(*)(Matrix2D<T, N>*);
	using Mat = Matrix2D<T, N>;

		//Default Constructor

        inline Matrix2D() noexcept(true);
		// Constructor allocates array of pointers type T*
	inline	Matrix2D(const size_t, const size_t) noexcept(false);

		//Constructor initializes Matrix2D object to 0 by default.
	inline	Matrix2D(const size_t, const size_t, const bool ) noexcept(false);

		//Constructor initializes Matrix2D object by user passed value
	inline	Matrix2D(const size_t, const size_t, const T&) noexcept(false);

	   //Constructor initializes Matrix2D object by user passed pointer T**
	inline Matrix2D(const size_t, const size_t, T**) noexcept(false);

	//Constructor initializes Matrix2D object by passing std::vector<std::vector<T>>& reference
	inline Matrix2D(const std::vector<std::vector<T>>&) noexcept(false);

	//Constructor intializes Matrix2D object from user passed reference to static array.
	inline Matrix2D( const T(&)[N][N]);

	//Constructor initialize Matrix2D object from user passed function pointer.
	inline Matrix2D(const size_t, const size_t ,T((*&)(T))) noexcept(false);

		//Copy constructor
	inline	Matrix2D(const Matrix2D&) noexcept(false);

		//Special Constructor for linearization data to flat array access
	//inline	Matrix2D(const Matrix2D&, bool  );

		//Move Constructor
	inline	Matrix2D( Matrix2D&&) noexcept(true);

		//Destructor
	inline	~Matrix2D();

		//Member functions
	inline T** get_m_data() const;

	//inline T*  get_flat_access_data() const;

	inline size_t get_rows() const;

	inline size_t get_cols() const;

	inline T& element_at(size_t, size_t);

	inline std::unique_ptr<Matrix2D<T,N>> slice(const size_t, const size_t, const size_t, const size_t);

	inline std::unique_ptr<Matrix2D<T, N>> test(const size_t, const size_t, const size_t, const size_t);

	inline T* linearize(const Matrix2D&);

	inline void resize(const size_t, const size_t);//Needs to be tested.
	
	//Member Operators
	inline T& operator()(const size_t, const size_t);

	inline Matrix2D& operator=(const Matrix2D&);

	inline Matrix2D& operator=(Matrix2D&&);

	inline Matrix2D& operator+=(const Matrix2D&);

	inline Matrix2D& operator-=(const Matrix2D&);

	inline Matrix2D& operator*=(const Matrix2D&);

	inline Matrix2D&  operator|=(const Matrix2D&);

	inline Matrix2D& operator&=(const Matrix2D&);

	inline Matrix2D& operator^=(const Matrix2D&);

	inline Matrix2D& operator%=(const Matrix2D&);

	inline Matrix2D& operator<<=(const Matrix2D&);

	inline Matrix2D& operator>>=(const Matrix2D&);

	inline    std::vector<std::pair<T,bool>>  operator==(const Matrix2D&);

	inline bool      operator!=(const Matrix2D&);

	inline Matrix2D& operator+(const T&);

	inline Matrix2D& operator*(const T&);

	inline Matrix2D& operator-(const T&);

	inline Matrix2D& operator/(const T&);

	inline Matrix2D& operator|(const T&);

	inline Matrix2D& operator&(const T&);

	inline Matrix2D& operator^(const T&);

	inline Matrix2D& operator%(const T&);

	inline Matrix2D& operator<<(const T&);

	inline Matrix2D& operator>>(const T&);

	inline bool      operator==(const T&);

	inline bool      operator!=(const T&);
	//friend operators and functions
	

	template<typename T, size_t N> friend std::ostream&    operator<<(std::ostream&, const Matrix2D<T, N>&);

	
	//Scalar versions
	template<typename T, size_t N> friend Matrix2D<T, N>&  Matrix2D_cosine(const T&, const size_t, const size_t);

	template<typename T, size_t N> friend Matrix2D<T, N>&  Matrix2D_sine(const T&, const size_t, const size_t);

	template<typename T, size_t N> friend Matrix2D<T, N>&  Matrix2D_tan(const T&, const size_t, const size_t);

	template<typename T, size_t N> friend Matrix2D<T, N>&  Matrix2D_asin(const T&, const size_t, const size_t);

	template<typename T, size_t N> friend Matrix2D<T, N>&  Matrix2D_acos(const T&, const size_t, const size_t);
	//Vector2D versions, std::vector<std::vector<T>>
	template<typename T, size_t N> friend Matrix2D<T, N>&  Matrix2D_acos_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend Matrix2D<T, N>&  Matrix2D_cosine_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend Matrix2D<T, N>&  Matrix2D_sine_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend Matrix2D<T, N>&  Matrix2D_tan_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend Matrix2D<T, N>&  Matrix2D_asin_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend std::unique_ptr < Matrix2D<T, N>> Matrix2D_asinh_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend std::unique_ptr <Matrix2D<T, N>> Matrix2D_acosh_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend std::unique_ptr < Matrix2D<T, N>> Matrix2D_atanh_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend std::unique_ptr <Matrix2D<T, N>> Matrix2D_sinh_vec(const std::vector<std::vector<T>>&);
	
	template<typename T, size_t N> friend std::unique_ptr <Matrix2D<T, N>> Matrix2D_cosh_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend std::unique_ptr<Matrix2D<T, N>, void(*)(Matrix2D<T, N>*)> Matrix2D_tanh_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend std::unique_ptr<Matrix2D<T, N>, void(*)(Matrix2D<T, N>*)> Matrix2D_erf_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend std::unique_ptr<Matrix2D<T, N>, void(*)(Matrix2D<T, N>*)> Matrix2D_erfc_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend std::unique_ptr<Matrix2D<T, N>, void(*)(Matrix2D<T, N>*)> Matrix2D_lgamma_vec(const std::vector<std::vector<T>>&);

	template<typename T, size_t N> friend std::unique_ptr<Matrix2D<T, N>, void(*)(Matrix2D<T, N>*)> Matrix2D_tgamma_vec(const std::vector<std::vector<T>>&);
	//TO DO Test and implement friend functions with unique_ptr return type.
	template<typename T, size_t N> friend std::unique_ptr < Matrix2D<T, N>> Matrix2D_test(const std::vector<std::vector<T>>&);
	//friend  Matrix2D<T,N>& operator+<>(const Matrix2D& , const Matrix2D& );

	//friend  Matrix2D<T,N>& operator+<>(const T& , const Matrix2D& );
	//class BadArgument{};
	//inline  Matrix2D& operator*=(const Matrix2D&, const Matrix2D&);

	//inline  Matrix2D operator+(const Matrix2D&, const Matrix2D&);
	private:

                size_t m_rows;
                size_t m_cols;
		T** m_data;
		
		
		//T* m_flat_access_data;

	inline	void allocate(const size_t, const size_t);

	inline void allocate_from(const size_t, const size_t, T**);

	inline void allocate_from_vec( const std::vector<std::vector<T>>&);

	inline void allocate_from_array(const T(&)[N][N]);

	inline void allocate_from_func_ptr(const size_t, const size_t, T((*&)(T)));

	inline	void deallocate();

	//inline void deleter();//Possible callee from unique_ptr lambda deleter function
	//inline  void deallocate2();

	inline	void fill(const T&);

	inline  void zero_fill();

	inline void copy_object(const Matrix2D&);

	inline bool check_size(const Matrix2D& );

	

	};

        //Non-member  operators

	  

	  template<typename T, size_t N>  Matrix2D<T, N>&  operator+(const Matrix2D<T, N>& , const Matrix2D<T, N>& );

	  template<typename T, size_t N>  Matrix2D<T, N>&  operator+(const T& , const Matrix2D<T, N>& );

	  template<typename T, size_t N>  Matrix2D<T, N>&  operator*(const Matrix2D<T, N>&, const Matrix2D<T, N>&);

	  template<typename T, size_t N>  Matrix2D<T, N>&  operator*(const T&, const Matrix2D<T, N>&);

	  template<typename T, size_t N>  Matrix2D<T, N>&  operator-(const Matrix2D<T, N>&, const Matrix2D<T, N>&);

	  template<typename T, size_t N>  Matrix2D<T, N>&  operator-(const T&, const Matrix2D<T, N>&);

	  template<typename T, size_t N>  Matrix2D<T, N>&  operator|(const Matrix2D<T, N>&, const Matrix2D<T, N>&);

	//template<typename T, size_t N> inline Matrix2D<T, N> operator+(const Matrix2D<T, N>&,  const Matrix2D<T, N>&);

	// template<typename T, size_t N> inline Matrix2D<T, N>&  operator+(const Matrix2D<T, N>& A, const Matrix2D<T, N>& B);
	//template<typename T, size_t N> inline Matrix2D<T, N> operator*=(const Matrix2D<T, N>&,   const Matrix2D<T, N>&);

	//Namespace memory allocation helper function.
	template<typename T> inline T** allocateMat2D(const size_t, const size_t);

	//Namespace memory deallocation helper function.
	template<typename T> inline void deallocateMat2D(T**&);
        
} // gms
























#endif /*__GMS_MATRIX2D_H__*/
