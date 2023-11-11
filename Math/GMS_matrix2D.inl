
template<typename T, size_t N> inline Matrix2D<T, N>::Matrix2D()

{
	
	this->m_rows = 0ULL;
	this->m_cols = 0ULL;
	this->m_data = nullptr;

}
template<typename T, size_t N> inline Matrix2D<T, N>::Matrix2D(const size_t rows, const size_t cols)
{
	allocate(rows, cols);
}

template<typename T, size_t N> inline Matrix2D<T, N>::Matrix2D(const size_t rows, const size_t cols, bool dummy )
{
	allocate(rows, cols);
	zero_fill();
}

template<typename T, size_t N> inline Matrix2D<T, N>::Matrix2D(const size_t rows, const size_t cols, const T& value)
{
	allocate(rows, cols);
	fill(value);
}

template<typename T, size_t N> inline Matrix2D<T, N>::Matrix2D(const Matrix2D& rhs)
{
	copy_object(rhs);
}

template<typename T, size_t N> inline Matrix2D<T, N>::Matrix2D(const size_t rows, const size_t cols, T** ptr)
{
	allocate_from(rows, cols, ptr);
}

//----------------------------------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>::Matrix2D(const  std::vector<std::vector<T>>& rhs)
	
{
	allocate_from_vec( rhs);

}
//---------------------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>::Matrix2D(const T(&ar)[N][N])
{
	allocate_from_array(ar);

}
//---------------------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>::Matrix2D(const size_t rows, const size_t cols, T((*&Func)(T)))
{
	allocate_from_func_ptr(rows, cols, Func);

}
//-----------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>::~Matrix2D()
{
	deallocate();

}

template<typename T, size_t N> inline Matrix2D<T, N>::Matrix2D( Matrix2D&& rhs)
{
	
	std::move(*this,rhs);
	rhs.m_cols = 0ULL;
	rhs.m_rows = 0ULL;
	rhs.m_data = nullptr;

}

//Implementation of private member functions
template<typename T, size_t N> inline void Matrix2D<T, N>::allocate(const size_t rows, const size_t cols) 
{
	using namespace std;// Do not pollute global space

		
	this->m_rows = rows;
	this->m_cols = cols;
	const size_t length = this->m_cols * this->m_rows;
	
        this->m_data = new(nothrow) T*[this->m_cols];
	if(nullptr==this->m_data) std::abort();
	this->m_data[0] = new(nothrow) T[length];
	if(nullptr==this->m_data[0]) std::abort();
	
	
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		size_t j = this->m_rows * i;
		this->m_data[i] = &this->m_data[0][j];
	}
}

template<typename T, size_t N> inline void Matrix2D<T, N>::allocate_from(const size_t rows, const size_t cols, T** ptr)
{

	this->m_rows = rows;
	this->m_cols = cols;
	const size_t length = this->m_rows * this->m_cols;
	this->m_data = new(nothrow) T*[this->m_rows];
	if(nullptr==this->m_data) std::abort();
	this->m_data[0] = new(nothrow) T[length];
	if(nullptr==this->m_data[0]) std::abort();
	
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		size_t j = this->m_rows * i;
		this->m_data[i] = &this->m_data[0][j];
	}
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] = ptr[i][j];
		}
	}
}

template<typename T, size_t N> inline void Matrix2D<T, N>::allocate_from_vec( const std::vector<std::vector<T>>& rhs)
{

	this->m_rows = rhs.size();
	this->m_cols = rhs.size();
	const size_t length = this->m_rows * this->m_cols;
	this->m_data = new(nothrow) T*[this->m_rows];
	if(nullptr==this->m_data) std::abort();
	this->m_data[0] = new(nothrow) T[length];
	if(nullptr==this->m_data[0]) std::abort();
	
  	for (register auto i = 0; i != this->m_rows; ++i)
	{
		size_t j0 = this->m_rows * i;
		this->m_data[i] = &this->m_data[0][j0];
	}
	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] = rhs[i][j];
		}
	}
}

template<typename T, size_t N> inline void Matrix2D<T, N>::allocate_from_array(const T(&ar)[N][N])
{
	this->m_rows = N;
	this->m_cols = N;
	const size_t length = this->m_rows * this->m_cols;
	this->m_data = new(nothrow) T*[this->m_rows];
	if(nullptr==this->m_data) std::abort();
	this->m_data[0] = new(nothrow) T[length];
	if(nullptr==this->m_data[0]) std::abort();
	for (register auto i = 0; i != this->m_rows; ++i)
	{
		size_t j0 = this->m_rows * i;
		this->m_data[i] = &this->m_data[0][j0];
	}

	for (register auto i = 0; i != this->m_rows; ++i)
	{

		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] = ar[i][j];
		}
	}

}
//-----------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------//

template<typename T, size_t N> inline void Matrix2D<T, N>::allocate_from_func_ptr(const size_t rows, const size_t cols,
	T((*&Func)(T)))
{



	this->m_rows = rows;
	this->m_cols = cols;
	const size_t length = this->m_rows * this->m_cols;
	this->m_data = new(nothrow) T*[this->m_rows];
	if(nullptr==this->m_data) std::abort();
	this->m_data[0] = new(nothrow) T[length];
        if(nullptr==this->m_data[0]) std::abort();
	
	for (register auto i = 0; i != this->m_rows; ++i)
	{
		size_t j0 = this->m_rows * i;
		this->m_data[i] = &this->m_data[0][j0];
	}
	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] = Func(T);
		}
	}
}
//----------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------//

template<typename T, size_t N> inline void Matrix2D<T, N>::deallocate()
{
	if (this->m_data)
	{
		
		delete[] this->m_data[0];
		delete[] this->m_data;
		this->m_data = nullptr;
	}
}



template<typename T, size_t N> inline void Matrix2D<T, N>::fill(const T& value)
{
	
	        for (register size_t i = 0; i != this->m_rows; ++i)
		{
			for (register size_t j = 0; j != this->m_cols; ++j)
			{
				this->m_data[i][j] = value;
			}
		}
	
}

template<typename T, size_t N> inline void Matrix2D<T, N>::zero_fill()
{
	
	        for (register size_t i = 0; i != this->m_rows; ++i)
		{
			for (register size_t j = 0; j != this->m_cols; ++j)
			{
				this->m_data[i][j] = { 0 };
			}
		}
	
}

template<typename T, size_t N> inline void Matrix2D<T, N>::copy_object(const Matrix2D& rhs)
{
    /*
     * This function is called from within Copy Constructor only.
     */

	this->m_rows = rhs.m_rows;
	this->m_cols = rhs.m_cols;
	const size_t length = this->m_cols * this->m_rows;
	this->m_data = new(nothrow) T*[this->m_cols];
	if(nullptr==this->m_data) std::abort();
	this->m_data[0] = new(nothrow) T[length];
	if(nullptr==this->m_data[0]) std::abort();
	
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		size_t j = this->m_rows * i;
		this->m_data[i] = &this->m_data[0][j];
	}

	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] = rhs.m_data[i][j];
		}
	}
}

template<typename T, size_t N> inline bool Matrix2D<T, N>::check_size(const Matrix2D& rhs)
{


	if ((this->m_rows != rhs.m_rows) || (this->m_cols != rhs.m_cols))
		return true;
	else
		return false;

}

//Member public functions:
	
template<typename T, size_t N>inline T* Matrix2D<T, N>::linearize(const Matrix2D& rhs)
{

	
	const size_t length = rhs.m_rows * rhs.m_cols;
	T* flat_access_ptr = nullptr;
	flat_access_ptr = new(nothrow) T[length];
	if(nullptr==flat_access_ptr) std::abort();
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			T temp = rhs.m_data[i][j];
			flat_access_ptr[i * length + j] = temp;
		}
	}

	return flat_access_ptr;
}

template<typename T, size_t N> inline void Matrix2D<T, N>::resize(const size_t rows_to_add, const size_t cols_to_add)

{

	size_t larger_rows_size = this->m_rows + rows_to_add;
	size_t larger_cols_size = this->m_cols + cols_to_add;
	T** temp = allocateMat2D<T>(this->m_rows, this->m_cols); //must be called with template argument<T>
	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			temp[i][j] = this->m_data[i][j];
		}
	}
	this->allocate_from(larger_rows_size, larger_cols_size,temp);
	deallocateMat2D(temp);
}

template<typename T, size_t N> inline T** Matrix2D<T, N>::get_m_data() const
{
	return m_data;
}

template<typename T, size_t N> inline size_t Matrix2D<T, N>::get_rows() const

{
	return m_rows;
}

template<typename T, size_t N> inline size_t Matrix2D<T, N>::get_cols() const
{
	return m_cols;
}

template<typename T, size_t N> inline T & Matrix2D<T, N>::element_at(const size_t x, const size_t y)
{

	return this->m_data[x][y];
}

template<typename T, size_t N> inline std::unique_ptr<Matrix2D<T,N>> Matrix2D<T, N>::slice(const size_t x0, const size_t x1, const size_t y0, const size_t y1)
{

	//Matrix2D<T, N> object = new Matrix2D<T, N>(this->m_rows, this->m_cols);
	std::unique_ptr<Matrix2D> safe_ptr;
	bool dummy = false;
	Matrix2D* ret_value = new Matrix2D(this->m_rows, this->m_cols,dummy);//Should be fixed : use , return Matrix2D(this->m_rows,this->m_cols);
	for (register size_t i = x0; i != x1; ++i)
	{
		for (register size_t j = y0; j != y1; ++j)
		{
			ret_value->m_data[i][j] = this->m_data[i][j];
		}
	}


	return safe_ptr = (ret_value );//Ownership transfered to unique_ptr;
}

//unique_ptr test initialization
template<typename T, size_t N> inline std::unique_ptr<Matrix2D<T,N>> Matrix2D<T, N>::test(const size_t x0, const size_t x1,
	const size_t y0, const size_t y1)
{
	Matrix2D* mat = new Matrix(10, 10);
	std::unique_ptr<Matrix2D> sp = (mat);
	return sp;
}
//Member Operators

template<typename T, size_t N> inline T& Matrix2D<T, N>::operator()(const size_t x, const size_t y)
{

	return this->m_data[x][y];
}

/*template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator=(const Matrix2D& rhs)
{


	if (this == &rhs) return (*this);
	

	
	this->m_rows = rhs.m_rows;
	this->m_cols = rhs.m_cols;
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] = rhs.m_data[i][j];
		}
	}

	return *this;
}
*/

template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator=(Matrix2D&& rhs)
{

        this->m_rows = rhs.m_rows;
	this->m_cols = rhs.m_cols;
	//std::copy(rhs.m_data, rhs.m_data + (rhs.m_rows*rhs.m_cols), this->m_data);
	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] = rhs.m_data[i][j];
		}
	}
	
	rhs.m_rows = 0;
	rhs.m_cols = 0;
	rhs.m_data = nullptr;
	return *this;
}

//----------------------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator*=(const Matrix2D& A)

{

	this->m_rows = A.m_rows;
	this->m_cols = A.m_cols;
	const size_t Z = this->m_rows;

	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			for (register size_t k = 0; k != Z; ++k)
			{
				this->m_data[i][j] += (this->m_data[i][k] * A.m_data[j][k]);
			}
		}
	}
	return *this;
}
//--------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator+(const T& value)
{

	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] += value;
		}
	}

	return *this;
}
//--------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator*(const T& value)

{
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] *= value;
		}
	}

	return *this;
}
//--------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator-(const T& value)
{
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] -= value;
		}
	}

	return *this;
}
//-------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator/(const T& value)
{

	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] /= value;
		}
	}
	return *this;
}
//---------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator|(const T& value)
{


	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (rgister size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] |= value;
		}
	}

	return *this;
}
//------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator&(const T& value)
{


	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] &= value;
		}
	}

	return *this;
}
//---------------------------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------------------------------------//

template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator^(const T& value)
{

	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] ^= value;
		}
	}
	
	return *this;
}

template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator%(const T& value)
{
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] %= value;
		}
	}

	return *this;
}

template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator<<(const T& value)
{

	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] <<= value;
		}
	}

	return *this;
}

//----------------------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator>>(const T& value)
{

	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] >>= value;
		}
	}

	return *this;
}
//---------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator|=(const Matrix2D& rhs)
{

	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] |= rhs.m_data[i][j];
		}
	}
	return *this;
}
//------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator&=(const Matrix2D& rhs)
{
	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] &= rhs.m_data[i][j];
		}
	}
	return *this;
}
//-------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator^=(const Matrix2D& rhs)
{

	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] ^= rhs.m_data[i][j];
		}
	}

	return *this;
}
//--------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------//

template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator%=(const Matrix2D& rhs)
{

	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] %= rhs.m_data[i][j];
		}
	}
	return *this;
}

template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator<<(const T& value)
{

	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] <<= value;
		}
	}

	return *this;
}

//----------------------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator>>(const T& value)
{

	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] >>= value;
		}
	}

	return *this;
}

//---------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator|=(const Matrix2D& rhs)
{

	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] |= rhs.m_data[i][j];
		}
	}
	return *this;
}

//------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator&=(const Matrix2D& rhs)
{

	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] &= rhs.m_data[i][j];
		}
	}
	return *this;
}

template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator^=(const Matrix2D& rhs)
{

	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] ^= rhs.m_data[i][j];
		}
	}

	return *this;
}

//-----------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator<<=(const Matrix2D& rhs)
{

	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] <<= rhs.m_data[i][j];
		}
	}

	return *this;
}

/--------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator>>=(const Matrix2D& rhs)
{

	for (register auto i = 0; i != this->m_rows; ++i)
	{
		for (register auto j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] >>= rhs.m_data[i][j];
		}
	}

	return *this;
}

template<typename T, size_t N> inline std::vector<std::pair<T,bool>> Matrix2D<T, N>::operator==(const Matrix2D& rhs)
{
	

	if (std::is_integral<T>::value == rhs.m_data[0][0])
	{
		std::pair<T, bool> init_value(0, false);
		std::vector<std::pair<T, bool>> result{ rhs.get_rows(), rhs.get_cols(), init_value };
		bool is_equal = false;
		for (register auto i = 0; i != rhs.get_rows(); ++i)
		{
			for (register auto j = 0; j != rhs.get_cols(); ++j)
			{
				if (this->m_data[i][j] == rhs.m_data[i][j])
				{
					is_equal = true;
					std::pair<T, bool> res(this->m_data[i][j], is_equal);
					result.push_back(res);
				}
			}
		
		}
		return result;
	}
	else
	{
		std::pair<T, bool> init_value(0.0, false);
		std::vector<std::pair<T, bool>> result{ rhs.get_rows(), rhs.get_cols(), init_value };
		bool almost_equal = false;
		for (register auto i = 0; i != rhs.get_rows(); ++i)
		{
			for (register auto j = 0; j != rhs.get_cols(); ++j)
			{
				if (std::fabs(this->m_data[i][j] - rhs.m_data[i][j]) < std::numeric_limits<T>::epsilon())
				{
					almost_equal = true;
					std::pair<T, bool> res(this->m_data[i][j], almost_equal);
					result.push_back(res);
				}
			}
		}
		return result;
	}
	
}

//---------------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator+=(const Matrix2D& rhs)
{
	
	
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] += rhs.m_data[i][j];
		}
	}

	return *this;
}

//-----------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N> operator+(const Matrix2D<T, N>& A, const Matrix2D<T, N>& B)
{


	Matrix2D<T, N> ret_value = Matrix2D<T, N>(A.get_rows(), A.get_cols(), false);
	for (register auto i = 0; i != A.get_rows(); ++i)
	{
		for (register auto j = 0; j != A.get_cols(); ++j)
		{
			ret_value.m_data[i][j] = A.m_data[i][j] + B.m_data[i][j];
		}
	}
	
	return ret_value;
}

//-------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N> operator+(const T& value, const Matrix2D<T, N>& A)
{
	Matrix2D<T, N> ret_value = Matrix2D<T, N>(A.get_rows(), A.get_cols(), false);

	for (register auto i = 0; i != A.m_rows; ++i)
	{
		for (register auto j = 0; j != A.m_cols; ++j)
		{
			ret_value.m_data[i][j] = A.m_data[i][j] + value;
		}
	}
	
	return ret_value;
}

//------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N> operator*(const Matrix2D<T, N>& A, const Matrix2D<T, N>& B)
{

	bool dummy = false;
	Matrix2D<T, N> ret_value = Matrix2D<T, N>(A.get_rows(), A.get_cols(), dummy);

	for (register auto i = 0; i != A.get_rows(); ++i)
	{
		for (register auto j = 0; j != A.get_cols(); ++j)
		{
			for (register auto k = 0; k != A.get_cols(); ++k)
			{
				ret_value.m_data[i][j] += (A.m_data[i][k] * B.m_data[k][j]);
			}
		}
	}
	return ret_value;
}

//-------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N> operator*(const T& value, const Matrix2D<T, N>& rhs)
{
	bool dummy = false;
	Matrix2D<T, N> ret_value = Matrix2D<T, N>(rhs.get_rows(), rhs.get_cols, dummy);

	for (register auto i = 0; i != rhs.get_rows(); ++i)
	{
		for (register auto j = 0; j != rhs.get_cols(); ++j)
		{
			ret_value.m_data[i][j] = rhs.m_data[i][j] * value;
		}
	}
	return ret_value;
}

//-----------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------//
template<typename T, size_t N> inline Matrix2D<T, N>& Matrix2D<T, N>::operator-=(const Matrix2D& rhs)
{
	
	for (register size_t i = 0; i != this->m_rows; ++i)
	{
		for (register size_t j = 0; j != this->m_cols; ++j)
		{
			this->m_data[i][j] -= rhs.m_data[i][j];
		}
	}

	return *this;
}

//-----------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------//
template<typename T, size_t N>  std::ostream& operator<<(std::ostream& os, const Matrix2D<T, N>& A)
{
	for (register auto i = 0; i != A.get_rows(); ++i)
	{
		for (register auto j = 0; j != A.get_cols(); ++j)
		{
			os << "index i = \t" << i << "\tindex j = \t" << j << "\tvalues = " << A.m_data[i][j] << std::endl;
		}
		os << "---------------------------------------------------------------------------------" << std::endl;
	}

	return os;
}

/*
 * Non-member functions
 */

//-----------------------------------------------------------------------------------------------//
template<typename T> inline T** allocateMat2D<T>(const size_t rows, const size_t cols)
{


	const size_t length = rows * cols;
	T** ptr = nullptr;
	ptr = new(nothrow) T*[length];
	if(nullptr==ptr) std::abort();
        ptr[0] = new(nothrow) T[rows];
	if(nullptr==ptr[0]) std::abort();
	
	for (register size_t i = 0; i != rows; ++i)
	{
		size_t j = rows * i;
		ptr[i] = &ptr[0][j];
		
	}
	for (register size_t i = 0; i != rows; ++i)
	{
		for (register size_t j = 0; j != cols; ++j)
		{
			ptr[i][j] = { 0 };
		}
	}

	return ptr;
}

//---------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------//
template<typename T> inline void deallocateMat2D(T**& ptr)
{
	if (ptr)
	{
		delete[] ptr[0];
		delete[] ptr;
		ptr = nullptr;
	}
}

