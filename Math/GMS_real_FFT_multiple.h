#ifndef __GMS_REAL_FFT_MULTIPLE_H__
#define __GMS_REAL_FFT_MULTIPLE_H__

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
/*
Real FFT forward/backward Direction Multiple  class  - definition.
Wrapper around FFTPACK 5.1 library
@file Real_FFT_Multiple.h
@author: Bernard Gingold
@version:  1.0  26/10/2015
@description F77_FFT_DECLARATIONS.h

*/

#include <vector>
#include "FFTPACK_51_DECLARATIONS.h"

namespace   gms
{

	


		/*
		@brief wrapper for F77 RFFTMI routine
		*/
		
		inline
		static
		void    real_fft_multiple_init(int,
			std::vector<float> &,
			int,
			int);


		/*
		@brief  wrapper for F77 RFFTMB routine.

		*/
		inline
		static
		void     real_fft_multiple_backward(int,
			int,
			int,
			int,
			std::vector<float> &,
			int,
			std::vector<float> &,
			int,
			std::vector<float> &,
			int,
			int);

		/*
		@brief wrapper for F77 RFFTMF routine.
		*/
		inline 
		static
		void     real_fft_multiple_forward(int,
			int,
			int,
			int,
			std::vector<float> &,
			int,
			std::vector<float> &,
			int,
			std::vector<float> &,
			int,
			int);


	

#include "GMS_real_FFT_multiple.inl"
}
#endif/*__GMS_REAL_FFT_MULTIPLE_H__*/
