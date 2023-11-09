#ifndef __GMS_REAL_FFT_2D_H__
#define __GMS_REAL_FFT_2D_H__ 261020151842

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
Real FFT forward/backward Direction 2D class  - definition.
@file:  Real_FFT_2D.h
Wrapper around FFTPACK 5.1 library
@author: Bernard Gingold
@version:  1.0  26/10/2015
@description: file  F77_FFT_DECLARATIONS.h

*/

#include<vector>
#include "FFTPACK_51_DECLARATIONS.h"

namespace gms
{



	

		/*
		@brief wrapper for F77 RFFT2I routine.
		For description of arguments go to F77_FFT_DECLARATIONS.h file
		*/
		void    real_fft_init2D(int ,
			int ,
			std::vector<float> & , 
			int , 
			int );

		/*
		@brief wrapper for F77 RFFT2B routine
		 4th argument should be vector<float> set to size of LDIMxM
		*/
		void    real_fft_backward2D(int ,
			 int ,
			 int ,
			 std::vector<float> &,
			 std::vector<float> &,
			 int ,
			 std::vector<float> &,
			 int ,
			 int );

		/*
		@brief wrapper for F77 RFFT2F routine
		4th argument should be set to size of LDIMxM
		*/
		void real_fft_forward2D(int ,
			int ,
			int ,
			std::vector<float> &,
			std::vector<float> &,
			int ,
			std::vector<float> &,
			int ,
			int );


#include "GMS_real_FFT_2D.inl"
}
#endif /*__GMS_REAL_FFT_2D_H__*/
