#ifndef __GMS_REAL_FFT_1D_H__
#define __GMS_REAL_FFT_1D_H__

#include <vector>
#include "GMS_FFTPACK_51_DECLARATIONS.h"

namespace  gms
{



	
                inline
		static	void    real_fft_init1D(int, 
		std::vector<float> &,
		int, 
		int);

                inline
		static  void    real_fft_backward1D(int,
			int, 
			std::vector<float> &, 
			int, 
			std::vector<float> &, 
			int, 
			std::vector<float> &, 
			int, 
			int);
                
                inline
		static  void    real_fft__forward1D(int,
			int,
			std::vector<float> &,
			int,
			std::vector<float> &,
			int, 
			std::vector<float> &, 
			int, 
			int);
	

#include "GMS_real_FFT_1D.inl"
}
#endif /*__GMS_REAL_FFT_1D_H__*/
