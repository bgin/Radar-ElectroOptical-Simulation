
/*
Real FFT forward/backward Direction 1D class  - wrapper implementation.
Wrapper around FFTPACK 5.1 library
@author: Bernard Gingold
@version:  1.0  26/10/2015

*/

void          gms::real_fft_init2D(int L, int M, std::vector<float> &WSAVE, int LENSAV, int IER)
{
	RFFT2I(&L, &M, &WSAVE[0], &LENSAV, &IER);
}

void          gms::real_fft_backward2D(int LDIM, int L, int M, std::vector<float> &R,
	std::vector<float> &WSAVE, int LENSAV, std::vector<float> &WORK, int LENWRK, int IER)
{
	RFFT2B(&LDIM, &L, &M, &R[0], &WSAVE[0], &LENSAV, &WORK[0], &LENWRK, &IER);
}

void         gms::real_fft_forward2D(int LDM, int L, int M, std::vector<float> &R,
	std::vector<float> &WSAVE, int LENSAV, std::vector<float> &WORK, int LENWRK, int IER)
{
	RFFT2F(&LDM, &L, &M, &R[0], &WSAVE[0], &LENSAV, &WORK[0], &LENWRK, &IER);
}
