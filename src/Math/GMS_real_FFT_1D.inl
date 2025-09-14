
void       gms::real_fft_init1D(int N, std::vector<float> &WSAVE, int LENSAV, int IER)
{
	RFFT1I(&N, &WSAVE[0], &LENSAV, &IER);
}

void       gms::real_fft_backward1D(int N, int INC, std::vector<float> &R, int LENR,
	std::vector<float> &WSAVE, int LENSAV, std::vector<float> &WORK, int LENWRK, int IER)
{
	RFFT1B(&N, &INC, &R[0], &LENR, &WSAVE[0], &LENSAV, &WORK[0], &LENWRK, &IER);
}

void        gms::real_fft__forward1D(int N, int INC, std::vector<float> &R, int LENR,
	std::vector<float> &WSAVE, int LENSAV, std::vector<float> &WORK, int LENWRK, int IER)
{
	RFFT1F(&N, &INC, &R[0], &LENR, &WSAVE[0], &LENSAV, &WORK[0], &LENWRK, &IER);
}
