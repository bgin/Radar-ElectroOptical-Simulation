




void      gms::real_fft_multiple_init(int N, std::vector<float> &WSAVE,
	int LENSAV, int IER)
{
	RFFTMI(&N, &WSAVE[0], &LENSAV, &IER);
}


void     gms::real_fft_multiple_backward(int LOT, int JUMP, int N, int INC,
	std::vector<float> &R, int LENR, std::vector<float> &WSAVE, int LENSAV,
	std::vector<float> &WORK, int LENWRK, int IER)
{
	RFFTMB(&LOT, &JUMP, &N, &INC, &R[0], &LENR, &WSAVE[0], &LENSAV, &WORK[0], &LENWRK, &IER);
}


void       gms::real_fft_multiple_forward(int LOT, int JUMP, int N, int INC,
	std::vector<float> &R, int LENR, std::vector<float> &WSAVE, int LENSAV,
	std::vector<float> &WORK, int LENWRK, int IER)
{
	RFFTMF(&LOT, &JUMP, &N, &INC, &R[0], &LENR, &WSAVE[0], &LENSAV, &WORK[0], &LENWRK, &IER);
}
