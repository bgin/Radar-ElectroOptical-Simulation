#include "CWCosine.h"


// Copyright (c) 2020, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
// Cosine Waveform signal class - implementation.


const double  radiolocation::CWCosineSignal::TWO_PI{ 6.28318530717958647692529 };



gms::radiolocation::CWCosineSignal::CWCosineSignal( struct CWCosineParams const& params)
{
	initialize(params);
}


gms::radiolocation::CWCosineSignal::CWCosineSignal( const double init_time,  const double interval,  const double sfreq,  const size_t n_samples)
{
	initialize(init_time, interval, sfreq, n_samples);
}


gms::radiolocation::CWCosineSignal::CWCosineSignal(const CWCosineSignal &rhs)
{
	initialize(rhs);
}

gms::radiolocation::CWCosineSignal::CWCosineSignal( CWCosineSignal &&rhs)
{
	initialize(rhs);
}


void     radiolocation::CWCosineSignal::initialize(  struct CWCosineParams const& params)
{

	
	this->m_duration = params.duration;
	this->m_frequency = params.cfreq;
	this->m_efrequency = params.efreq;
	this->m_init_time = params.start_time;
	this->m_interval = params.interval;
	this->m_polarization = JonesVector(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 1.0));
	this->m_samples = params.n_samples;

		
		double a_cfreq = this->m_frequency;
		double a_efreq = this->m_efrequency;
		size_t a_samples = this->m_samples;
		double a_stime = this->m_init_time;
		double a_interval = this->m_interval;
		size_t i;
		std::vector<double> a_envelope(a_samples);
		std::vector<std::pair<double, double>> a_cossignal(a_samples);
		std::vector<double> a_phase(a_samples);
		double delta{ 0.0 }; double t{ 0.0 }; double t2{ 0.0 };
		double inv_samples{ 1.0 / static_cast<double>(a_samples) };

		
		for (i = 0; i < a_samples; ++i)
		{
			a_stime += a_interval;
			delta = static_cast<double>(i)  * inv_samples;
			t = a_stime * delta;
			t2 = a_stime + delta;
			//double arg = ((TWO_PI * a_efreq * t2) + (0.25 * TWO_PI));
			a_envelope.operator[](i) = params.envelope((TWO_PI * a_efreq * t2) + (0.25 * TWO_PI));
			//std::pair<double, double> values(a_stime, a_esignal.operator[](i) * ::cos((TWO_PI* a_cfreq * t) + (0.10 * TWO_PI)));
			a_phase.operator[](i) = ((TWO_PI * params.cfreq * t) + 0.10 * TWO_PI);
			a_cossignal.operator[](i).operator=({ t, a_envelope[i] * ::cos(a_phase.operator[](i)) });
		}

		this->m_envelope = std::vector<double>(std::move(a_envelope));
		this->m_phase = std::vector<double>(std::move(a_phase));
		this->m_cos_signal = std::vector<std::pair<double, double>>(std::move(a_cossignal));


	
}

void      gms::radiolocation::CWCosineSignal::initialize( const double init_time, 
                                                          const double interval,  
                                                          const double sfreq,  
                                                          const size_t n_samples)
{


	this->m_duration = 0.0;
	this->m_efrequency = 0.0;
	this->m_samples = n_samples;
	this->m_frequency = sfreq;
	this->m_init_time = init_time;
	this->m_interval = interval;
	this->m_polarization = JonesVector();
	const double HALF_PI{ 1.57079632679489661923132 };
	_Field_size_(this->m_samples) this->m_envelope = std::vector<double>(this->m_samples, 0.0);
	_Field_size_(this->m_samples) this->m_phase = std::vector<double>(this->m_samples, HALF_PI);
	_Field_size_(this->m_samples) this->m_cos_signal = std::vector<std::pair<double, double>>(this->m_samples);

	for (std::size_t i{ 0 }; i != this->m_samples; ++i)
		this->m_cos_signal.operator[](i).operator=({ static_cast<double>(i), this->m_envelope[i] * ::cos(this->m_phase[i]) });
}
		
	



void      gms::radiolocation::CWCosineSignal::initialize(const CWCosineSignal &rhs)
{
	this->m_duration = rhs.m_duration;
	this->m_frequency = rhs.m_frequency;
	this->m_efrequency = rhs.m_efrequency;
	this->m_init_time = rhs.m_init_time;
	this->m_interval = rhs.m_interval;
	this->m_polarization = JonesVector(rhs.m_polarization);
	this->m_samples = rhs.m_samples;
	this->m_envelope = std::vector<double>(rhs.m_envelope);
	this->m_phase = std::vector<double>(rhs.m_phase);
	this->m_cos_signal = std::vector<std::pair<double, double>>(rhs.m_cos_signal);
}

void      gms::radiolocation::CWCosineSignal::initialize( CWCosineSignal &&rhs)
{
	this->m_duration = std::move(rhs.m_duration);
	this->m_efrequency = std::move(rhs.m_efrequency);
	this->m_frequency = std::move(rhs.m_frequency);
	this->m_init_time = std::move(rhs.m_init_time);
	this->m_samples = std::move(rhs.m_samples);
	this->m_polarization = JonesVector(std::move(rhs.m_polarization));
	this->m_interval = std::move(rhs.m_interval);
	this->m_envelope = std::vector<double>(std::move(rhs.m_envelope));
	this->m_phase = std::vector<double>(std::move(rhs.m_phase));
	this->m_cos_signal = std::vector<std::pair<double, double>>(std::move(rhs.m_cos_signal));
}

 gms::radiolocation::CWCosineSignal &  
 radiolocation::CWCosineSignal::operator=( const CWCosineSignal &rhs)
{
	if (this == &rhs) return *this;

	this->m_duration = rhs.m_duration;
	this->m_efrequency = rhs.m_efrequency;
	this->m_frequency = rhs.m_frequency;
	this->m_init_time = rhs.m_init_time;
	this->m_interval = rhs.m_interval;
	this->m_samples = rhs.m_samples;
	this->m_polarization.operator=(rhs.m_polarization);
	this->m_envelope.operator=( rhs.m_envelope);
	this->m_phase.operator=(rhs.m_phase);
	this->m_cos_signal.operator=(rhs.m_cos_signal);
	
	return *this;
	
}

gms::radiolocation::CWCosineSignal &        
radiolocation::CWCosineSignal::operator=( CWCosineSignal &&rhs)
{
	if (this == &rhs) return *this;

	this->m_duration = std::move(rhs.m_duration);
	this->m_efrequency = std::move(rhs.m_efrequency);
	this->m_frequency = std::move(rhs.m_frequency);
	this->m_init_time = std::move(rhs.m_init_time);
	this->m_interval = std::move(rhs.m_interval);
	this->m_polarization.operator=( std::move(rhs.m_polarization));
	this->m_samples = std::move(rhs.m_samples);
	this->m_phase.operator=(std::move(rhs.m_phase));
	this->m_envelope.operator=(std::move(rhs.m_envelope));
	this->m_cos_signal.operator=(std::move(rhs.m_cos_signal));

	return *this;
}

// Perform IQ Decomposition. Employ OPENMP multi-threading for IQ decomposition computation.
void        gms::radiolocation::CWCosineSignal::quadrature_components_extraction(std::vector<std::pair<double, double>> &IQ)
{
	
		
		// Define automatic variables in order to use OPENMP on class members.
		std::size_t a_samples{ this->m_samples };
		std::vector<double> a_cos_part(a_samples);
		std::vector<double> a_sin_part(a_samples);
		std::vector<std::pair<double, double>> a_cosine(this->m_cos_signal);
		double a_interval{ this->m_interval };
		double a_cfreq{ this->m_frequency };
		double inv_samples{ static_cast<double>((1 / a_samples)) };
		double step{ 0.0 }; double delta{ 0.0 }; double t{ 0.0 };
		std::size_t i;
		

		for (i = 0; i < a_samples; ++i)
		{
			step += a_interval;
			delta = static_cast<double>(i)* inv_samples;
			t = step * delta;
			a_cos_part.operator[](i) = 2.0 * ::cos(TWO_PI * a_cfreq * t);
			a_sin_part.operator[](i) = -2.0 * ::sin(TWO_PI * a_cfreq * t);

			IQ.operator[](i).operator=({ a_cosine.operator[](i).second * a_cos_part.operator[](i),
				a_cosine.operator[](i).second * a_sin_part.operator[](i) });
		}

		
}

// Compute complex envelope.
void        gms::radiolocation::CWCosineSignal::complex_envelope(std::vector<std::pair<double, double>> &IQComponents, 
                                                            std::vector<double> &cmplx_envlp)
{
	
		double j_imag{ j().imag() };
		// Not much work here in order to justify OPENMP multi-threading
		for (size_t i = 0; i != this->m_samples; ++i)
			cmplx_envlp[i] = IQComponents.operator[](i).first + (j_imag * IQComponents.operator[](i).second));
	
	
}

void        gms::radiolocation::CWCosineSignal::analytic_signal(const std::vector<double> &cmplx_envlp)
{
		
		std::size_t a_samples = this->m_samples;
		std::vector<double> a_sin_part(a_samples);
		std::vector<double> a_cos_part(a_samples);
		std::vector<std::pair<double, double>> a_cosine(this->m_cos_signal);
		double a_interval = this->m_interval;
		double a_cfreq = this->m_frequency;
		double step{ 0.0 }; double delta{ 0.0 }; double t{ 0.0 };
		
		size_t i;
		double j_imag{ j().imag() };
		

		for (i = 0; i < a_samples; ++i)
		{
			step += a_interval;
			delta = static_cast<double>(i) / static_cast<double>(a_samples);
			t = step * delta;
			a_cos_part.operator[](i) = ::cos(TWO_PI * a_cfreq * t);
			a_sin_part.operator[](i) = j_imag * ::sin(TWO_PI * a_cfreq * t);
			a_cosine.operator[](i).second = (cmplx_envlp.operator[](i) * a_cos_part.operator[](i) +
				cmplx_envlp.operator[](i) * a_sin_part.operator[](i));
		}
		this->m_cos_signal.operator=(a_cosine);

}





void        gms::radiolocation::CWCosineSignal::debug_print() const
{
	std::printf("CWCosineSignal::debug_print\n");
	std::printf("CWCosineSignal object content:\n\n");
	std::printf("&CWCosineSignal=%p\n", this);
	std::printf("&this->m_frequency=%p,this->m_frequency=%.15f\n", &this->m_frequency,this->m_frequency);
	std::printf("&this->m_efrequency=%p,this->m_efreqency=%.15f\n", &this->m_efrequency,this->m_efrequency);
	std::printf("&this->m_duration=%p,this->m_duration=%.15f\n", &this->m_duration,this->m_duration);
	std::printf("&this->m_samples=%p,this->m_samples=%ull\n", &this->m_samples,this->m_samples);
	std::printf("&this->m_init_time=%p,this->m_init_time=%.15f\n", &this->m_init_time,this->m_init_time);
	std::printf("&this->m_interval=%p,this->m_interval=%.15f\n", &this->m_interval,this->m_interval);
	std::printf("r(t): |   s(t):   |      phi(t):    t:\n");
	for (size_t i = 0; i != this->m_samples; ++i)
	{
		std::printf(" %.15f, %.15f, %.15f,   %.9f\n", this->m_envelope[i], this->m_cos_signal[i].second,this->m_phase[i], this->m_cos_signal[i].first);
	}
	std::printf("End of CWCosineSignal object dump\n");
}



void        gms::radiolocation::CWCosineSignal::save_to_file(const char* fname1,const char* fname2)
{


	FILE *fp1, *fp2;
	if (fopen_s(&fp1, fname1, "wt") != 0)
	{
		std::printf("Failed to open file:%s, ...returning\n", fname1);
		return;
	}
	else
	{
		std::printf("Started writing dump of this->m_cos_signal to file:%s\n ", fname1);
		for (size_t i = 0; i != this->m_samples; ++i)
			fprintf(fp1, "%.9f,%.9f\n", this->m_cos_signal[i].first, this->m_cos_signal[i].second);

		std::printf("Finished writing dump of this->m_cos_signal to file:%s\n", fname1);
		fclose(fp1);
	}

	if (fopen_s(&fp2, fname2, "wt") != 0)
	{
		std::printf("Failed to open file:%s, ...returning\n", fname2);
		return;
	}
	else
	{
		std::printf("Started writing dump of this->m_envelope to file:%s\n", fname2);
		for (size_t i = 0; i != this->m_samples; ++i)
			fprintf(fp2, "%.9f,%.9f\n", this->m_cos_signal[i].first, this->m_envelope[i]);

		std::printf("Finished writing dump of this->m_envelope to file:%s\n", fname2);
		fclose(fp2);
	}
}


