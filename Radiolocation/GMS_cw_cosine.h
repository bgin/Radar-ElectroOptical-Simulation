#ifndef __GMS_CW_COSINE_H__
#define __GMS_CW_COSINE_H__

//  For testing purpose only
// Copyright (c) 2020, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
// Cosine Waveform signal class.

#include "WaveformInterface.h"
//#include "WaveformPolarization.h"
#include "Polarization.h"



namespace  gms {
      namespace  radiolocation {


	class CWCosineSignal : public  Waveform
	{

	public:

		CWCosineSignal() = default;
		
		CWCosineSignal( struct CWCosineParams const&);

		CWCosineSignal( const double,  const double,  const double,  const size_t);

		CWCosineSignal( const CWCosineSignal &);

		CWCosineSignal( CWCosineSignal &&);

		virtual ~CWCosineSignal()
		{

		}

		
		CWCosineSignal &  operator=( const CWCosineSignal &);

		CWCosineSignal &  operator=( CWCosineSignal &&);


	        virtual  	void                                    save_to_file(const char*,const char*);


		virtual         void                                    quadrature_components_extraction(std::vector<std::pair<double, double>> &) override;

		virtual         void                                    complex_envelope(std::vector<std::pair<double, double>> &, std::vector<double> &) override;

		virtual         void                                    analytic_signal(const std::vector<double> &) override;

		virtual         std::vector<std::pair<double, double>>   pulse_samples() const override;

		virtual         std::size_t                              pulse_samples_count() const override;


		virtual         void                                    debug_print() const override;

	public:

	        std::vector<std::pair<double, double>> m_cos_signal;

		std::vector<double> m_envelope; // signal envelope.

	        std::vector<double> m_phase; // phase argument.

		double  m_frequency; // carrier frequency

		double  m_efrequency; // envelope signal frequency

		

		double  m_interval; // samples interval

		double  m_init_time; // starting time-point of samples measurement

		size_t  m_samples; // number of std::pair<double,double> samples

		double  m_duration; // signal duration in ns.

		JonesVector m_polarization; // signal linear polarization.

		const   static double  TWO_PI;

		

		
			
		void    initialize( struct  CWCosineParams const&);

		void    initialize( const double,  const double,  const double,  const size_t);

		void    initialize( const CWCosineSignal &);

		void    initialize( CWCosineSignal &&);



	};

	struct CWCosineParams
	{
		 std::function<double(double)> envelope;
		 double cfreq;
		 double efreq;
		 double duration;
		 double interval;
		 double start_time;
		 std::size_t n_samples;
		
	};

    }
}


#endif  __GMS_CW_COSINE_H__

