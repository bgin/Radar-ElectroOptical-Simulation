
#ifndef __GMS_STAT_DISTRIBUTIONS_H__
#define __GMS_STAT_DISTRIBUTIONS_H__

#include <vector>
#include <time.h>
#include <random>
#include <type_traits>
#include <tchar.h>
#include <tuple>
//#include <Windows.h>
//#include "fdlibcppwrapper_utils.h"
namespace gms
{

	template<typename Distribution, typename Engine, typename Value> class UniformDistributions
	{
	public:

		
		/*
		 *  uni_real_distribution(std::vector<Value>&, Distribution, Engine)
		 *  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		 *  Returns:  std::vector<Value>& filled with randomly generated values
		 */

		inline static  std::enable_if<std::is_floating_point<Value>::value, Value>   uni_real_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  uni_int_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/

		inline  static std::enable_if<std::is_integral<Value>::value, Value>   uni_int_distribution(std::vector<Value>&, Distribution , Engine );


	};

	template<typename Distribution, typename Engine, typename Value> class BernoulliDistributions
	{
	public:

		/*
		*  bernoulli_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_integral<Value>::value, Value> bernoulli_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  binomial_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static  std::enable_if<std::is_integral<Value>::value, Value> binomial_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  geometric_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_integral<Value>::value, Value>  geometric_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  neg_binomial_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static  std::enable_if<std::is_integral<Value>::value, Value> neg_binomial_distribution(std::vector<Value>&, Distribution , Engine );
	};

	template<typename Distribution, typename Engine, typename Value> class PoissonDistributions
	{
	public:

		/*
		*  poisson_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_integral<Value>::value, Value>  poisson_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  exponential_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_floating_point<Value>::value, Value> exponential_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  gamma_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_floating_point<Value>::value, Value> gamma_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  weibull_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_floating_point<Value>::value, Value> weibull_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  extrme_value_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_floating_point<Value>::value, Value> extreme_value_distribution(std::vector<Value>&, Distribution , Engine );
	};

	template<typename Distribution, typename Engine, typename Value> class NormalDistributions
	{
	public:

		/*
		*  normal_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_floating_point<Value>::value, Value> normal_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  lognormal_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_floating_point<Value>::value, Value> lognormal_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  chi_squared_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_floating_point<Value>::value, Value> chi_squared_distribution(std::vector<Value>&, Distribution, Engine);

		/*
		*  cauchy_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_floating_point<Value>::value, Value> cauchy_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  fisher_f_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_floating_point<Value>::value, Value> fisher_f_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  student_t_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_floating_point<Value>::value, Value> student_t_distribution(std::vector<Value>&, Distribution , Engine );
	};

	template<typename Distribution, typename Engine, typename Value> class SamplingDistribution
	{
	public:

		/*
		*  discrete_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static std::enable_if<std::is_integral<Value>::value, Value> discrete_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  piecewise_const_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static  std::enable_if<std::is_floating_point<Value>::value, Value> piecewise_const_distribution(std::vector<Value>&, Distribution , Engine );

		/*
		*  piecewise_linear_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, std::uniform_distribution, Engine e
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		inline static  std::enable_if<std::is_floating_point<Value>::value, Value> piecewise_linear_distribution(std::vector<Value>&, Distribution , Engine );
	};


	// UniDistributionImpl class
	template<typename Value> class UniDistributions
	{
	public: 
		/*
		*  uni_int_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, const Value (lower bound), const Value (upper bound), size_t (init seed), const int (switch var
		*  to choose between various engines)
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_integral<T>::value, T> uni_int_distribution(std::vector<Value>&,
		const Value, const Value, size_t, const int);

		/*
		*  uni_real_distribution(std::vector<Value>&, Distribution, Engine)
		*  Arguments: std::vector<Value>&, const Value (lower bound), const Value (upper bound), size_t (init seed), const int (switch var
		*  to choose between various engines)
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> uni_real_distribution(std::vector<Value>&,
			const Value, const Value, size_t, const int);


	};

	// BernoulliDistributionImpl
	template<typename Value> class BernoulliDistributionImpl
	{
	public:

		enum  class RanNumGenerators : unsigned int
		{
			DEFAULT_RANDOM_ENGINE,
			LINEAR_CONGRUENTIAL_ENGINE,
			MINSTD_RANDOM,
			MINSTD_RANDOM0,
			MT199937,
			MT199937_64,
			RAN_LUX24,
			RAN_LUX48,
			KNUTH_B

		};
		
		/*
		*  bernoulli_distribution(std::vector<Value>&, const double, size_t, const int);
		*  Arguments: std::vector<Value>&, const double probability, size_t (init seed), const int (switch var
		*  to choose between various engines)
		*  Returns:  std::vector<Value>& filled with randomly generated values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_integral<T>::value, T> bernoulli_distribution(std::vector<Value>&,
		const double, size_t, const int);

		/*
		*  binomial_distribution(std::vector<Value>&,Value,const double, size_t,const int);
		*  Arguments: std::vector<Value>& values, const Value t, const double probability, size_t (init seed), const int (switch var
		*  to choose between various engines)
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_integral<T>::value, T> binomial_distribution(std::vector<Value>&,
			 Value, const double, size_t, const int);

		/*
		*  geometric_distribution(std::vector<Value>&,Value,const double, size_t,const int);
		*  Arguments: std::vector<Value>& values, const double probability, size_t (init seed), const int (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_integral<T>::value, T> geometric_distribution(std::vector<Value>&,
			const double, size_t, const int);

		/*
		*  negative_binomial_distribution(std::vector<Value>&,Value,const double, size_t,const int);
		*  Arguments: std::vector<Value>& values, const double probability, size_t (init seed), const int (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_integral<T>::value, T> negative_binomial_distribution(std::vector<Value>&,
			const Value, const double, size_t, const int);

		//ENUM CLASS TEST.
		inline static void test22(RanNumGenerators);
		//ENUM CLASS TEST.
	};

	//Poisson distribution implementation class

	template<typename Value> class PoissonDistributionImpl
	{

	public:

		enum class RandomNumberEngines : unsigned int
		{
			DEFAULT_RANDOM_ENGINE,
			LINEAR_CONGRUENTIAL_ENGINE,
			MINSTD_RANDOM,
			MINSTD_RANDOM0,
			MT199937,
			MT199937_64,
			RAN_LUX24,
			RAN_LUX48,
			KNUTH_B
		};

		/*
		* poisson_distribution(std::vector<Value>&,Value,const double, size_t, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values, const Value, size_t (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_integral<T>::value, T> poisson_distribution(std::vector<Value>&,
			const Value, long long, RandomNumberEngines);

		/*
		* exponential_distribution(std::vector<Value>&,Value,const double, size_t, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values, const double probability, size_t (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> exponential_distribution(std::vector<Value>&,
			const Value, unsigned long long, RandomNumberEngines);

		/*
		* gamma_distribution(std::vector<Value>&,Value,const double,const double, unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values, const double alpha,const double beta, unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> gamma_distribution(std::vector<Value>&,
			const Value, const Value, unsigned long long, RandomNumberEngines);

		/*
		* weibull_distribution(std::vector<Value>&,Value,const double,const double, unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values, const double alpha,const double beta, unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> weibull_distribution(std::vector<Value>&,
			const Value, const Value, unsigned long long, RandomNumberEngines);

		/*
		* extreme_value_distribution(std::vector<Value>&,Value,const double,const double, unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values, const double alpha,const double beta, unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> extreme_value_distribution(std::vector<Value>&,
			const Value, const Value, unsigned long long, RandomNumberEngines);
	};

	// Normal Distribution implementation class

	template<typename Value> class NormalDistributionImpl
	{

	public:

		enum class RandomNumberEngines : unsigned int
		{
			DEFAULT_RANDOM_ENGINE,
			LINEAR_CONGRUENTIAL_ENGINE,
			MINSTD_RANDOM,
			MINSTD_RANDOM0,
			MT199937,
			MT199937_64,
			RAN_LUX24,
			RAN_LUX48,
			KNUTH_B
		};

		/*
		* normal_distribution(std::vector<Value>&,Value,const double,const double, unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values, const Value mikron,const Value sigma, unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> normal_distribution(std::vector<Value>&,
			const Value, const Value, unsigned long long, RandomNumberEngines);
		

		/*
		* lognormal_distribution(std::vector<Value>&,Value,const double,const double, unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values, const Value mikron,const Value sigma, unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> lognormal_ditribution(std::vector<Value>&,
			const Value, const Value, unsigned long long, RandomNumberEngines);


		/*
		* chi_squared_distribution(std::vector<Value>& ,const Value, unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values, const Value n, unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> chi_squared_distribution(std::vector<Value>&,
			const Value, unsigned long long, RandomNumberEngines);




		/*
		* cauchy_distribution(std::vector<Value>&,const Value,const Value, unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values, const Value a , const Value b,unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> cauchy_distribution(std::vector<Value>&,
			const Value, const Value, unsigned long long, RandomNumberEngines);

		/*
		* fisher_f_distribution(std::vector<Value>&,const Value,const Value, unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values, const Value m , const Value n,unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> fisher_f_distribution(std::vector<Value>&,
			const Value, const Value, unsigned long long, RandomNumberEngines);

		/*
		* student_t_distribution(std::vector<Value>&,const Value, unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values,  const Value n,unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline static std::enable_if<std::is_floating_point<T>::value, T> student_t_distribution(std::vector<Value>&,
			const Value, unsigned long long, RandomNumberEngines);
	};

	// Sampling DistributionImpl class
	template<typename Value> class SamplingDistributionImpl
	{
	public:

		enum class RandomNumberEngines : unsigned int
		{
			DEFAULT_RANDOM_ENGINE,
			LINEAR_CONGRUENTIAL_ENGINE,
			MINSTD_RANDOM,
			MINSTD_RANDOM0,
			MT199937,
			MT199937_64,
			RAN_LUX24,
			RAN_LUX48,
			KNUTH_B
		};

		/*
		* discrete_distribution(std::vector<Value>&,initializer_list<double> , unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values,  initializer_list<double> ilist,unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		template<typename T = Value> typename inline std::enable_if<std::is_integral<T>::value, T> discrete_distribution(std::vector<Value>&,
			std::initializer_list<double>, unsigned long long, RandomNumberEngines);

		/*
		* piecewise_constant_distribution(std::vector<Value>&,initializer_list<double> , unsigned long long, enum RandomNumberEngines);
		*  Arguments: std::vector<Value>& values,  initializer_list<double> ilist,unsigned long long (init seed), RandomNumberEngines engines (switch var
		*  to choose between various engines).
		*  Returns:  std::vector<Value>& filled with randomly generated and binomially distributed values
		*/
		//TO DO implement it.
	};

	// Test only.
	//extern "C" void test_binomial_distribution_wrapper(unsigned int *, const int);
	// Test only.

	//TO DO:
	// Test implementation.

	template<typename _Ty> class UniformDistributionWrapper
	{

	public:

		// overloaded uniform_int_distribution_wrapper()
		// Arguments: template parameter _Ty * data, const int data_length, const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to uniform_in_distribution, true on success , false on failure.
		inline static bool uniform_int_distribution_wrapper(_Ty *, const int, const int);

		// overloaded uniform_int_distribution_wrapper()
		// Arguments: template parameter _Ty * data, const int data_length, const _Ty a, const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to uniform_in_distribution, true on success , false on failure.
		inline static bool uniform_int_distribution_wrapper(_Ty *, const int, const _Ty, const int);

		// overloaded uniform_int_distribution_wrapper()
		// Arguments: template parameter _Ty * data, const int data_length, const _Ty a,const _Ty b, const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to uniform_in_distribution, true on success , false on failure.
		inline static bool uniform_int_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, const int);


		// overloaded uniform_int_distribution_wrapper()
		// Arguments: template parameter _Ty * data, const int data_length, const _Ty a,const _Ty b, std::size_t seed, const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to uniform_in_distribution, true on success , false on failure.
		inline static bool uniform_int_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, size_t, const int);


		// overloaded uniform_real_distribution_wrapper()
		// Arguments: template parameter _Ty * data, const int data_length, const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to uniform_in_distribution, true on success , false on failure.
		inline static bool uniform_real_distribution_wrapper(_Ty *, const int, const int);

		// overloaded uniform_real_distribution_wrapper()
		// Arguments: template parameter _Ty * data, const int data_length,const _Ty a, const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to uniform_in_distribution, true on success , false on failure.
		inline static bool uniform_real_distribution_wrapper(_Ty *, const int, const _Ty , const int);

		// overloaded uniform_real_distribution_wrapper()
		// Arguments: template parameter _Ty * data, const int data_length,const _Ty a, const _Ty b, const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to uniform_in_distribution, true on success , false on failure.
		inline static bool uniform_real_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, const int);

		// overloaded uniform_real_distribution_wrapper()
		// Arguments: template parameter _Ty * data, const int data_length,const _Ty a, const _Ty b, std::size_t seed , const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to uniform_in_distribution, true on success , false on failure.
		inline static bool uniform_real_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, std::size_t, const int);
	};


	template<typename Value> class BernoulliDistributionsWrapper
	{
	public:
		
		// overloaded bernoulli_distribution_wrapper()
		// Arguments: template parameter Value * data, const int data_length, const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to bernoulli_distribution, true on success , false on failure.
		 inline static bool bernoulli_distribution_wrapper(Value  * , const int, const int);

		// overloaded bernoulli_distribution_wrapper()
		// Arguments: template parameter Value * data, const int data_length, const double p ,const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to bernoulli_distribution, true on success , false on failure.
		 inline static bool bernoulli_distribution_wrapper(Value *, const int, const double, const int);

		// overloaded bernoulli_distribution_wrapper()
		// Arguments: template parameter Value * data, const int data_length, const double p ,std::size_t seed, const int num_of_engine.
		// Returns:  as a argument user passed array of type Value filled with the results of call to bernoulli_distribution, true on success , false on failure.
		 inline static bool bernoulli_distribution_wrapper(Value *, const int, const double, std::size_t, const int);

		  // overloaded binomial_distribution_wrapper()
		 //  Arguments: template parameter Value * data , const int data_length, switch argument const int num_of_engine.
		 //  Returns:   user passed array of type Value filled with results of call to binomial_distribution, true on success , false on failure.
		  inline static bool binomial_distribution_wrapper(Value *, const int,const int);

		  
		  // overloaded binomial_distribution_wrapper()
		  //  Arguments: template parameter Value * data , const int data_length,Value t, switch argument const int num_of_engine.
		  //  Returns:  as a argument user passed array of type Value filled with results of call to binomial_distribution, true on success , false on failure.
		 
		  inline static bool binomial_distribution_wrapper(Value *, const int, Value, const int);


		  // overloaded binomial_distribution_wrapper()
		  //  Arguments: template parameter Value * data , const int data_length,Value t,const double p , const int num_of_engine.
		  //  Returns:  as an argument user passed array of type Value filled with results of call to binomial_distribution , true on success fasle on failure.
		  inline static bool binomial_distribution_wrapper(Value *, const int, Value, const double, const int);

		  // overloaded binomial_distribution_wrapper()
		  //  Arguments: template parameter Value * data , const int data_length,Value t,const double p ,size_t seed, const int num_of_engine.
		  //  Returns:  as an argument user passed array of type Value filled with results of call to binomial_distribution, true on success false on failure.
		  inline static bool binomial_distribution_wrapper(Value *, const int, Value, const double, std::size_t, const int);

		  // overloaded geometric_distribution_wrapper()
		  // Arguments: template parameter Value * data, const int data_length, const int num_of_engine.
		  // Returns:  as a argument user passed array of type Value filled with the results of call to bernoulli_distribution, true on success.
		  inline static bool geometric_distribution_wrapper(Value *, const int, const int);

		  // overloaded geometric_distribution_wrapper()
		  // Arguments: template parameter Value * data, const int data_length,const double p, const int num_of_engine.
		  // Returns:  as a argument user passed array of type Value filled with the results of call to bernoulli_distribution, true on success.
		  inline static bool geometric_distribution_wrapper(Value *, const int, const double, const int);

		  // overloaded geometric_distribution_wrapper()
		  // Arguments: template parameter Value * data, const int data_length,const double p,std::size_t seed, const int num_of_engine.
		  // Returns:  as a argument user passed array of type Value filled with the results of call to bernoulli_distribution, true on success.
		  inline static bool geometric_distribution_wrapper(Value *, const int, const double, std::size_t, const int);

		  // overloaded negative_binomial_distribution_wrapper()
		  // Arguments: template parameter Value * data, const int data_length,const int num_of_engine.
		  // Returns:  as a argument user passed array of type Value filled with the results of call to bernoulli_distribution, true on success.
		  inline static bool negative_binomial_distribution_wrapper(Value *, const int, const int);

		  // overloaded negative_binomial_distribution_wrapper()
		  // Arguments: template parameter Value * data, const int data_length,const Value, const int num_of_engine.
		  // Returns:  as a argument user passed array of type Value filled with the results of call to bernoulli_distribution, true on success.
		  inline static bool negative_binomial_distribution_wrapper(Value *, const int, const Value, const int);

		  // overloaded negative_binomial_distribution_wrapper()
		  // Arguments: template parameter Value * data, const int data_length,const Value,const double p, const int num_of_engine.
		  // Returns:  as a argument user passed array of type Value filled with the results of call to bernoulli_distribution, true on success.
		  inline static bool negative_binomial_distribution_wrapper(Value *, const int, const Value, const double, const int);

		  // overloaded negative_binomial_distribution_wrapper()
		  // Arguments: template parameter Value * data, const int data_length,const Value,const double p,std::size_t seed, const int num_of_engine.
		  // Returns:  as a argument user passed array of type Value filled with the results of call to bernoulli_distribution, true on success.
		  inline static bool negative_binomial_distribution_wrapper(Value *, const int, const Value, const double, std::size_t, const int);


	private:
		const static unsigned int MaxUintVal = 2147483647U;
		//TO DO: to be implemented

		inline static void test_assertion(Value *, const int, const int);

		inline static bool check_arguments(Value *, const int, const int);
	};

	// PoissonDiStributionWrapper class.

	template<typename _Ty> class PoissonDistributionsWrapper
	{

	public:

		 /*
          *  overloaded poisson_distribution_wrapper().
		  *  Arguments: template parameter Value * data, const int data_length, PoissonDistributionImpl::RandomNumberEngines engine.
		  *  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
         */
		inline static bool poisson_distribution_wrapper(_Ty *, const int, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded poisson_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length,_Ty m, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool poisson_distribution_wrapper(_Ty *, const int, _Ty, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded poisson_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length,_Ty m, unsigned long long seed, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool poisson_distribution_wrapper(_Ty *, const int, _Ty, unsigned long long, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded exponential_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool exponential_distribution_wrapper(_Ty *, const int, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded exponential_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length,_Ty lambda, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool exponential_distribution_wrapper(_Ty *, const int, _Ty, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded exponential_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length,_Ty lambda, unsigned  long long seed, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool exponential_distribution_wrapper(_Ty *, const int, _Ty, unsigned long long, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded gamma_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool gamma_distribution_wrapper(_Ty *, const int, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded gamma_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length,_Ty alpha, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool gamma_distribution_wrapper(_Ty *, const int, const _Ty, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded gamma_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length,_Ty alpha,_Ty beta, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool gamma_distribution_wrapper(_Ty *, const int, const _Ty,  const _Ty, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded gamma_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length,_Ty alpha,_Ty beta, unsigned long long seed, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool gamma_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, unsigned long long, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded weibull_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool weibull_distribution_wrapper(_Ty *, const int, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded weibull_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length,_Ty a, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool weibull_distribution_wrapper(_Ty *, const int, const _Ty, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded weibull_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length,_Ty a,_Ty b, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		
		*/
		inline static bool weibull_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded weibull_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length,_Ty a,_Ty b, unsigned long long seed,PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool weibull_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, unsigned long long, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded extreme_value_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool extreme_value_distribution_wrapper(_Ty *, const int, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);


		/*
		*  overloaded extreme_value_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length, const _Ty a, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool extreme_value_distribution_wrapper(_Ty *, const int, _Ty, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded extreme_value_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length, const _Ty a, const _Ty b, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool extreme_value_distribution_wrapper(_Ty *, const int, _Ty, _Ty, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded extreme_value_distribution_wrapper().
		*  Arguments: template parameter Value * data, const int data_length, const _Ty a, const _Ty b, unsigned long long seed, PoissonDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool extreme_value_distribution_wrapper(_Ty *, const int, _Ty, _Ty, unsigned long long, typename PoissonDistributionImpl<_Ty>::RandomNumberEngines);
	private:
		const static unsigned int MaxUintVal = 2147483647U;
		//unsigned long long total_mem = 0ULL;
		//bool ret_val = ::GetPhysicallyInstalledSystemMemory(&total_mem);

		inline static void test_assertion(_Ty *, const int);

		inline static bool report_total_memory(unsigned long long );
	};

	template<typename _Ty> class NormalDistributionsWrapper 
	{
	public:

		/*
		*  overloaded normal_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool normal_distribution_wrapper(_Ty *, const int, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded normal_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length,const _Ty m, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool normal_distribution_wrapper(_Ty *, const int, const _Ty, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded normal_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length,const _Ty m,_Ty s, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool normal_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded normal_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length,const _Ty m,_Ty s, unsigned long long seed, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool normal_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, unsigned long long,
			typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded lognormal_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		*/
		inline static bool lognormal_distribution_wrapper(_Ty *, const int, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded lognormal_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, const _Ty m, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.
		
		*/
		inline static bool lognormal_distribution_wrapper(_Ty *, const int, const _Ty, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded lognormal_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, const _Ty m, const _Ty s, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool lognormal_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded lognormal_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, const _Ty m, const _Ty s,unsigned long long seed, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool lognormal_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, unsigned long long, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded chi_squared_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length,  NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool chi_squared_distribution_wrapper(_Ty *, const int, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded chi_squared_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length,_Ty n,  NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool chi_squared_distribution_wrapper(_Ty *, const int, const _Ty , typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded chi_squared_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length,_Ty n, unsigned long long seed, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool chi_squared_distribution_wrapper(_Ty *, const int, const _Ty, unsigned long long, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded cauchy_squared_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool cauchy_distribution_wrapper(_Ty *, const int, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded cauchy_squared_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length,const _Ty a, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool cauchy_distribution_wrapper(_Ty *, const int, const _Ty, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded cauchy_squared_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length,const _Ty a, const _Ty b, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool cauchy_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded cauchy_squared_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length,const _Ty a, const _Ty b, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool cauchy_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, unsigned long long, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded fisher_f_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool fisher_f_distribution_wrapper(_Ty *, const int, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);
		
		/*
		*  overloaded fisher_f_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, const _Ty m, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool fisher_f_distribution_wrapper(_Ty *, const int, const _Ty, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded fisher_f_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, const _Ty m, const _Ty n, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool fisher_f_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded fisher_f_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, const _Ty m, const _Ty n, unsigned long long seed, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool fisher_f_distribution_wrapper(_Ty *, const int, const _Ty, const _Ty, unsigned long long seed, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded student_t_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length,  NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool student_t_distribution_wrapper(_Ty *, const int, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded student_t_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, const _Ty n, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool student_t_distribution_wrapper(_Ty *, const int, const _Ty, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);

		/*
		*  overloaded student_t_distribution_wrapper().
		*  Arguments: template parameter _Ty * data, const int data_length, const _Ty n,unsigned long long seed, NormalDistributionImpl::RandomNumberEngines engine.
		*  Returns: as an argument user passed array of type Value filled by the result of call to poisson_distribution, true on success, false on failure.

		*/
		inline static bool student_t_distribution_wrapper(_Ty *, const int, const _Ty, unsigned long long, typename NormalDistributionImpl<_Ty>::RandomNumberEngines);
	private:


		inline static void test_assertion(_Ty *, const int);

		inline static bool report_total_memory(unsigned long long);

	};
#include "GMS_stat_distributions.inl"
}
#endif /*__GMS_STAT_DISTRIBUTIONS_H__*/
