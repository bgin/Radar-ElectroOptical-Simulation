
template<typename Distribution, typename Engine, typename Value> inline  std::enable_if<std::is_integral<Value>::value, Value>
 gms::UniformDistributions<Distribution, Engine, Value>::uni_int_distribution(std::vector<Value>& values,
Distribution d, Engine e)
{
  
	 for (auto& elems : values)
		 elems = d(e);
	 
	 std::enable_if<std::is_integral<Value>::value, Value> s;
	 return s;
}
//------------------------------------------------------------------------------------------------------------------//
 template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
	 gms::UniformDistributions<Distribution, Engine, Value>::uni_real_distribution(std::vector<Value>& values,
		 Distribution d, Engine e)
	 {
		 for (auto& elems : values)
			 elems = d(e);
	 }
//-----------------------------------------------------------------------------------------------------------------//
template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_integral<Value>::value, Value>
	gms::BernoulliDistributions<Distribution, Engine, Value>::bernoulli_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);

		std::enable_if<std::is_integral<Value>::value, Value> s;
		return s;
	}
//---------------------------------------------------------------------------------------------------------------//
template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_integral<Value>::value, Value>
	gms::BernoulliDistributions<Distribution, Engine, Value>::binomial_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);
	}
//-------------------------------------------------------------------------------------------------------------//
	template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_integral<Value>::value, Value>
		gms::BernoulliDistributions<Distribution, Engine, Value>::geometric_distribution(std::vector<Value>& values,
			Distribution d, Engine e)
		{
			for (auto& elems : values)
				elems = d(e);
		}
//-----------------------------------------------------------------------------------------------------------//
template<typename Distribution, typename Engine , typename Value> inline std::enable_if<std::is_integral<Value>::value, Value>
	gms::BernoulliDistributions<Distribution, Engine, Value>::neg_binomial_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);
	}
//----------------------------------------------------------------------------------------------------------//
template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_integral<Value>::value, Value>
	gms::PoissonDistributions<Distribution, Engine, Value>::poisson_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);
	}
//--------------------------------------------------------------------------------------------------------//
	template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point <Value>::value, Value>
	gms::PoissonDistributions<Distribution, Engine, Value>::exponential_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);
	}
//------------------------------------------------------------------------------------------------------//
template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
	gms::PoissonDistributions<Distribution, Engine, Value>::gamma_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);
	}
//----------------------------------------------------------------------------------------------------//
	template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
		gms::PoissonDistributions<Distribution, Engine, Value>::weibull_distribution(std::vector<Value>& values,
			Distribution d, Engine e)
		{
			for (auto& elems : values)
				elems = d(e);
		}
//---------------------------------------------------------------------------------------------------//
template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
	gms::PoissonDistributions<Distribution, Engine, Value>::extreme_value_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);
	}
//-------------------------------------------------------------------------------------------------//
template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
	gms::NormalDistributions<Distribution, Engine, Value>::cauchy_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);
	}
//------------------------------------------------------------------------------------------------//
template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
	gms::NormalDistributions<Distribution, Engine, Value>::chi_squared_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);
	}
//------------------------------------------------------------------------------------------------//
template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
	gms::NormalDistributions<Distribution, Engine, Value>::fisher_f_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);
	}
//----------------------------------------------------------------------------------------------//
template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
	gms::NormalDistributions<Distribution, Engine, Value>::lognormal_distribution(std::vector<Value>& values,
		Distribution d, Engine e)
	{
		for (auto& elems : values)
			elems = d(e);
	}
//---------------------------------------------------------------------------------------------//
	template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
		gms::NormalDistributions<Distribution, Engine, Value>::normal_distribution(std::vector<Value>& values,
			Distribution d, Engine e)
		{
			for (auto& elems : values)
				elems = d(e);
		}
//--------------------------------------------------------------------------------------------//
		template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
			gms::NormalDistributions<Distribution, Engine, Value>::student_t_distribution(std::vector<Value>& values,
				Distribution d, Engine e)
			{
				for (auto& elems : values)
					elems = d(e);
			}
//-------------------------------------------------------------------------------------------//
	template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_integral<Value>::value, Value>
		gms::SamplingDistribution<Distribution, Engine, Value>::discrete_distribution(std::vector<Value>& values,
			Distribution d, Engine e)
		{
			for (auto& elems : values)
				elems = d(e);
		}
//-------------------------------------------------------------------------------------------//
		template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
			gms::SamplingDistribution<Distribution, Engine, Value>::piecewise_const_distribution(std::vector<Value>& values,
				Distribution d, Engine e)
			{
				for (auto& elems : values)
					elems = d(e);
			}
//-------------------------------------------------------------------------------------------//
			template<typename Distribution, typename Engine, typename Value> inline std::enable_if<std::is_floating_point<Value>::value, Value>
				gms::SamplingDistribution<Distribution, Engine, Value>::piecewise_linear_distribution(std::vector<Value>& values,
					Distribution d, Engine e)
				{
					for (auto& elems : values)
						elems = d(e);
				}


// template<typename Value> class UniDistribution

	  template<typename Value>	template<typename T = Value>   inline std::enable_if<std::is_integral<T>::value,T> gms::UniDistributions<Value>::uni_int_distribution(
					std::vector<Value>& values, const Value min_val, const Value max_val, size_t init_seed, const int choose_from_engine)
				{

					

					long long mod_int_seed = clock();
					init_seed = init_seed ^ mod_int_seed;

					switch (choose_from_engine)
					{
					case 0:
					{
							  std::default_random_engine dre(mod_int_seed);
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(dre);
					}
						break;
					case 1:
					{

							  const unsigned int val = 10;
							  const  unsigned int val2 = 65553;
							  const  unsigned int val3 = 2147483647;
							  std::linear_congruential_engine<Value, val, val2, val3> lce(mod_int_seed);
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(lce);
					}
						break;
					case 2:
					{

							  std::minstd_rand mre;
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(mre);
					}
						break;
					case 3:
					{
							  std::minstd_rand0 mre0(mod_int_seed);
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(mre0);

					}
						break;
					case 4:
					{

							  std::mt19937 mte(mod_int_seed);
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(mte);
					}
						break;
					case 5:
					{

							  std::mt19937_64 mte_64(mod_int_seed);
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(mte_64);
					}
						break;
					case 6:
					{
							  std::ranlux24 rl24(mod_int_seed);
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(rl24);
					}
						break;
					case 7:
					{            //TO DO:
							     // Implement in different function. 
							  /*std::tuple<std::ranlux24, std::ranlux24_base, std::ranlux48> ran_engines;
							  
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							 // std::get<0>(ran_engines);
							  for (auto& elems : values)
								  elems = distribution(std::get<0>(ran_engines));
							 
								  for (auto j = 0; j != values.size(); ++j)
								  {
									  values[j] = distribution(std::get<0>(ran_engines));
								  }*/
							  std::knuth_b kbe(mod_int_seed);
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(kbe);
							 
					}
						break;
					case 8:
					{
							  std::ranlux24_base rl24_base(mod_int_seed);
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(rl24_base);
					}
						break;
					case 9:
					{
							  std::ranlux48 rl48(mod_int_seed);
							  std::uniform_int_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(rl48);
					}
						break;
					default:
					{
							  // throw std::runtime_error(std::basic_string<TCHAR>("Invalid switch argument has been passed, only values [0:9] are accepted\n"));
							   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:9] are accepted\n"));
					}

					}

					
				}

				template<typename Value> template<typename T = Value> inline std::enable_if<std::is_floating_point<T>::value, T> gms::UniDistributions<Value>::uni_real_distribution(
					std::vector<Value>& values, const Value min_val, const Value max_val, size_t init_seed, const int choose_from_engine)
				{

					long long seed = clock();
					init_seed = init_seed ^ seed;

					switch (choose_from_engine)
					{
					case 0:
					{
							  std::default_random_engine dre(init_seed);
							  std::uniform_real_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(dre);
					}
						break;
					case 1:
					{
							  const unsigned int val1 = 56U;
							  const unsigned int val2 = 65536U;
							  const unsigned int val3 = 2147483647U;
							  std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
							  std::uniform_real_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(lce);
							  /*for (auto i = 0; i != values.size(); ++i)
							  {
								  //unsigned int rand_val = lce.operator()();
								  values[i] = distribution(lce);
							  }*/
					}
						break;
					case 2:
					{
							  std::minstd_rand msr(init_seed);
							  std::uniform_real_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(msr);
					}
						break;
					case 3:
					{
							  std::minstd_rand0 msr0(init_seed);
							  std::uniform_real_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(msr0);
					}
						break;
					case 4:
					{
							  std::mt19937 mte(init_seed);
							  std::uniform_real_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(mte);
					}
						break;
					case 5:
					{
							  std::mt19937_64 mte64(init_seed);
							  std::uniform_real_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(mte64);
					}
						break;
					case 6:
					{
							  std::ranlux24 rl24(init_seed);
							  std::uniform_real_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(rl24);
					}
						break;
					case 7:
					{
							  std::ranlux24_base rl24base(init_seed);
							  std::uniform_real_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(rl24base);
					}
						break;
					case 8:
					{
							  std::knuth_b kbe(init_seed);
							  std::uniform_real_distribution<Value> distribution(min_val, max_val);
							  for (auto& elems : values)
								  elems = distribution(kbe);
					}
						break;
					default:
					{
							  // throw std::runtime_error(std::basic_string<TCHAR>("Invalid switch argument has been passed, only values [0:9] are accepted\n"));
							   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
					}

					}
				}

				template<typename Value> template<typename T = Value> inline std::enable_if<std::is_integral<T>::value, T> gms::BernoulliDistributionImpl<Value>::bernoulli_distribution(
					std::vector<Value>& values, const double probability, size_t init_seed, const int choose_from_engine)
				{

					long long seed = clock();
					init_seed = init_seed ^ seed;

					switch (choose_from_engine)
					{
					case 0:
					{
							  std::default_random_engine dre(init_seed);
							  std::bernoulli_distribution bd{ probability };
							  for (auto& elems : values)
								  elems = bd(dre);
					}
						break;
					case 1:
					{
							  const unsigned int val1 = 56U;
							  const unsigned int val2 = 65536U;
							  const unsigned int val3 = 2147483647U;
							  std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
							  std::bernoulli_distribution bd{ probability };
							  for (auto& elems : values)
								  elems = bd(lce);
					}
						break;
					case 2:
					{
							  std::minstd_rand msr(init_seed);
							  std::bernoulli_distribution bd{ probability };
							  for (auto& elems : values)
								  elems = bd(msr);
					}
						break;
					case 3:
					{
							  std::minstd_rand0 msr0(init_seed);
							  std::bernoulli_distribution bd{ probability };
							  for (auto& elems : values)
								  elems = bd(msr0);
					}
						break;
					case 4:
					{
							  std::mt19937 mt(init_seed);
							  std::bernoulli_distribution bd{ probability };
							  for (auto& elems : values)
								  elems = mt(bd);
					}
						break;
					case 5:
					{
							  std::mt19937_64 mt_64(init_seed);
							  std::bernoulli_distribution bd{ probability };
							  for (auto& elems : values)
								  elems = mt_64(bd);
					}
						break;
					case 6:
					{
							  std::ranlux24 rl24(init_seed);
							  std::bernoulli_distribution bd{ probability };
							  for (auto& elems : values)
								  elems = bd(rl24);
					}
						break;
					case 7:
					{
							  std::ranlux24_base rl24_base(init_seed);
							  std::bernoulli_distribution bd{ probability };
							  for (auto& elems : values)
								  elems = bd.operator()(rl24_base);
					}
						break;
					case 8:
					{
							  std::knuth_b kb(init_seed);
							  std::bernoulli_distribution bd{ probability };
							  for (auto& elems : values)
								  elems = bd.operator()(kb);
					}
						break;
					default:
					{
							   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
					}

					}
				}

					template<typename Value> template<typename T = Value> inline std::enable_if<std::is_integral<T>::value, T> BernoulliDistributionImpl<Value>::binomial_distribution(std::vector<Value>& values,
						Value t, const double probablility, size_t init_seed, const int choose_from_engine)
					{

						long long rand_seed = clock();
						init_seed = init_seed ^ rand_seed;
						switch (choose_from_engine)
						{
						case 0:
						{
								  // binomial distribution initialized by default random engine.
								  std::default_random_engine dre(init_seed);
								  std::binomial_distribution<Value,double> bd(t, probablility);
								  for (auto& elems : values)
									  elems = bd.operator()(dre); //operator call, should speed compiling? Needs to be tested.
						}
							break;
						case 1:
						{
								  // Cannot use constexpr because of weird bug : 'constexpr' is undefined.
								  // This bug only persists when either ICC or MSVC compiles MS Native CPP Test Framework.
								  const unsigned int val1 = 56U;
								  const unsigned int val2 = 65536U;
								  const unsigned int val3 = 2147483647U;
								  std::linear_congruential_engine<Value, val1, val2, val3> lce(init_seed);
								  std::binomial_distribution<Value, double> bd(t, probablility);
								  for (auto& elems : values)
									  elems = bd.operator()(lce);
						}
							break;
						case 2:
						{
								  // std::binomial_distribution initialized by std::minstd_rand engine adapter.
								  std::minstd_rand  msr(init_seed);
								  std::binomial_distribution<Value, double> bd(t, probablility);
								  for (auto& elems : values)
									  elems = bd.operator()(msr);
						}
							break;
						case 3:
						{          // std::binomial_distribution initialized by std::minstd_rand0 engine adapter.
								  std::minstd_rand0 msr0(init_seed);
								  std::binomial_distribution<Value, double> bd(t, probablility);
								  for (auto& elems : values)
									  elems = bd.operator()(msr0);
						}
							break;
						case 4:
						{
								  //std::binomial_distribution initialized by std::mt19937 engine adapter.
								  std::mt19937 mt(init_seed);
								  std::binomial_distribution<Value, double> bd(t, probablility);
								  for (auto& elems : values)
									  elems = bd.operator()(mt);
						}
							break;
						case 5:
						{
								  //std::binomial_distribution initialized by std::mt19937_64 engine adapter.
								  std::mt19937_64 mt64(init_seed);
								  std::binomial_distribution<Value, double> bd(t, probablility);
								  for (auto& elems : values)
									  elems = bd.operator()(mt64);
						}
							break;
						case 6:
						{
								  //std::binomial_distribution initialized by std::ranlux24 engine adapter.
								  std::ranlux24 rl24(init_seed);
								  std::binomial_distribution<Value, double> bd(t, probablility);
								  for (auto& elems : values)
									  elems = bd.operator()(rl24);

						}
							break;
						case 7:
						{
                                 //std::binomial_distribution initialized by std::ranlux48 engine adapter.
								  std::ranlux48 rl48(init_seed);
								  std::binomial_distribution<Value, double> bd(t, probablility);
								  for (auto& elems : values)
									  elems = bd.operator()(rl48);
						}
							break;
						case 8:
						{
								  //std::binomial_distribution initialized by std::knuth_b engine adapter.
								  std::knuth_b kbe(init_seed);
								  std::binomial_distribution<Value, double> bd(t, probablility);
								  for (auto& elems : values)
									  elems = bd.operator()(kbe);
						}
							break;
						default:
						{
								   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
						}

						}

					}

					template<typename Value> template<typename T = Value> inline std::enable_if<std::is_integral<T>::value, T> BernoulliDistributionImpl<Value>::geometric_distribution(std::vector<Value>& values,
						const double probability, size_t init_seed, const int choose_from_engine)
					{

						long long seed = clock();
						init_seed = init_seed ^ seed;
						switch (choose_from_engine)
						{
						case 0:
						{
								  // std::geometric_distribution iitialized by std::default_random_engine engine.
								  std::default_random_engine dre(init_seed);
								  std::geometric_distribution<Value, double> gd(probability);
								  for (auto& elems : values)
									  elems = gd.operator()(dre);
						}
							break;
						case 1:
						{
								  const unsigned int val1 = 56U;
								  const unsigned int val2 = 65536U;
								  const unsigned int val3 = 2147483647U;
								  std::linear_congruential_engine<Value, val1, val2, val3> lce(init_seed);
								  std::geometric_distribution<Value, double> gd(probability);
								  for (auto& elems : values)
									  elems = gd.operator()(lce);
						}
							break;
						case 2:
						{
								  // std::geometric_distribution initialized by std::minstd_rand engine adapter.
								  std::minstd_rand msr(init_seed);
								  std::geometric_distribution<Value, double> gd(probability);
								  for (auto& elems : values)
									  elems = gd.operator()(msr);
						}
							break;
						case 3:
						{
							 // std::geometric_distribution initialized by std::minstd_rand0 engine adapter.
								  std::minstd_rand0 msr0(init_seed);
								  std::geometric_distribution<Value, double> gd(probability);
								  for (auto& elems : values)
									  elems = gd.operator()(msr0);
						}
							break;
						case 4:
						{
								  // std::geometric_distribution initialized by std::m19937 engine adapter.
								  std::mt19937 mt(init_seed);
								  std::geometric_distribution<Value, double> gd(probability);
								  for (auto& elems : values)
									  elems = gd.operator()(mt);
						}
							break;
						case 5:
						{
								  // std::geometric_distribution initialized by std::mt19937_64 engine adapter.
								  std::mt19937_64 mt64(init_seed);
								  std::geometric_distribution<Values, double> gd(probability);
								  for (auto& elems : values)
									  elems = gd.operator()(mt64);
						}
							break;
						case 6:
						{
								  // std::geometric_distribution initialized by std::ranlux24 engine adapter.
								  std::ranlux24 rl24(init_seed);
								  std::geometric_distribution<Value, double> gd(probability);
								  for (auto& elems : values)
									  elems = gd.operator()(rl24);
						}
							break;
						case 7:
						{
								  // std::geometric_distribution initialized by std::ranlux48 engine adapter.
								  std::ranlux48 rl48(init_seed);
								  std::geometric_distribution<Value, double> gd(probability);
								  for (auto& elems : values)
									  elems = gd.operator()(rl48);
						}
							break;
						case 8:
						{
								  // std::geometric_distribution initialized by std::knuth_b engine adaptor.
								  std::knuth_b kbe(init_seed);
								  std::geometric_distribution<Value, double> gd(probability);
								  for (auto& elems : values)
									  elems = gd.operator()(kbe);
						}
							break;
						default:
						{
								   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
						}

						}

					}

					template<typename Value> template<typename T = Value> inline std::enable_if<std::is_integral<T>::value, T> BernoulliDistributionImpl<Value>::negative_binomial_distribution(std::vector<Value>& values,
						const Value k, const double probability, size_t init_seed, const int choose_from_engine)
					{

						long long seed = 0ULL;
						init_seed = init_seed ^ (seed = clock());
						switch (choose_from_engine)
						{
						case 0:
						{
								  std::default_random_engine dre(init_seed);
								  std::negative_binomial_distribution<Value> nbd(k, probability);
								  for (auto& elems : values)
									  elems = nbd.operator()(dre);
						}
							break;
						case 1:
						{
								 const unsigned int val1 = 56U;
								  const unsigned int val2 = 65536U;
								  const unsigned int val3 = 2147483647U;
								 
								  std::linear_congruential_engine<Value, val1, val2, val3 > lce(init_seed);
								  std::negative_binomial_distribution<Value> nbd(k, probability);
								  for (auto& elems : values)
									  elems = nbd.operator()(lce);

						}
							break;
						case 2:
						{
								  std::minstd_rand msr(init_seed);
								  std::negative_binomial_distribution<Value> nbd(k, probability);
								  for (auto& elems : values)
									  elems = nbd.operator()(msr);
						}
							break;
						case 3:
						{
								  std::minstd_rand0 msr0(init_seed);
								  std::negative_binomial_distribution<Value> nbd(k, probability);
								  for (auto& elems : values)
									  elems = nbd.operator()(msr0);
						}
							break;
						case 4:
						{
								  std::mt19937 mt(init_seed);
								  std::negative_binomial_distribution<Value> nbd(k, probability);
								  for (auto& elems : values)
									  elems = nbd.operator()(mt);
						}
							break;
						case 5:
						{
								  std::mt19937_64 mt64(init_seed);
								  std::negative_binomial_distribution<Value> nbd(k, probability);
								  for (auto& elems : values)
									  elems = nbd.operator()(mt64);
						}
							break;
						case 6:
						{
								  std::ranlux24 rl24(init_seed);
								  std::negative_binomial_distribution<Value> nbd(k, probability);
								  for (auto& elems : values)
									  elems = nbd.operator()(rl24);
						}
							break;
						case 7:
						{
								  std::ranlux48 rl48(init_seed);
								  std::negative_binomial_distribution<Value> nbd(k, probability);
								  for (auto& elems : values)
									  elems = nbd.operator()(rl48);
						}
							break;
						case 8:
						{
								  std::knuth_b kbe(init_seed);
								  std::negative_binomial_distribution<Value> nbd(k, probability);
								  for (auto& elems : values)
									  elems = nbd.operator()(kbe);
						}
							break;
						default:
						{
								   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
						}

						}

					}

					template<typename Value> template<typename T = Value> inline std::enable_if<std::is_integral<T>::value, T> PoissonDistributionImpl<Value>::poisson_distribution(std::vector<Value>& values,
						const Value m, long long init_seed, RandomNumberEngines engines)
					{

						long long seed = 0LL;
						init_seed = init_seed ^ (seed = clock());
						switch (engines)
						{
						case PoissonDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
						{
					   std::default_random_engine dre(init_seed);
			                   std::poisson_distribution<Value, double> pd(m);
					   for (auto& elems : values)
  					        elems = pd.operator()(dre);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
						{
						const unsigned int val1 = 56U;
  				                const unsigned int val2 = 65536U;const unsigned int val3 = 2147483647U;
     																					 
						std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
   					        std::poisson_distribution<Value, double> pd(m);
						for (auto& elems : values)
    						     elems = pd.operator()(lce);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
						{
						   std::minstd_rand msr(init_seed);
						   std::poisson_distribution<Value, double> pd(m);
						   for (auto& elems : values)
   						        elems = pd.operator()(msr);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
						{
						     std::minstd_rand0 msr0(init_seed);
						     std::poisson_distribution<Value, double> pd(m);
						     for (auto& elems : values)
  							 elems = pd.operator()(msr0);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MT199937:
						{
						  std::mt19937 mt(init_seed);
						  std::poisson_distribution<Value, double> pd(m);
						  for (auto& elems : values)
						      elems = pd.operator()(mt);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
						{
						 std::mt19937_64 mt64(init_seed);
						 std::poisson_distribution<Value, double> pd(m);
						 for (auto& elems : values)
						     elems = pd.operator()(mt64);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
						{
						 std::ranlux24 rl24(init_seed);
					         std::poisson_distribution<Value, double> pd(m);
						 for (auto& elems : values)
						     elems = pd.operator()(rl24);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
						{
					          std::ranlux48 rl48(init_seed);
						  std::poisson_distribution<Value, double> pd(m);
						  for (auto& elems : values)
				 		       elems = pd.operator()(rl48);
						}
						default:
						{
								   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
						}

						}

					}


					template<typename Value> template<typename T = Value> inline std::enable_if<std::is_floating_point<T>::value, T> PoissonDistributionImpl<Value>::exponential_distribution(std::vector<Value>& values,
						const Value lambda, unsigned long long init_seed, RandomNumberEngines engines)
					{

						unsigned long long seed = 0ULL;
						seed = ::clock();
						init_seed = init_seed + seed;
						switch (engines)
						{
						case PoissonDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
						{
						   std::default_random_engine dre(init_seed);
						   std::exponential_distribution<Value> ed(lambda);
				                   for (auto& elems : values)
                                                        elems = ed.operator()(dre);
																										   ////
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
						{
						  const unsigned int val1 = 56U;
                                                  const unsigned int val2 = 65536U;
                                                  const unsigned int val3 = 2147483647U;
                                                  std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
                                                  std::exponential_distribution<Value> ed(lambda);
                                                  for (auto& elems : values)
                                                  elems = ed.operator()(lce);																													
																												//
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
						{
						   std::minstd_rand msr(init_seed);
						   std::exponential_distribution<Value> ed(lambda);
						   for (auto& elems : values)
						       elems = ed.operator()(msr);
						   //std::cout << "std::minstd_rand msr(%lu) - called" << init_seed << std::endl;

						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
						{
						   std::minstd_rand0 msr0(init_seed);
						   std::exponential_distribution<Value> ed(lambda);
						   for (auto& elems : values)
							elems = ed.operator()(msr0);

						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MT199937:
						{
						  std::mt19937 mt(init_seed);
						  std::exponential_distribution<Value> ed(lambda);
						  for (auto& elems : values)
							  elems = ed.operator()(mt);

						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
						{
						 std::mt19937_64 mt64(init_seed);
						 std::exponential_distribution<Value> ed(lambda);
						 for (auto& elems : values)
							 elems = ed.operator()(mt64);

						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
						{
  					          std::ranlux24 rl24(init_seed);
						  std::exponential_distribution<Value> ed(lambda);
					          for (auto& elems : values)
							   elems = ed.operator()(rl24);

						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
						{
						   std::ranlux48 rl48(init_seed);
						   std::exponential_distribution<Value> ed(lambda);
						   for (auto& elems : values)
							   elems = ed.operator()(rl48);

						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
						{
						 std::knuth_b kbe(init_seed);
						 std::exponential_distribution<Value> ed(lambda);
						 for (auto& elems : values)
							 elems = ed.operator()(kbe);

						}
							break;
						default:
						{
						     throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
						}

						} 

					}

					template<typename Value> template<typename T = Value> 
					inline std::enable_if<std::is_floating_point<T>::value, T>   PoissonDistributionImpl<Value>::gamma_distribution(std::vector<Value>& values,
					const Value alpha, const Value beta, unsigned long long init_seed, RandomNumberEngines engines)
					{

						unsigned long long seed = 0ULL;
						init_seed = init_seed ^ (seed = clock());
						switch (engines)
						{
						case PoissonDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
						{
						   std::default_random_engine dre(init_seed);
						   std::gamma_distribution<Value> gd(alpha, beta);
						   for (auto& elems : values)
							   elems = gd.operator()(dre);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
						{
						  const unsigned int val1 = 56U;
						  const unsigned int val2 = 65536U;
						  const unsigned int val3 = 2147483647U;
						  std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
						  std::gamma_distribution<Value> gd(alpha, beta);
						  for (auto& elems : values)
							  elems = gd.operator()(lce);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
						{
						   std::minstd_rand msr(init_seed);
						   std::gamma_distribution<Value> gd(alpha, beta);
						   for (auto& elems : values)
							   elems = gd.operator()(msr);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
						{
						  std::minstd_rand0 msr0(init_seed);
						  std::gamma_distribution<Value> gd(alpha, beta);
						  for (auto& elems : values)
							  elems = gd.operator()(msr0);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MT199937:
						{
						  std::mt19937 mt(init_seed);
						  std::gamma_distribution<Value> gd(alpha, beta);
						  for (auto& elems : values)
							  elems = gd.operator()(mt);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
						{
						 std::mt19937_64 mt64(init_seed);
  					         std::gamma_distribution<Value> gd(alpha, beta);
						 for (auto& elems : values)
							 elems = gd.operator()(mt64);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
						{
						 std::ranlux24 rl24(init_seed);
						 std::gamma_distribution<Value> gd(alpha, beta);
						 for (auto& elems : values)
							 elems = gd.operator()(rl24);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
						{
					         std::ranlux48 rl48(init_seed);
					         std::gamma_distribution<Value> gd(alpha, beta);
						 for (auto& elems : values)
							   elems = gd.operator()(rl48);
						}
							break;
						case PoissonDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
						{
						 std::knuth_b kbe(init_seed);
						 std::gamma_distribution<Value> gd(alpha, beta);
						 for (auto& elems : values)
							 elems = gd.operator()(kbe);
						}
							break;
						default:
						{
						    throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
						}

						}
					}

						template<typename Value> template<typename T = Value> 
						inline std::enable_if<std::is_floating_point<T>::value, T> PoissonDistributionImpl<Value>::weibull_distribution(std::vector<Value>& values,
					        const Value a, const Value b, unsigned long long init_seed, RandomNumberEngines engines)
						{

							unsigned long long seed = 0ULL;
							init_seed = init_seed ^ (seed = clock());
							switch (engines)
							{
							case PoissonDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
							{
							   std::default_random_engine dre(init_seed);
							   std::weibull_distribution<Value> wb(a, b);
							   for (auto& elems : values)
								   elems = wb.operator()(dre);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
							{
							  const unsigned int val1 = 56U;
							  const unsigned int val2 = 65536U;
							  const unsigned int val3 = 2147483647U;
							  std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
							  std::weibull_distribution<Value> wd(a, b);
							  for (auto& elems : values)
								elems = wd.operator()(lce);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
							{
							   std::minstd_rand msr(init_seed);
							   std::weibull_distribution<Value> wd(a, b);
							   for (auto& elems : values)
								   elems = wd.operator()(msr);
							
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
							{
							  std::minstd_rand0 msr0(init_seed);
							  std::weibull_distribution<Value> wd(a, b);
							  for (auto& elems : values)
								elems = wd.operator()(msr0);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::MT199937:
							{
							  std::mt19937 mt(init_seed);
							  std::weibull_distribution<Value> wd(a, b);
						          for (auto& elems : values)
								  elems = wd.operator()(mt);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
							{
							 std::mt19937_64 mt64(init_seed);
							 std::weibull_distribution<Value> wd(a, b);
							 for (auto& elems : values)
								 elems = wd.operator()(mt64);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
							{
							   std::ranlux24 rl24(init_seed);
							   std::weibull_distribution<Value> wd(a, b);
						           for (auto& elems : values)
								   elems = wd.operator()(rl24);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
							{
							   std::ranlux48 rl48(init_seed);
							   std::weibull_distribution<Value> wd(a, b);
							   for (auto& elems : values)
								   elems = wd.operator()(rl48);
																								   
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
							{
							 std::knuth_b kbe(init_seed);
							 std::weibull_distribution<Value> wd(a, b);
							 for (auto& elems : values)
								 elems = wd.operator()(kbe);

							}
								break;
							default:
							{
									   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
							}

							}
						}

                        
						template<typename Value> template<typename T = Value> 
						inline std::enable_if<std::is_floating_point<T>::value, T> 
						PoissonDistributionImpl<Value>::extreme_value_distribution(std::vector<Value>& values,
							const Value a, const Value b, unsigned long long init_seed, RandomNumberEngines engines)
						{

							unsigned long long seed = 0ULL;
							init_seed = init_seed ^ (seed = clock());
							switch (engines)
							{
							case PoissonDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
							{
							   std::default_random_engine dre(init_seed);
							   std::extreme_value_distribution<Value> evd(a, b);
							   for (auto& elems : values)
								   elems = evd.operator()(dre);

							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
							{
							const unsigned int val1 = 56U;
							const unsigned int val2 = 65536U;
							const unsigned int val3 = 2147483647U;
							std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
							std::extreme_value_distribution<Value> evd(a, b);
							for (auto& elems : values)
								elems = evd.operator()(lce);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
							{
						           std::minstd_rand msr(init_seed);
							   std::extreme_value_distribution<Value> evd(a, b);
							   for (auto& elems : values)
								   elems = evd.operator()(msr);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
							{
							  std::minstd_rand0 msr0(init_seed);
							  std::extreme_value_distribution<Value> evd(a, b);
							  for (auto& elems : values)
								elems = evd.operator()(msr0);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::MT199937:
							{
						          std::mt19937 mt(init_seed);
							  std::extreme_value_distribution<Value> evd(a, b);
							  for (auto& elems : values)
  							         elems = evd.operator()(mt);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
							{
							 std::mt19937_64 mt64(init_seed);
							 std::extreme_value_distribution<Value> evd(a, b);
							 for (auto& elems : values)
								 elems = evd.operator()(mt64);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
							{
							   std::ranlux24 rl24(init_seed);
							   std::extreme_value_distribution<Value> evd(a, b);
							   for (auto& elems : values)
								   elems = evd.operator()(rl24);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
							{
							   std::ranlux48 rl48(init_seed);
							   std::extreme_value_distribution<Value> evd(a, b);
							   for (auto& elems : values)
								   elems = evd.operator()(rl48);
							}
								break;
							case PoissonDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
							{
							 std::knuth_b kbe(init_seed);
							 std::extreme_value_distribution<Value> evd(a, b);
							 for (auto& elems : values)
								 elems = evd.operator()(kbe);
							}
								break;
							default:
							{
									   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
							}

							}

						}

						// class NormalDistributionImpl
						template<typename Value> template<typename T = Value> 
						inline std::enable_if<std::is_floating_point<T>::value, T> 
						NormalDistributionImpl<Value>::normal_distribution(std::vector<Value>& values,
							const Value mikron, const Value sigma, unsigned long long init_seed, RandomNumberEngines engines)
						{

							unsigned long long seed = 0ULL;
							init_seed = init_seed ^ (seed = clock());
							switch (engines)
							{
							case NormalDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
							{
							  std::default_random_engine dre(init_seed);
							  std::normal_distribution<Value> nd(mikron, sigma);
   						          for (auto& elems : values)
								  elems = nd.operator()(dre);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
							{

							   const unsigned int val1 = 56U;
							   const unsigned int val2 = 65536U;
							   const unsigned int val3 = 2147483647U;
						           std::linear_congruential_engine<unsigned int,val1,val2,val3> lce(init_seed);
							   std::normal_distribution<Value> nd(mikron, sigma);
							   for (auto& elems : values)
								   elems = nd.operator()(lce);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
							{
							  std::minstd_rand msr(init_seed);
							  std::normal_distribution<Value> nd(mikron, sigma);
							  for (auto& elems : values)
							        elems = nd.operator()(msr);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
							{
							   std::minstd_rand0 msr0(init_seed);
							   std::normal_distribution<Value> nd(mikron, sigma);
							   for (auto& elems : values)
								   elems = nd.operator()(msr0);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937:
							{
							 std::mt19937 mt(init_seed);
							 std::normal_distribution<Value> nd(mikron, sigma);
							 for (auto& elems : values)
								 elems = nd.operator()(mt);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
							{
							std::mt19937_64 mt64(init_seed);
							std::normal_distribution<Value> nd(mikron, sigma);
							for (auto& elems : values)
								elems = nd.operator()(mt64);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
							{
							  std::ranlux24 rl24(init_seed);
							  std::normal_distribution<Value> nd(mikron, sigma);
							  for (auto& elems : values)
								  elems = nd.operator()(rl24);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
							{
							  std::ranlux48 rl48(init_seed);
							  std::normal_distribution<Value> nd(mikron, sigma);
							  for (auto& elems : values)
								  elems = nd.operator()(rl48);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
							{
							std::knuth_b kbe(init_seed);
							std::normal_distribution<Value> nd(mikron, sigma);
							for (auto& elems : values)
								elems = nd.operator()(kbe);
							}
								break;
							default:
							{
									   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
							}
							}
						}

						template<typename Value> template<typename T = Value> 
						inline std::enable_if<std::is_floating_point<T>::value, T> 
						NormalDistributionImpl<Value>::lognormal_ditribution(std::vector<Value>& values,
							const Value mikron, const Value sigma, unsigned long long init_seed, RandomNumberEngines engines)
						{

							unsigned long long seed = 0ULL;
							init_seed = init_seed ^ (seed = clock());
							switch (engines)
							{
							case NormalDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
							{
							  std::default_random_engine dre(init_seed);
							  std::lognormal_distribution<Value> ld(mikron, sigma);
							  for (auto& elems : values)
								  elems = ld.operator()(dre);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
							{
							   const unsigned int val1 = 56U;
							   const unsigned int val2 = 65536U;
							   const unsigned int val3 = 2147483647U;
							   std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
							   std::lognormal_distribution<Value> ld(mikron, sigma);
							   for (auto& elems : values)
								   elems = ld.operator()(lce);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
							{
							  std::minstd_rand msr(init_seed);
							  std::lognormal_distribution<Value> ld(mikron, sigma);
   						          for (auto& elems : values)
								  elems = ld.operator()(msr);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
							{
							   std::minstd_rand0 msr0(init_seed);
							   std::lognormal_distribution<Value> ld(mikron, sigma);
							   for (auto& elems : values)
								   elems = ld.operator()(msr0);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937:
							{
							 std::mt19937 mt(init_seed);
							 std::lognormal_distribution<Value> ld(mikron, sigma);
							 for (auto& elems : values)
								 elems = ld.operator()(mt);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
							{
							 std::mt19937_64 mt64(init_seed);
							 std::lognormal_distribution<Value> ld(mikron, sigma);
							 for (auto& elems : values)
							 elems = ld.operator()(mt64);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
							{
							  std::ranlux24 rl24(init_seed);
							  std::lognormal_distribution<Value> ld(mikron, sigma);
							  for (auto& elems : values)
								  elems = ld.operator()(rl24);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
							{
							  std::ranlux48 rl48(init_seed);
							  std::lognormal_distribution<Value> ld(mikron, sigma);
							  for (auto& elems : values)
								  elems = ld.operator()(rl48);
							}
								break;
							case NormalDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
							{
							std::knuth_b kbe(init_seed);
							std::lognormal_distribution<Value> ld(mikron, sigma);
							for (auto& elems : values)
								elems = ld.operator()(kbe);
							}
								break;
							default:
							{
									   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
							}

							}
						}

						//chi_squared_distribution()

                                            template<typename Value> template<typename T = Value> 
                                            inline std::enable_if<std::is_floating_point<T>::value, T> 
                                            NormalDistributionImpl<Value>::chi_squared_distribution(std::vector<Value>& values,
	                                    const Value n, unsigned long long init_seed, RandomNumberEngines engines)
                                            {


						unsigned long long seed = 0ULL;
						init_seed = init_seed ^ (seed = clock());
						switch (engines)
						{
						case NormalDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
						{
						  std::default_random_engine dre(init_seed);
						  std::chi_squared_distribution<Value> csd(n);
						  for (auto& elems : values)
							  elems = csd.operator()(dre);
						}
							break;
						case NormalDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
						{
						   const unsigned int val1 = 56U;
						   const unsigned int val2 = 65536U;
						   const unsigned int val3 = 2147483647U;
						   std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
						   std::chi_squared_distribution<Value> csd(n);
						   for (auto& elems : values)
							   elems = csd.operator()(lce);
						}
							break;
						case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
						{
						  std::minstd_rand msr(init_seed);
						  std::chi_squared_distribution<Value> csd(n);
						  for (auto& elems : values)
							  elems = csd.operator()(msr);
						}
							break;
						case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
						{
						   std::minstd_rand0 msr0(init_seed);
						   std::chi_squared_distribution<Value> csd(n);
						   for (auto& elems : values)
							   elems = csd.operator()(msr0);
						}
							break;
						case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937:
						{
						 std::mt19937 mt(init_seed);
						 std::chi_squared_distribution<Value> csd(n);
						 for (auto& elems : values)
							 elems = csd.operator()(mt);
						}
							break;
						case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
						{
						std::mt19937_64 mt64(init_seed);
						std::chi_squared_distribution<Value> csd(n);
						for (auto& elems : values)
							elems = csd.operator()(mt64);
						}
							break;
						case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
						{
						  std::ranlux24 rl24(init_seed);
						  std::chi_squared_distribution<Value> csd(n);
						  for (auto& elems : values)
							  elems = csd.operator()(rl24);
						}
							break;
						case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
						{
						  std::ranlux48 rl48(init_seed);
						  std::chi_squared_distribution<Value> csd(n);
						  for (auto& elems : values)
							  elems = csd.operator()(rl48);
						}
							break;
						case NormalDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
						{
						std::knuth_b kbe(init_seed);
						std::chi_squared_distribution<Value> csd(n);
						for (auto& elems : values)
							elems = csd.operator()(kbe);
						}
							break;
						default:
						{
								   throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
						}

						}
                                        }

					//cauchy_distribution

                          template<typename Value> template<typename T = Value> 
                          inline std::enable_if<std::is_floating_point<T>::value, T> 
                          NormalDistributionImpl<Value>::cauchy_distribution(std::vector<Value>& values,
	                   const Value a, const Value b, unsigned long long init_seed, RandomNumberEngines engines)
                          {

	         
	                                    unsigned long long seed = 0ULL;
										init_seed = init_seed ^ (seed = clock());
										switch (engines)
										{
										case NormalDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
										{
										  std::default_random_engine dre(init_seed);
										  std::cauchy_distribution<Value> cd(a, b);
										  for (auto& elems : values)
											  elems = cd.operator()(dre);
										}
											break;
										case NormalDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
										{
										   const unsigned int val1 = 56U;
										   const unsigned int val2 = 65536U;
										   const unsigned int val3 = 2147483647U;
									           std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
										   std::cauchy_distribution<Value> cd(a, b);
										   for (auto& elems : values)
											   elems = cd.operator()(lce);
										}
											break;
										case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
										{
										  std::minstd_rand msr(init_seed);
										  std::cauchy_distribution<Value> cd(a, b);
  									          for (auto& elems : values)
											  elems = cd.operator()(msr);
										}
											break;
										case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
										{
										   std::minstd_rand0 msr0(init_seed);
										   std::cauchy_distribution<Value> cd(a, b);
										   for (auto& elems : values)
											   elems = cd.operator()(msr0);
										}
											break;
										case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937:
										{
											 std::mt19937 mt(init_seed);
											 std::cauchy_distribution<Value> cd(a, b);
											 for (auto& elems : values)
												 elems = cd.operator()(mt);
										}
											break;
										case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
										{
										std::mt19937_64 mt64(init_seed);
										std::cauchy_distribution<Value> cd(a, b);
										for (auto& elems : values)
 										     elems = cd.operator()(mt64);
										}
											break;
										case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
										{
										  std::ranlux24 rl24(init_seed);
										  std::cauchy_distribution<Value> cd(a, b);
										  for (auto& elems : values)
											  elems = cd.operator()(rl24);
										}
											break;
										case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
										{
										  std::ranlux48 rl48(init_seed);
										  std::cauchy_distribution<Value> cd(a, b);
										  for (auto& elems : values)
											  elems = cd.operator()(rl48);
										}
											break;
										case NormalDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
										{
										std::knuth_b kbe(init_seed);
										std::cauchy_distribution<Value> cd(a, b);
										for (auto& elems : values)
											elems = cd.operator()(kbe);
										}
											break;
										default:
										{
									  throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
										}

										}

}

                                        template<typename Value> template<typename T = Value> 
                                        inline std::enable_if<std::is_floating_point<T>::value, T> NormalDistributionImpl<Value>::fisher_f_distribution(std::vector<Value>& values,
	                                 const Value m, const Value n, unsigned long long init_seed, RandomNumberEngines engines)
                                       {

											unsigned long long seed = 0ULL;
											init_seed = init_seed ^ (seed = clock());
											switch (engines)
											{
											case NormalDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
											{
											  std::default_random_engine dre(init_seed);
											  std::fisher_f_distribution<Value> fd(m, n);
											  for (auto& elems : values)
												  elems = fd.operator()(dre);
											}
												break;
											case NormalDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
											{
											   const unsigned int val1 = 56U;
											   const unsigned int val2 = 65536U;
											   const unsigned int val3 = 2147483647U;
											   std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
											   std::fisher_f_distribution<Value> fd(m, n);
											   for (auto& elems : values)
											   elems = fd.operator()(lce);
											}
												break;
											case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
											{
											  std::minstd_rand msr(init_seed);
											  std::fisher_f_distribution<Value> fd(m, n);
											  for (auto& elems : values)
												  elems = fd.operator()(msr);
											}
												break;
											case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
											{
											   std::minstd_rand0 msr0(init_seed);
										           std::fisher_f_distribution<Value> fd(m, n);
											   for (auto& elems : values)
												   elems = fd.operator()(msr0);
											}
												break;
											case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937:
											{
											 std::mt19937 mt(init_seed);
											 std::fisher_f_distribution<Value> fd(m, n);
											 for (auto& elems : values)
												 elems = fd.operator()(mt);
											}
												break;
											case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
											{
											std::mt19937_64 mt64(init_seed);
											std::fisher_f_distribution<Value> fd(m, n);
											for (auto& elems : values)
												elems = fd.operator()(mt64);
											}
												break;
											case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
											{
											  std::ranlux24 rl24(init_seed);
											  std::fisher_f_distribution<Value> fd(m, n);
											  for (auto& elems : values)
												  elems = fd.operator()(rl24);
											}
												break;
											case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
											{
											  std::ranlux48 rl48(init_seed);
											  std::fisher_f_distribution<Value> fd(m, n);
											  for (auto& elems : values)
												  elems = fd.operator()(rl48);
											}
												break;
											case NormalDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
											{
											std::knuth_b kbe(init_seed);
											std::fisher_f_distribution<Value> fd(m, n);
											for (auto& elems : values)
												elems = fd.operator()(kbe);
											}
												break;
											default:
											{
					throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
											}

											}
}

										//student_t_distribution()

                                         template<typename Value> template<typename T = Value> 
                                         inline std::enable_if<std::is_floating_point<T>::value, T> NormalDistributionImpl<Value>::student_t_distribution(std::vector<Value>& values,
	                                  const Value m, unsigned long long init_seed, RandomNumberEngines engines)
{

											 unsigned long long seed = 0ULL;
											 init_seed = init_seed ^ (seed = clock());
											 switch (engines)
											 {
											 case NormalDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
											 {
											   std::default_random_engine dre(init_seed);
											   std::student_t_distribution<Value> sd(m);
											   for (auto& elems : values)
												   elems = sd.operator()(dre);
											 }
												 break;
											 case NormalDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
											 {
												const unsigned int val1 = 56U;
												const unsigned int val2 = 65536U;
												const unsigned int val3 = 2147483647U;
												std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
												std::student_t_distribution<Value> sd(m);
												for (auto& elems : values)
													elems = sd.operator()(lce);
											 }
												 break;
											 case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
											 {
										           std::minstd_rand msr(init_seed);
											   std::student_t_distribution<Value> sd(m);
											   for (auto& elems : values)
											        elems = sd.operator()(msr);
											 }
												 break;
											 case NormalDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
											 {
											std::minstd_rand0 msr0(init_seed);
											std::student_t_distribution<Value> sd(m);
											for (auto& elems : values)
												elems = sd.operator()(msr0);
											 }
												 break;
											 case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937:
											 {
											  std::mt19937 mt(init_seed);
											  std::student_t_distribution<Value> sd(m);
											  for (auto& elems : values)
												  elems = sd.operator()(mt);
											 }
												 break;
											 case NormalDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
											 {
											 std::mt19937_64 mt64(init_seed);
											 std::student_t_distribution<Value> sd(m);
											 for (auto& elems : values)
												 elems = sd.operator()(mt64);
											 }
												 break;
											 case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
											 {
											   std::ranlux24 rl24(init_seed);
											   std::student_t_distribution<Value> sd(m);
											   for (auto& elems : values)
												   elems = sd.operator()(rl24);
											 }
												 break;
											 case NormalDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
											 {
										           std::ranlux48 rl48(init_seed);
											   std::student_t_distribution<Value> sd(m);
											   for (auto& elems : values)
												   elems = sd.operator()(rl48);
											 }
												 break;
											 case NormalDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
											 {
											 std::knuth_b kbe(init_seed);
											 std::student_t_distribution<Value> sd(m);
										         for (auto& elems : values)
												 elems = sd.operator()(kbe);
											 }
												 break;
											 default:
											 {
							throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
											 }

											 }

}

template<typename Value> template<typename T = Value> inline std::enable_if<std::is_integral<T>::value, T> 
SamplingDistributionImpl<Value>::discrete_distribution(std::vector<Value>& values,
				std::initializer_list<double> ilist, unsigned long long init_seed, RandomNumberEngines engines)
										 {
#if defined (_DEBUG)
											 for (auto& elems : ilist)
												 _ASSERTE(fabs(elems) > 0.0);
											
#else
											 for (auto p = ilist.begin(); p != ilist.end(); p++)
											 {
												 if (fabs(*p) <= 0.0)
												 {
	 throw std::runtime_error(std::string("ERROR in SamplingDistributionImpl::discrete_distribution() - Invalid Argument Value: Out of Range\n"));
												 }
											 }

#endif
											 unsigned long long seed = 0ULL;
											 init_seed = init_seed ^ (seed = clock());
											 switch (engines)
											 {
											 case SamplingDistributionImpl<Value>::RandomNumberEngines::DEFAULT_RANDOM_ENGINE:
											 {
											 std::default_random_engine dre(init_seed);
											 std::discrete_distribution<Value> dd(ilist);
											 for (auto& elems : values)
												 elems = dd.operator()(dre);
											 }
												 break;
											 case SamplingDistributionImpl<Value>::RandomNumberEngines::LINEAR_CONGRUENTIAL_ENGINE:
											 {
											  const unsigned int val1 = 56U;
											  const unsigned int val2 = 65536U;
											  const unsigned int val3 = 2147483647U;
											  std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
											  std::discrete_distribution<Value> dd(ilist);
											  for (auto& elems : values)
												  elems = dd.operator()(lce);
											 }
												 break;
											 case SamplingDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM:
											 {
										          std::minstd_rand mrs(init_seed);
											 std::discrete_distribution<Value> dd(ilist);
											 for (auto& elems : values)
												 elems = dd.operator()(mrs);
																														
											 }
												 break;
											 case SamplingDistributionImpl<Value>::RandomNumberEngines::MINSTD_RANDOM0:
											 {
											  std::minstd_rand0 mrs0(init_seed);
											  std::discrete_distribution<Value> dd(ilist);
										          for (auto& elems : values)
												  elems = dd.operator()(mrs0);
											 }
												 break;
											 case SamplingDistributionImpl<Value>::RandomNumberEngines::MT199937:
											 {
											std::mt19937 mt(init_seed);
											std::discrete_distribution<Value> dd(ilist);
											for (auto& elems : values)
												elems = dd.operator()(mt);
											 }
												 break;
											 case SamplingDistributionImpl<Value>::RandomNumberEngines::MT199937_64:
											 {
											   std::mt19937_64 mt64(init_seed);
											   std::discrete_distribution<Value> dd(ilist);
											   for (auto& elems : values)
												   elems = dd.operator()(mt64);
											 }
												 break;
											 case SamplingDistributionImpl<Value>::RandomNumberEngines::RAN_LUX24:
											 {
											 std::ranlux24 rl24(init_seed);
											 std::discrete_distribution<Value> dd(ilist);
											 for (auto& elems : values)
												 elems = dd.operator()(rl24);
											 }
												 break;
											 case SamplingDistributionImpl<Value>::RandomNumberEngines::RAN_LUX48:
											 {
											 std::ranlux48 rl48(init_seed);
											 std::discrete_distribution<Value> dd(ilist);
											 for (auto& elems : values)
												 elems = dd.operator()(rl48);
											 }
												 break;
											 case SamplingDistributionImpl<Value>::RandomNumberEngines::KNUTH_B:
											 {
											   std::knuth_b kbe(init_seed);
											   std::discrete_distribution<Value> dd(ilist);
											   for (auto& elems : values)
												   elems = dd.operator()(kbe);
																												  
											 }
												 break;
											 default:
											 {
									throw std::runtime_error(std::string("Invalid switch argument has been passed, only values [0:8] are accepted\n"));
											 }
											 }
										 }


                         //For Testing enum class implementation.
				     template<typename Value>    inline  void BernoulliDistributionImpl<Value>::test22(RanNumGenerators generators)
					{
						switch (generators)
						{
						case BernoulliDistributionImpl<Value>::RanNumGenerators::DEFAULT_RANDOM_ENGINE:
						{
						  std::cout << "case DEFAULT_RANDOM_ENGINE\n" << std::endl;
						}
							break;
						case BernoulliDistributionImpl<Value>::RanNumGenerators::KNUTH_B:
						{
							std::cout << "case KNUTH_B\n" << std::endl;
						}
							break;
						default:
						{
							 std::cout << "No valid choice has been made\n" << std::endl;
						}

						}// End of switch statement
						
					}
					 //End of test22 function.

                       // To be called by C wrapper.
					/*void gms::test_binomial_distribution_wrapper(unsigned int * ptr_data, const int length_data)
					{
						std::vector<unsigned int> vec(length_data);
						unsigned int t = 5;
						const double p = 0.5;
						size_t seed = clock();
			                        std::enable_if<std::is_integral<unsigned int>::value, unsigned int> s1 = 
			                        gms::BernoulliDistributionImpl<unsigned int>::binomial_distribution(vec, t, p, seed, 2);
						for (auto i = 0; i != length_data; ++i)
							ptr_data[i] = vec[i];
					}*/
					//
					//Implementation  of BernoulliDistributionsWrapper class
					
			template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::binomial_distribution_wrapper(Value * data, const int data_length, const int num_of_engine)
					{
#if defined (_DEBUG)
						_ASSERTE(num_of_engine >= 0 && num_of_engine <= 8);
						_ASSERTE((data != nullptr) || (data_length > 0 && data_length < MaxUintVal));
						
#else
						if((data == nullptr) || (data_length < 0 || data_length > std::numeric_limits<unsigned int>::max()) || (num_of_engine < 0 || num_of_engine > 8))
						{
							return false;
						}
#endif
						std::vector<Value> vec(data_length);
						Value t = 5;
						const double p = 0.5;
						size_t seed = clock() ;
						auto s1 = gms::BernoulliDistributionImpl<Value>::binomial_distribution(vec, t, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						
						return true;
					}


					//Helper functions , private.
					template<typename Value> inline void BernoulliDistributionsWrapper<Value>::test_assertion(Value * data, const int data_length,
						const int num_of_engine)
					{
						_ASSERTE(num_of_engine >= 0 && num_of_engine <= 8);
						_ASSERTE((data != nullptr) || (data_length > 0 && data_length < MaxUintVal));
					}


					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::check_arguments(Value * data, const int data_length,
						const int num_of_engine)
					{
						if ((data == nullptr) || (data_length < 0 || data_length > MaxUintVal) || (num_of_engine < 0 || num_of_engine > 8))
						{
							return false;
						}

						return true;
					}


template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::binomial_distribution_wrapper(Value * data, const int data_length, Value t, const int num_of_engine)
					{
#if defined (_DEBUG)
						//_ASSERTE(num_of_engine >= 0 && num_of_engine <= 8);
						//_ASSERTE((data != nullptr) || (data_length > 0 && data_length < std::numeric_limits<unsigned int>::max()));
						test_assertion(data, data_length, num_of_engine);
#else
						if((data == nullptr) || (data_length < 0 || data_length > std::numeric_limits<unsigned int>::max()) || (num_of_engine < 0 || num_of_engine > 8))
						{
							return false;
						}
#endif
						std::vector<Value> vec(data_length);
						const double p = 0.5;
						size_t seed{ clock() };
						auto s1 = gms::BernoulliDistributionImpl<Value>::binomial_distribution(vec, t, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::binomial_distribution_wrapper(Value * data, const int data_length, Value t, const double p, const int num_of_engine)
					{
#if defined (_DEBUG)
						_ASSERTE(num_of_engine >= 0 && num_of_engine <= 8);
						_ASSERTE((data != nullptr) || (data_length > 0 && data_length < MaxUintVal));
#else
						if((data == nullptr) || (data_length < 0 || data_length > std::numeric_limits<unsigned int>::max()) || (num_of_engine < 0 || num_of_engine > 8))
						{
							return false;
						}
#endif
						std::vector<Value> vec(data_length);
						size_t seed{ ::clock() };
						auto s1 = gms::BernoulliDistributionImpl<Value>::binomial_distribution(vec, t, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}

						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::binomial_distribution_wrapper(Value * data, const int data_length,
					 Value t, const double p, std::size_t seed, const int num_of_engine	)
					{
#if defined (_DEBUG)
						
						_ASSERTE((data != nullptr) || (data_length > 0 && data_length < MaxUintVal));
#else
						if ((data == nullptr) || (data_length < 0 || data_length > std::numeric_limits<unsigned int>::max()))
						{
							return false;
						}
#endif
						std::vector<Value> vec(data_length);
						auto s1 = gms::BernoulliDistributionImpl<Value>::binomial_distribution(vec, t, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::bernoulli_distribution_wrapper(Value * data, const int data_length, const int num_of_engine)
					{
#if defined (_DEBUG)
						_ASSERTE(num_of_engine >= 0 && num_of_engine <= 8);
						_ASSERTE((data != nullptr) || (data_length > 0 && data_length < MaxUintVal));
#else
						if((data == nullptr) || (data_length < 0 || data_length > std::numeric_limits<unsigned int>::max()))
						{
							return false;
						}
#endif
						std::vector<Value> vec( data_length);
						const double p = 0.5;
						std::size_t seed{ clock() };
						auto s1 = gms::BernoulliDistributionImpl<Value>::bernoulli_distribution(vec, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::bernoulli_distribution_wrapper(Value * data, const int data_length,
						const double p, const int num_of_engine)
					{
#if defined (_DEBUG)
						
						_ASSERTE((data != nullptr) || (data_length > 0 && data_length < MaxUintVal));
#else
						if ((data == nullptr) || (data_length < 0 || data_length > std::numeric_limits<unsigned int>::max()) )
						{
							return false;
						}
#endif
						std::vector<Value> vec( data_length);
						std::size_t seed{ ::clock() };
						auto s1 = gms::BernoulliDistributionImpl<Value>::bernoulli_distribution(vec, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
							data[i] = vec[i];

						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::bernoulli_distribution_wrapper(Value * data, const int data_length,
						const double p, std::size_t seed, const int num_of_engine)
					{
#if defined (_DEBUG)
						
						_ASSERTE((data != nullptr) || (data_length > 0 && data_length < MaxUintVal));
#else
						if ((data == nullptr) || (data_length < 0 || data_length > std::numeric_limits<unsigned int>::max()))
						{
							return false;
						}
#endif
						std::vector<Value> vec( data_length);
						auto s1 = gms::BernoulliDistributionImpl<Value>::bernoulli_distribution(vec, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::geometric_distribution_wrapper(Value * data, const int data_length,
						const int num_of_engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length, num_of_engine);
#else
						if ((data == nullptr) || (data_length < 0 || data_length > std::numeric_limits<unsigned int>::max())) 
						{
							return false;
						}
#endif
						std::vector<Value> vec( data_length);
						const double p{ 0.5 };
						std::size_t seed{ ::clock() };
						auto s1 = gms::BernoulliDistributionImpl<Value>::geometric_distribution(vec, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::geometric_distribution_wrapper(Value * data, const int data_length,
						const double p, const int num_of_engine)
					{
#if defined (_DEBUG)

						test_assertion(data, data_length, num_of_engine);
#else
						if ((data == nullptr) || (data_length < 0 || data_length > std::numeric_limits<unsigned int>::max()))
						{
							return false;
						}
#endif
						std::vector<Value> vec( data_length);
#if defined USE_WIN_API
						std::size_t seed{ ::GetTickCount64() };
#else
						std::size_t seed{ ::clock() };
#endif
						auto s1 = gms::BernoulliDistributionImpl<Value>::geometric_distribution(vec, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::geometric_distribution_wrapper(Value * data, const int data_length,
						const double p, std::size_t seed, const int num_of_engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length, num_of_engine);
#else
						if ((data == nullptr) || (data_length < 0 || data_length > std::numeric_limits<unsigned int>::max()))
						{
							return false;
						}
#endif
						std::vector<Value> vec( data_length);
						auto s1 = gms::BernoulliDistributionImpl<Value>::geometric_distribution(vec, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::negative_binomial_distribution_wrapper(Value * data, const int data_length,
						const int num_of_engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length, num_of_engine);
#else
						if ((data == nullptr) || (data_length < 0 || data_length > MaxUintVal) )
						{
							return false;
						}
#endif
						std::vector<Value> vec( data_length);
#if defined USE_WIN_API
						std::size_t seed{ ::GetTickCount64() };
#else
						std::size_t seed{ ::clock() };
#endif
						const double p{ 0.5 };
						const Value t{ 5 };
						auto s1 = gms::BernoulliDistributionImpl<Value>::negative_binomial_distribution(vec, t, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::negative_binomial_distribution_wrapper(Value * data, const int data_length,
						const Value k, const int num_of_engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length, num_of_engine);
#else
						if ((data == nullptr) || (data_length < 0 || data_length > MaxUintVal))
						{
							return false;
						}
#endif
						std::vector<Value> vec( data_length);
#if defined (USE_WIN_API)
						std::size_t seed{ ::GetTickCount64() };
#else
						std::size_t seed{::clock()};
#endif
						const double p{ 0.5 };
						auto s1 = gms::BernoulliDistributionImpl<Value>::negative_binomial_distribution(vec, k, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::negative_binomial_distribution_wrapper(Value * data, const int data_length,
						const Value k, const double p, const int num_of_engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length, num_of_engine);
#else
						if ((data == nullptr) || (data_length < 0 || data_length > MaxUintVal))
						{
							return false;
						}
#endif
						std::vector<Value> vec( data_length);
#if defined (USE_WIN_API)
						std::size_t seed{ ::GetTickCount64() };
#else
						std::size_t seed{::clock()};
#endif
						auto s1 = gms::BernoulliDistributionImpl<Value>::negative_binomial_distribution(vec, k, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename Value> inline bool BernoulliDistributionsWrapper<Value>::negative_binomial_distribution_wrapper(Value * data, const int data_length,
						const Value k, const double p, std::size_t seed, const int num_of_engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length, num_of_engine);
#else
						if ((data == nullptr) || (data_length < 0 || data_length > MaxUintVal) )
						{
							return false;
						}
#endif
						std::vector<Value> vec{ data_length };
						auto s1 = gms::BernoulliDistributionImpl<Value>::negative_binomial_distribution(vec, k, p, seed, num_of_engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					// Implementation of PoissonDistributionWrapper class


					// Helper functions.

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::report_total_memory(unsigned long long phys_mem)
					{
						bool b_status = false;
						if ( b_status = ::GetPhysicallyInstalledSystemMemory(&phys_mem))
						{
#if defined (_DEBUG)
							std::printf("::GetPhysicallyInstalledSystemMemory() - Succedded, installed memory:[%2uKBs]\n\n", phys_mem);
#endif
							return true;
						}
						else
						{
#if defined (_DEBUG)
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with an error code %2d\n", ::GetLastError());
#endif
							return false;
						}
					}


					template<typename _Ty> inline void PoissonDistributionsWrapper<_Ty>::test_assertion(_Ty * data, const int data_length)
						
					{
#if defined (_DEBUG)
						unsigned long long available_mem = 0ULL;
						bool result = false;
						_ASSERTE(result = report_total_memory(available_mem));

						
						_ASSERTE((data != nullptr) || (data_length > 0 && data_length < (int)available_mem));
#endif

					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::poisson_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
						
#endif
						std::vector<_Ty> vec( data_length);
#if defined USE_WIN_API
					unsigned long long seed{ ::GetTickCount64() };
#else
					unsigned long long seed{::clock()};
#endif
						_Ty m{ 5 };
						auto s1 = gms::PoissonDistributionImpl<_Ty>::poisson_distribution(vec, m, seed,engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}

						return true;
						
					}


					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::poisson_distribution_wrapper(_Ty * data, const int data_length,
						_Ty m, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						
						auto s1 = gms::PoissonDistributionImpl<_Ty>::poisson_distribution(vec, m, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::poisson_distribution_wrapper(_Ty * data, const int data_length,
						_Ty m, unsigned long long seed, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
						auto s1 = gms::PoissonDistributionImpl<_Ty>::poisson_distribution(vec, m, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::exponential_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif

						std::vector<_Ty> vec( data_length);

					const	_Ty lambda{ 2.0 };


#if defined (USE_WIN_API)
						unsigned long long seed(  ::GetTickCount64()) ;
#else
						unsigned long long seed{::clock()};
#endif


						auto s1 = gms::PoissonDistributionImpl<_Ty>::exponential_distribution(vec, lambda, seed, engine);

						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}

						return true;
					}


					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::exponential_distribution_wrapper(_Ty * data, const int data_length,
						_Ty lambda, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem) )
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						auto s1 = gms::PoissonDistributionImpl<_Ty>::exponential_distribution(vec, lambda, seed, engines);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::exponential_distribution_wrapper(_Ty * data, const int data_length,
						_Ty lambda, unsigned long long seed, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem) )
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
						auto s1 = gms::PoissonDistributionImpl<_Ty>::exponential_distribution(vec, lambda, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::gamma_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						const _Ty alpha{ 2.0 };
						const _Ty beta{ 3.0 };
						auto s1 = gms::PoissonDistributionImpl<_Ty>::gamma_distribution(vec, alpha, beta, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::gamma_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty alpha, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem) )
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						const _Ty beta{ 3.0 };
						auto s1 = gms::PoissonDistributionImpl<_Ty>::gamma_distribution(vec, alpha, beta, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::gamma_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty alpha, const _Ty beta, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)

						
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem) )
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						auto s1 = gms::PoissonDistributionImpl<_Ty>::gamma_distribution(vec, alpha, beta, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::gamma_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty alpha, const _Ty beta, unsigned long long seed, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem) )
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
						auto s1 = gms::PoissonDistributionImpl<_Ty>::gamma_distribution(vec, alpha, beta, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::weibull_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem) )
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						const _Ty a = 2.0;
						const _Ty b = 3.0;
						auto s1 = gms::PoissonDistributionImpl<_Ty>::weibull_distribution(vec, a, b, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}

						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::weibull_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty a, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem) )
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						const _Ty b{ 3.0 };
						auto s1 = gms::PoissonDistributionImpl<_Ty>::weibull_distribution(vec, a, b, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::weibull_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty a, const _Ty b, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						
						test_assertion(data,data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem) )
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						auto s1 = gms::PoissonDistributionImpl<_Ty>::weibull_distribution(vec, a, b, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::weibull_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty a, const _Ty b, unsigned long long seed, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem) )
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
						auto s1 = gms::PoissonDistributionImpl<_Ty>::weibull_distribution(vec, a, b, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::extreme_value_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem) )
							{
								return false;
							}
						}
#endif
						// allocate vector for random data.
						std::vector<_Ty> vec( data_length);
						// Create seed.
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						// Initialize a and b variables
						const _Ty a{ 2.0 };
						const _Ty b{ 3.0 };
						auto s1 = gms::PoissonDistributionImpl<_Ty>::extreme_value_distribution(vec, a, b, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							// copy random values to client array.
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::extreme_value_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty a, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						const _Ty b{ 3.0 };
						auto s1 = gms::PoissonDistributionImpl<_Ty>::extreme_value_distribution(vec, a, b, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::extreme_value_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty a, const _Ty b, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						auto s1 = gms::PoissonDistributionImpl<_Ty>::extreme_value_distribution(vec, a, b, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool PoissonDistributionsWrapper<_Ty>::extreme_value_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty a, const _Ty b, unsigned long long seed, typename gms::PoissonDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
						auto s1 = gms::PoissonDistributionImpl<_Ty>::extreme_value_distribution(vec, a, b, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					// Implementation of NormalDistributionsWrapper class
					// Helper functions.

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::report_total_memory(unsigned long long phys_mem)
					{
						bool b_status = false;
						if (b_status = GetPhysicallyInstalledSystemMemory(&phys_mem))
						{
#if PRINT_DATA
							std::printf("::GetPhysicallyInstalledSystemMemory() - Succedded, installed memory:[%2u]", phys_mem);
#endif
							return true;
						}
						else
						{
#if PRINT_DATA
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with an error code %2d\n", GetLastError());
#endif
							return false;
						}
					}

					template<typename _Ty> inline void NormalDistributionsWrapper<_Ty>::test_assertion(_Ty * data, const int data_length)

					{
#if defined (_DEBUG)
						unsigned long long available_mem = 0ULL;
						bool result = false;
						_ASSERTE(result = report_total_memory(available_mem));


						_ASSERTE((data != nullptr) || (data_length > 0 && data_length < (int)available_mem));
#endif

					}


					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::normal_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						const _Ty m{ 2.0 };
						const _Ty s{ 5.0 };
						auto s1 = gms::NormalDistributionImpl<_Ty>::normal_distribution(vec, m, s, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::normal_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty m, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
						std::vector<_Ty> vec( data_length);
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{ ::clock() };
#endif
						const _Ty s{ 5.0 };
						auto s1 = gms::NormalDistributionImpl<_Ty>::normal_distribution(vec, m, s, seed, engine);
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::normal_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty m, const _Ty s, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						// check for nullptr and for data_length, on true evaluation call CrtDebugBreak()
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;

						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
						// Create vector<_Ty> vec
						std::vector<_Ty> vec{ data_length };
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
						// Call NormalDistributionImpl<_Ty>::normal_distribution(vec,m,s,seed,engine)
						auto s1 = gms::NormalDistributionImpl<_Ty>::normal_distribution(vec, m, s, seed, engine);
						//Copy resulting normal distribution random  values to user supplied _Ty * data.
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						// std::vector<_Ty>::~vector<_Ty> called here!!
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::normal_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty m, const _Ty s, unsigned long long seed, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						// check for nullptr and for data_length, on true evaluation call CrtDebugBreak()
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
						// Create vector<_Ty> vec
						std::vector<_Ty> vec{ data_length };
						// Call NormalDistributionImpl<_Ty>::normal_distribution(vec,m,s,seed,engine)
						auto s1 = gms::NormalDistributionImpl<_Ty>::normal_distribution(vec, m, s, seed, engine);
						//Copy resulting normal distribution random  values to user supplied _Ty * data.
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}

						// std::vector<_Ty>::~vector<_Ty> called here!!
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::lognormal_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						// check for nullptr and for data_length, on true  call CrtDebugBreak()
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec( data_length);

#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Initializing mikron and sigma variables, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty m{ 2.0 };
						const _Ty s{ 5.0 };
						// Call NormalDistributionImpl<_Ty>::lognormal_distribution(vec,m,s,seed,engine)
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::lognormal_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::lognormal_ditribution(vec, m, s, seed, engine);
						//Copy resulting normal distribution random  values to user supplied _Ty * data.
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						// std::vector<_Ty>::~vector<_Ty> called here!!
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::lognormal_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty m, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2d, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec( data_length);

#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2d, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif

#if defined (_DEBUG)
						std::printf("Initializing sigma variable, in Function: %2d, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty s{ 3.0 };
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::lognormal_distribution(), in Function: %2d, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::lognormal_ditribution(vec, m, s, seed, engine);

#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2d, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::lognormal_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty m, const _Ty s, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::lognormal_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::lognormal_ditribution(vec, m, s, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::lognormal_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty m, const _Ty s, unsigned long long seed, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec( data_length);
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::lognormal_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::lognormal_ditribution(vec, m, s, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::chi_squared_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Initializing sigma variable, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty n{ 3.0 };
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::chi_squared_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::chi_squared_distribution(vec, n, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::chi_squared_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty n, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::chi_squared_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::chi_squared_distribution(vec, n, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::chi_squared_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty n, unsigned long long seed, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);

#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::chi_squared_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::chi_squared_distribution(vec, n, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::cauchy_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Initializing a and b variables, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty a{ 2.0 };
						const _Ty b{ 5.0 };
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::cauchy_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::cauchy_distribution(vec, a, b, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::cauchy_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty a, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Initializing a and b variables, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						
						const _Ty b{ 45.5 };
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::cauchy_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::cauchy_distribution(vec, a, b, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::cauchy_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty a, const _Ty b, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::cauchy_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::cauchy_distribution(vec, a, b, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::cauchy_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty a, const _Ty b, unsigned long long seed, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::cauchy_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::cauchy_distribution(vec, a, b, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::fisher_f_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned  long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Initializing m and n variables, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty m{ 5.0 };
						const _Ty n{ 15.0 };
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::fisher_f_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::fisher_f_distribution(vec, m, n, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::fisher_f_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty m, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Initializing  n variable, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty n{ 10.5 };
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::fisher_f_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::fisher_f_distribution(vec, m, n, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::fisher_f_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty m, const _Ty n, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::fisher_f_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::fisher_f_distribution(vec, m, n, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::fisher_f_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty m, const _Ty n, unsigned long long seed, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::fisher_f_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::fisher_f_distribution(vec, m, n, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::student_t_distribution_wrapper(_Ty * data, const int data_length,
						typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if PRINT_DATA
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if PRINT_DATA
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if PRINT_DATA
						std::printf("Initializing  n variable, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty n{ 5.6 };
#if PRINT_DATA
						std::printf("Calling NormalDistributionImpl::student_t_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::student_t_distribution(vec, n, seed, engine);
#if PRINT_DATA
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::student_t_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty n, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						unsigned long long seed{ ::GetTickCount64() };
#else
						unsigned long long seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::student_t_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::student_t_distribution(vec, n, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}

					template<typename _Ty> inline bool NormalDistributionsWrapper<_Ty>::student_t_distribution_wrapper(_Ty * data, const int data_length,
						const _Ty n, unsigned long long seed, typename gms::NormalDistributionImpl<_Ty>::RandomNumberEngines engine)
					{
#if defined (_DEBUG)
						test_assertion(data, data_length);
#else
						unsigned long long available_mem = 0ULL;
						// runtime checks
						if (!report_total_memory(available_mem){
							std::printf("::GetPhysicallyInstalledSystemMemory() - Failed with na error %2d\n", ::GetLastError());
							return false;
						}
						else
						{     // check for nullptr or invalid data length when true, return false.
							if ((data == nullptr) || (data_length < 0 || data_length > available_mem))
							{
								return false;
							}
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Calling NormalDistributionImpl::student_t_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::NormalDistributionImpl<_Ty>::student_t_distribution(vec, n, seed, engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool UniformDistributionWrapper<_Ty>::uniform_int_distribution_wrapper(_Ty * data,
						const int data_length, const int num_of_engine)
					{
#if defined (_DEBUG)
						_ASSERTE((data != nullptr) || (data_length > 0) );
#else
						if((data == nullptr) || (data_length <= 0) )
						{
							return false;
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						std::size_t seed{ ::GetTickCount64() };
#else
						std::size_t seed{clock()};
#endif
#if defined (_DEBUG)
						std::printf("Initializing  a and b variables(fixed length interval), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty a{ 1 };
						const _Ty b{ 100 };
#if defined (_DEBUG)
						std::printf("Calling UniformDistributionImpl::uni_int_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::UniDistributions<_Ty>::uni_int_distribution(vec, a, b, seed, num_of_engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}

						return true;
					}

					template<typename _Ty> inline bool UniformDistributionWrapper<_Ty>::uniform_int_distribution_wrapper(_Ty * data,
						const int data_length, const _Ty a, const int num_of_engine)
					{
#if defined (_DEBUG)
						_ASSERTE((data != nullptr) || (data_length > 0) );
#else
						if((data == nullptr) || (data_length <= 0))
						{
							return false;
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						std::size_t seed{ ::GetTickCount64() };
#else
						std::size_t seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Initializing  b  variables(fixed length interval), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty b{ 100 };
#if defined (_DEBUG)
						std::printf("Calling UniformDistributionImpl::uni_int_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::UniDistributions<_Ty>::uni_int_distribution(vec, a, b, seed, num_of_engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool UniformDistributionWrapper<_Ty>::uniform_int_distribution_wrapper(_Ty * data,
						const int data_length, const _Ty a, const _Ty b, const int num_of_engine)
					{
#if defined (_DEBUG)
						_ASSERTE((data != nullptr) || (data_length > 0));
#else
						if ((data == nullptr) || (data_length <= 0))
						{
							return false;
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						std::size_t seed{ ::GetTickCount64() };
#else
						std::size_t seed{::clock()};
#endif
#if defined (_DEBUG)
						std::printf("Calling UniformDistributionImpl::uni_int_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::UniDistributions<_Ty>::uni_int_distribution(vec, a, b, seed, num_of_engine);
#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool UniformDistributionWrapper<_Ty>::uniform_int_distribution_wrapper(_Ty * data,
						const int data_length, const _Ty a, const _Ty b, std::size_t seed, const int num_of_engine)
					{
#if defined (_DEBUG)
						_ASSERTE((data != nullptr) || (data_length > 0));
#else
						if ((data == nullptr) || (data_length <= 0))
						{
							return false;
						}
#endif
#if defined (_DEBUG)
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if defined (_DEBUG)
						std::printf("Calling UniformDistributionImpl::uni_int_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::UniDistributions<_Ty>::uni_int_distribution(vec, a, b, seed, num_of_engine);

#if defined (_DEBUG)
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool UniformDistributionWrapper<_Ty>::uniform_real_distribution_wrapper(_Ty * data,
						const int data_length, const int num_of_engine)
					{
#if defined (_DEBUG)
						_ASSERTE((data != nullptr) || (data_length > 0));
#else
						if((data == nullptr) || (data_length <= 0))
						{
							return true;
						}
#endif
#if PRINT_DATA
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function at address: %2d, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if PRINT_DATA
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function at address: %2d, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						std::size_t seed{ ::GetTickCount64() };
#else
						std::size_t seed{::clock()};
#endif
#if PRINT_DATA
						std::printf("Initializing  a and b variables(fixed length interval), in Function at address: %2d, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty a{ 1.0 };
						const _Ty b{ 100.0 };
#if PRINT_DATA
						std::printf("Calling UniformDistributionImpl::uni_real_distribution(), in Function at address: %2d, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::UniDistributions<_Ty>::uni_real_distribution(vec, a, b, seed, num_of_engine);
#if PRINT_DATA
						std::printf("Copying random values to user supplied array type _Ty in Function at address: %2d, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}

						return true;
					}


					template<typename _Ty> inline bool UniformDistributionWrapper<_Ty>::uniform_real_distribution_wrapper(_Ty * data,
						const int data_length, const _Ty a, const int num_of_engine)
					{
#if (_DEBUG)
						_ASSERTE((data != nullptr) || (data_length > 0));
#else
						if ((data == nullptr) || (data_length <= 0))
						{
							return true;
						}
#endif
#if   PRINT_DATA
						std::printf("Instantiating std::vector<_Ty> vec{data_length} object, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						std::vector<_Ty> vec(data_length);
#if   PRINT_DATA
						std::printf("Initializing seed value: seed{::GetTickCount64()}, in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
#if defined (USE_WIN_API)
						std::size_t seed{ ::GetTickCount64() };
#else
						std::size_t seed{::clock()};
#endif
#if   PRINT_DATA
						std::printf("Initializing  b  variables(fixed length interval), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						const _Ty b{ 100.0 };
#if   PRINT_DATA
						std::printf("Calling UniformDistributionImpl::uni_real_distribution(), in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);
#endif
						auto s1 = gms::UniDistributions<_Ty>::uni_real_distribution(vec, a, b, seed, num_of_engine);
#if PRINT_DATA
						std::printf("Copying random values to user supplied array type _Ty in Function: %2s, at line of code: %2d\n", __FUNCTIONW__, __LINE__);

#endif
						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}


					template<typename _Ty> inline bool UniformDistributionWrapper<_Ty>::uniform_real_distribution_wrapper(_Ty * data,
						const int data_length, const _Ty a, const _Ty b ,const int num_of_engine)
					{



						std::vector<_Ty> vec(data_length);


						std::size_t seed{::clock()};
						auto s1 = gms::UniDistributions<_Ty>::uni_real_distribution(vec, a, b, seed, num_of_engine);

						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return;
					}


					template<typename _Ty> inline bool UniformDistributionWrapper<_Ty>::uniform_real_distribution_wrapper(_Ty * data,
						const int data_length, const _Ty a, const _Ty b, std::size_t seed, const int num_of_engine)
					{


						std::vector<_Ty> vec(data_length);

						auto s1 = gms::UniDistributions<_Ty>::uni_real_distribution(vec, a, b, seed, num_of_engine);

						for (auto i = 0; i != vec.size(); ++i)
						{
							data[i] = vec[i];
						}
						return true;
					}
