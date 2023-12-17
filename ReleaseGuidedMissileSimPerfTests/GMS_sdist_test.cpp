#include "stdafx.h"
#include "CppUnitTest.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace fdlibcppwrapper;

namespace StatDistWrappersTest
{		

	// TestUniformDistributionEquality(const int, Distribution, Engine, const T, const T)
	template<typename Distribution, typename Engine, typename T> void TestUniformDistributionEquality(const int vec_length,
		Distribution d, Engine e, const T min_val, const T max_val)
	{
		Logger::WriteMessage("Begin Test of uni_int_distribution() - Testing for Inequality \n");

		_ASSERTE(min_val < max_val);
		
		std::vector<T> actual_rand_values(vec_length);
		std::vector<T> expected_rand_values(vec_length);
		/*unsigned int init_seed = 0xFFFFFFFF;
		unsigned int status = 0xFFFFFFFF;
		if ((status = _rdrand32_step(&init_seed)) == 0)
		{
			
			throw std::runtime_error (std::string("RDRAND Failed \n"));
		}*/
		
		
		//std::knuth_b kbe(init_seed);
		//std::uniform_int_distribution<T> ud(min_val, max_val);
		Logger::WriteMessage("Calling uniform distribution functions\n");
		std::enable_if<std::is_integral<T>::value, T> s1 = UniformDistributions<Distribution, Engine, T>::uni_int_distribution(actual_rand_values,
			d,e);
		std::enable_if<std::is_integral<T>::value, T> s2 = UniformDistributions < Distribution, Engine, T>::uni_int_distribution(expected_rand_values,
			d,e);

		Logger::WriteMessage("Executing Assert::AreNotEqual statement\n");
		for (auto i = 0; i != vec_length; ++i)
		{
			Assert::AreNotEqual(expected_rand_values[i], actual_rand_values[i], L"Uniform Distribution Equality Test Failed", LINE_INFO());
		}

		Logger::WriteMessage("End of uni_int_distribution() Testing for Inequality \n");
	}

//--------------------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------------------//
	template<typename Distribution, typename Engine, typename Value> void TestUniformDistributionInequality(const int vec_length,
		Distribution d, Engine e, const Value min_val, const Value max_val)
	{
		Logger::WriteMessage("Begin Test of uni_int_distribution() - Testing for Equality \n");
		_ASSERTE(min_val < max_val);

		std::vector<T> actual_rand_values(vec_length);
		std::vector<T> expected_rand_values(vec_length);

		Logger::WriteMessage("Calling uniform distribution functions\n");
		std::enable_if<std::is_integral<Value>::value, Value> s1 = UniformDistributions<Distribution, Engine, Value>::uni_int_distribution(actual_rand_values,
			d, e);
		std::enable_if<std::is_integral<Value>::value, Value> s2 = UniformDistributions<Distributions, Engine, Value>::uni_int_distribution(expected_rand_values,
			d, e);
		Logger::WriteMessage("Executing Assert::AreEqual statement\n");
		for (auto i = 0; i != vec_length; ++i)
		{
			Assert::AreEqual(expected_rand_values[i], actual_rand_values[i], L"Uniform Distribution InEquality Failed", LINE_INFO());
		}
		Logger::WriteMessage("End of uni_int_distribution() Testing for Equality\n")
	}
//-------------------------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------------------------//

	// TestBernoulliDistributionEquality(const int, Distribution, Engine)
	template<typename Distribution, typename Engine, typename Value> void TestBernoulliDistributionEquality(const int vec_length,
		Distribution d, Engine e)
	{
		//_ASSERTE(fabs(probability) >= 0.0 && fabs(probability) <= 1.0);

		Logger::WriteMessage("Begin Test of bernoulli_distribution() - Testing for Equality\n");
		Logger::WriteMessage("Allocating data vectors\n");
		std::vector<Value> actual_rand_values(vec_length);
		std::vector<Value> expected_rand_values(vec_length);
		Logger::WriteMessage("Calling bernoulli distribution functions\n");
		std::enable_if<std::is_integral<Value>::value, Value> s1 = BernoulliDistributions<Distribution, Engine, Value>::bernoulli_distribution(actual_rand_values,
			d, e);
		std::enable_if<std::is_integral<Value>::value, Value> s2 = BernoulliDistributions<Distribution, Engine, Value>::bernoulli_distribution(expected_rand_values,
			d, e);
		Logger::WriteMessage("Executing Assert::AreNotEqual statement\n");
		for (auto i = 0; i != vec_length; ++i)
		{
			Assert::AreNotEqual(expected_rand_values[i], actual_rand_values[i], L"Bernoulli Distribution Inequality Test Failed", LINE_INFO());
		}
		Logger::WriteMessage("End of bernoulli_distribution() Test\n");
	}
//-------------------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------------------//
	template<typename Distribution, typename Engine, typename Value> void TestBernoulliDistributionInEquality(const int vec_length,
		Distribution d, Engine e)
	{
		Logger::WriteMessage("Begin Test of bernoulli_distribution() - Testing for InEquality\n");
		Logger::WriteMessage("Allocating data vectors\n");
		std::vector<Value> actual_rand_values(vec_length);
		std::vector<Value> expected_rand_values(vec_length);
		Logger::WriteMessage("Calling bernoulli distribution functions\n");
		std::enable_if<std::is_integral<Value>::value, Value> s1 = BernoulliDistributions<Distribution, Engine, Value>::bernoulli_distribution(actual_rand_values,
			d, e);
		std::enable_if<std::is_integral<Value>::value, Value> s2 = BernoulliDistributions<Distribution, Engine, Value>::bernoulli_distribution(expected_rand_values,
			d, e);
		Logger::WriteMessage("Executing Assert::AreEqual statement\n");
		for (auto i = 0; i != vec_length; ++i)
		{
			Assert::AreEqual(expected_rand_values[i], actual_rand_values[i], L"Bernoulli Distribution Equality Test Failed", LINE_INFO());
		}
		Logger::WriteMessage("End of bernoulli_distribution() Test\n");
	}
//------------------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------------------//
	TEST_CLASS(FDLIB_CPP_WRAPPER_TESTS)
		
	{
	public:
		
		TEST_METHOD(UniforDistributionInEqualityTest)
		{
			// TODO: Your test code here
			/*std::vector<unsigned int> rand_vec(100), rand_vec2(100); std::vector<double> rand_vec3(100);
			size_t seed = 0xFFFFFFFF;
			unsigned int status = _rdrand32_step(&seed);
			unsigned int minval = 10; unsigned int maxval = 65553;
			const int choose = 3;
			std::uniform_int_distribution<unsigned int> uid(10, 20);
			std::knuth_b ke;
		std::enable_if<std::is_integral<unsigned int>::value,unsigned int> s =	UniformDistributions<std::uniform_int_distribution<unsigned int>, std::knuth_b, unsigned int>::uni_int_distribution(rand_vec, uid, ke);
			UniDistributions<unsigned int>::uni_int_distribution(rand_vec2, minval, maxval, seed, 3);
			for (auto i = 0; i != rand_vec.size(); ++i)
				std::cout << rand_vec2[i] << std::endl;

			UniDistributions<double>::uni_real_distribution(rand_vec3, 1.0, 10.0, seed, 1);
			std::cout << "uniform_real_distribution\n" << std::endl;
			for (auto i = 0; i != rand_vec3.size(); ++i)
				std::cout << rand_vec3[i] << std::endl;*/

			unsigned int init_seed = 0xFFFFFFFF;
			unsigned int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{

			throw std::runtime_error (std::string("RDRAND Failed \n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::knuth_b kbe(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>,std::knuth_b, unsigned int>(vec_length, ud, kbe, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::knuth_b adapter engine.
		*/
		TEST_METHOD(UniforDistribution_knuth_b_Engine_Test)
		{
			Logger::WriteMessage("std::knuth_b Engine Test\n");
			unsigned int init_seed = 0xFFFFFFFF;
			unsigned int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{

				throw std::runtime_error(std::string("RDRAND Failed \n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::knuth_b kbe(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			Logger::WriteMessage("Calling Test function\n");
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::knuth_b, unsigned int>(vec_length, ud, kbe, min_val, max_val);
		}

		/*
		 *  TEST_METHOD for Uniform Distribution initialized by std::default_random_engine.
		 */
		TEST_METHOD(UniformDistribution_default_random_engine_Engine_Test)
		{
			Logger::WriteMessage("default_random_engine Engine Test\n");
			unsigned int init_seed = 0xFFFFFFFF;
			unsigned int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{

				throw std::runtime_error(std::string("Hardware RDRAND Failed \n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::default_random_engine dre(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			Logger::WriteMessage("Calling Test function\n");
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::default_random_engine, unsigned int>(vec_length, ud, dre, min_val, max_val);
		}

		/*
		 *  TEST_METHOD for Uniform Distribution initialized by std::minstd_rand engine.
		 */
		TEST_METHOD(UniformDistribution_std_minstd_rand_Engine_Test)
		{
			Logger::WriteMessage("std::minstd_rand Engine Test\n");
			unsigned int init_seed = 0xFFFFFFFF;
			unsigned int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string(" Hardware RDRAND Failed \n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::minstd_rand msr(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			Logger::WriteMessage("Calling Test function\n");
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::minstd_rand, unsigned int>(vec_length, ud, msr, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::minstd_rand0 engine.
		*/
		TEST_METHOD(UniformDistribution_std_minstd_rand0_Engine_Test)
		{
			Logger::WriteMessage("std::minstd_rand0 Engine Test\n");
			unsigned int init_seed = 0xFFFFFFFF;
			unsigned int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed \n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::minstd_rand0 msr0(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			Logger::WriteMessage("Calling Test function\n");
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::minstd_rand0, unsigned int>(vec_length, ud, msr0, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::mt19937 engine.
		*/
		TEST_METHOD(UniformDistribution_std_mt19937_Engine_Test)
		{
			Logger::WriteMessage("std::mt19937 Engine Test\n");
			unsigned int init_seed = 0xFFFFFFFF;
			unsigned int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::mt19937 mt(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			Logger::WriteMessage("Calling Test Function\n");
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::mt19937, unsigned int>(vec_length, ud, mt, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::mt19937_64 engine.
		*/

		TEST_METHOD(UniformDistribution_std_mt19937_64_Engine_Test)
		{
			Logger::WriteMessage("std::mt19937_64 Engine Test\n");
			unsigned int init_seed = 0xFFFFFFFF;
			unsigned int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::mt19937_64 mt64(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			Logger::WriteMessage("Calling Test Function\n");
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::mt19937_64, unsigned int>(vec_length, ud, mt64, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::ranlux24_base engine.
		*/

		TEST_METHOD(UniformDistribution_std_ranlux24_base_Engine_Test)
		{
			Logger::WriteMessage("std::ranlux24_base Engine Test\n");
			unsigned int init_seed = 0xFFFFFFFF;
			unsigned int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::ranlux24_base rl24_base(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			Logger::WriteMessage("Calling Test Function\n");
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::ranlux24_base, unsigned int>(vec_length, ud, rl24_base, min_val, max_val);
		}


		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::ranlux48_base engine.
		*/

		TEST_METHOD(UniformDistribution_std_ranlux48_base_Engine_Test)
		{
			Logger::WriteMessage("std::ranlux48_base Engine Test\n");
			unsigned int init_seed = 0x0;
			 int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::ranlux48_base rl48_base(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::ranlux48_base, unsigned int>(vec_length, ud, rl48_base, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::ranlux24 engine.
		*/

		TEST_METHOD(UniformDistribution_std_ranlux24_Engine_Test)
		{
			Logger::WriteMessage("std::ranlux24 Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::ranlux24 rl24(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::ranlux24, unsigned int>(vec_length, ud, rl24, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::linear_congruential_engine.
		*/

		TEST_METHOD(UniformDistribution_std_linear_congruential_engine_Engine_Test)
		{
			Logger::WriteMessage("std::linear_congruential_engine Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			const unsigned int val1 = 64U;
		    const	unsigned int val2 = 65536U;
			const unsigned int val3 = 2147483648U;
			std::linear_congruential_engine<unsigned int, val1, val2, val3> lce(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::linear_congruential_engine<unsigned int,val1,val2,val3>, unsigned int>(vec_length, ud, lce, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::marsenne_twister_engine.
		*/

		TEST_METHOD(UniformDistribution_std_marsenne_twister_engine_Engine_Test)
		{
			Logger::WriteMessage("std::marsenne_twister_engine Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100; //TO DO , Correct std::marsenne_twister_engine arguments!!
			const unsigned int val1 = 32U;
			const unsigned int val2 = 624U;
			const unsigned int val3 = 397U;
			const unsigned int val4 = 31U;
			const unsigned int val5 = 2567483615U;
			const unsigned int val6 = 11U;
			const unsigned int val7 = 4294967295U;
			const unsigned int val8 = 7U;
			const unsigned int val9 = 2636928640U;
			const unsigned int val10 = 15U;
			const unsigned int val11 = 4022730752U;
			const unsigned int val12 = 18U;
			const unsigned int val13 = 1812433253U;
			//std::mt19937 m;
			std::mersenne_twister_engine<unsigned int, val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12,val13> mte(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::mersenne_twister_engine<unsigned int, val1, val2, val3, val4,
				val5, val6, val7, val8, val9, val10, val11, val12, val13>, unsigned int>(vec_length, ud, mte, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::subtract_with_carry.
		*/

		/*TEST_METHOD(UniformDistribution_std_subtract_with_carry_Engine_Test)
		{
			Logger::WriteMessage("std::subtract_with_carry Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			const unsigned int val1 = 24U;
			const unsigned int val2 = 12U;
			const unsigned int val3 = 65U;
			
			std::subtract_with_carry<unsigned int, val1, val2, val3> swc(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::subtract_with_carry<unsigned int, val1, val2, val3>, unsigned int>(vec_length,
				ud, swc, min_val, max_val);
		}*/

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::subtract_with_carry_01.
		*/

		/*TEST_METHOD(UniformDistribution_std_subtract_with_carry_01_Engine_Test)
		{
			Logger::WriteMessage("std::subtract_with_carry_01 Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			const unsigned int val1 = 24U;
			const unsigned int val2 = 10U;
			const unsigned int val3 = 24U;
			
			std::subtract_with_carry_01<unsigned int, val1, val2, val3> swc(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::subtract_with_carry_01<unsigned int, val1, val2, val3>, unsigned int>(vec_length,
				ud, swc, min_val, max_val);
		}*/

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::ranlux3.
		*/

		TEST_METHOD(UniformDistribution_std_ranlux3_Engine_Test)
		{
			Logger::WriteMessage("std::ranlux3 Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			//std::ranlux3 rl3(init_seed);
			//std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			//TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::ranlux3, unsigned int>(vec_length, ud, rl3, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::ranlux3_01.
		*/
		TEST_METHOD(UniformDistribution_std_ranlux3_01_Engine_Test)
		{
			Logger::WriteMessage("std::ranlux3_01 Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			//std::ranlux3_01 rl3(vec_length);
			//std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			//TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::ranlux3_01, unsigned int>(vec_length, ud, rl3, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::ranlux4.
		*/

		TEST_METHOD(UniformDistribution_std_ranlux4_Engine_Test)
		{
			Logger::WriteMessage("std::ranlux4 Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			//std::ranlux4 rl4(init_seed);
			//std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			//TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::ranlux4, unsigned int>(vec_length, ud, rl4, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::ranlux4_01.
		*/

		TEST_METHOD(UniformDistribution_std_ranlux4_01)
		{
			Logger::WriteMessage("std::ranlux4_01 Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			//std::ranlux4_01 rl4(init_seed);
			//std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			//TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::ranlux4_01, unsigned int>(vec_length, ud, rl4, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::ranlux64_base_01.
		*/

		TEST_METHOD(UniformDistribution_std_ranlux64_base_01)
		{
			Logger::WriteMessage("std::ranlux64_base_01 Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			std::ranlux64_base_01 rl64b01(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::ranlux64_base_01, unsigned int>(vec_length, ud, rl64b01, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::discard_block_engine<knuth_b>.
		*/

		TEST_METHOD(UniformDistribution_std_discard_block_engine_Engine_Test)
		{
			Logger::WriteMessage("std::discard_block_engine<std::knuth_b> Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			const unsigned int val1 = 256U;
			const unsigned int val2 = 14U;
			std::discard_block_engine<std::knuth_b, val1, val2> dbe(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::discard_block_engine<std::knuth_b, val1, val2>, unsigned int>(vec_length, ud, dbe, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::discard_block_engine<minstd_rand>.
		*/

		TEST_METHOD(UniformDistribution_std_discard_block_engine_minstd_rand_Engine_Test)
		{
			Logger::WriteMessage("std::discard_block_engine<minstd_rand> Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			const unsigned int val1 = 65U;
			const unsigned int val2 = 14U;
			std::discard_block_engine<std::minstd_rand, val1, val2> dbe(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality < std::uniform_int_distribution<unsigned int>, std::discard_block_engine<std::minstd_rand, val1, val2>, unsigned int>(vec_length,
				ud, dbe, min_val, max_val);
		}
		
		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::discard_block_engine<minstd_rand>.
		*/

		TEST_METHOD(UniformDistribution_std_discarcd_block_engine_minstd_rand0_Engine_Test)
		{
			Logger::WriteMessage("std::discard_block_engine<std::minstd_rand0> Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			const unsigned int val1 = 64U;
			const unsigned int val2 = 12U;
			std::discard_block_engine<std::minstd_rand0, val1, val2> dbe(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::discard_block_engine<std::minstd_rand0, val1, val2>, unsigned int>(vec_length,
				ud, dbe, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::independent_bits_engine<std::linear_congruential_engine>.
		*/

		TEST_METHOD(UniformDistribution_std_independent_bits_engine_Engine_Test)
		{
			Logger::WriteMessage("std::independent_bits_engine Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			const unsigned int val1 = 14U;
			const unsigned int val2 = 25U;
			const unsigned int val3 = 128U;
			const unsigned int num_bits = 30U;
			std::independent_bits_engine < std::linear_congruential_engine < unsigned int, val1, val2, val3>, num_bits,unsigned int> ibe(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality < std::uniform_int_distribution<unsigned int>, std::independent_bits_engine<std::linear_congruential_engine<unsigned int, val1, val2, val3>, num_bits,unsigned int>,unsigned int>(vec_length,
				ud, ibe, min_val, max_val);
		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::independent_bits_engine<std::minstd_rand>.
		*/

		TEST_METHOD(UniformDistribution_std_independent_bits_engine_minstd_rand_Engine_Test)
		{
			Logger::WriteMessage("std::independent_bits_engine<std::minstd_rand> Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			const unsigned int num_bits = 30U;
			std::independent_bits_engine<std::minstd_rand, num_bits, unsigned int> ibe(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::independent_bits_engine<std::minstd_rand, num_bits, unsigned int>, unsigned int>(vec_length,
				ud, ibe, min_val, max_val);

		}

		/*
		*  TEST_METHOD for Uniform Distribution initialized by std::independent_bits_engine<std::minstd_rand0>.
		*/

		TEST_METHOD(UniformDistribution_std_independent_bits_engine_minstd_rand0_Engine_Test)
		{
			Logger::WriteMessage("std::independent_bits_engine<std::minstd_rand0> Engine Test\n");
			unsigned int init_seed = 0x0;
			int status = 0xFFFFFFFF;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{
				throw std::runtime_error(std::string("Hardware RDRAND Failed\n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators and distribution functors started\n");
			const unsigned int min_val = 10;
			const unsigned int max_val = 100;
			const unsigned int vec_length = 100;
			const unsigned int num_bits = 25U;
			std::independent_bits_engine<std::minstd_rand0, num_bits, unsigned int> ibe(init_seed);
			std::uniform_int_distribution<unsigned int> ud(min_val, max_val);
			TestUniformDistributionEquality<std::uniform_int_distribution<unsigned int>, std::independent_bits_engine<std::minstd_rand0, num_bits, unsigned int>, unsigned int>(vec_length,
				ud, ibe, min_val, max_val);
		}

		TEST_METHOD(BernoulliDistributionInEqualityTest)
		{
			unsigned int init_seed = 0xFFFFFFFF;
			unsigned int status = 0xFFFFFFFF;
			const double probability = 0.5;
			const unsigned int vec_length = 100;
			if ((status = _rdrand32_step(&init_seed)) == 0)
			{

				throw std::runtime_error(std::string("RDRAND Failed \n"));
			}
			Logger::WriteMessage("Initialization of seed completed successfully\n");
			Logger::WriteMessage("Initialization of random number generators started\n");
			std::knuth_b kbe(init_seed);
			std::bernoulli_distribution bd{ probability };
			TestBernoulliDistributionEquality<std::bernoulli_distribution, std::knuth_b, unsigned int>(vec_length, bd, kbe);
		}

	};
}
