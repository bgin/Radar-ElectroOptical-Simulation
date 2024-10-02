#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <immintrin.h>
#include "GMS_stats_functions.h"

bool    arith_mean(double ** population, const int len_i, const int len_j, double * mean)
{
	bool b_result;
	if ((population == NULL) && (len_i <= 0) && (len_j <= 0) )
	{
		return b_result = false;
	}

	for (int i = 0; i < len_i; ++i)
	{
		for (int j = 0; j < len_j; ++j)
		{
			*mean += population[i][j];
		}
	}
	*mean /= (len_i * len_j);
	return b_result = true;
}

bool    geo_mean(double ** population, const int len_i, const int len_j, double * mean)
{
	bool b_result;
	if ((population == NULL) && (len_i <= 0) && (len_j <= 0) )
	{
		return b_result = false;
	}

	for (int i = 0; i < len_i; ++i)
	{
		for (int j = 0; j < len_j; ++j)
		{
			*mean *= population[i][j];
		}
	}
	*mean = pow(*mean, 1.0 / (len_i * len_j));
	return b_result = true;
}

bool    std_deviation(double ** population, const int len_i, const int len_j, double * mean)
{
	bool b_result;
	if ((population == NULL) && (len_i <= 0) && (len_j <= 0) )
	{
		return b_result = false;
	}

	double average = 0.0L;
	for (int i = 0; i < len_i; ++i)
	{
		for (int j = 0; j < len_j; ++j)
		{
			average += population[i][j];
		}
	}

	average /= (len_i * len_j);

	for (int i = 0; i < len_i; ++i)
	{
		for (int j = 0; j < len_j; ++j)
		{
			population[i][j] = (population[i][j] - average) * (population[i][j] - average);
			*mean += population[i][j];
		}
	}
	is_population_modified = true;
	*mean = sqrt(*mean / (len_i * len_j));

	return b_result = true;
}

bool   geomet_median(double **population, const int len_i, const int len_j, double * mean)
{
	bool b_result;
	if ((population == NULL) && (len_i <= 0) && (len_j <= 0) )
	{
		return b_result = false;
	}

	if (!is_population_modified)
	{
		double oldmean = 0.0L;
		double nom = 0.0L;
		double denom = 0.0L;
		for (int i = 0; i < len_i; ++i)
		{
			for (int j = 0; j < len_j; ++j)
			{
				nom += population[i][j] / (sqrt((population[i][j] - oldmean) * (population[i][j] - oldmean)));
				denom += 1.0 / (sqrt((population[i][j] - oldmean) * (population[i][j] - oldmean)));
				if ((*mean = fabs(nom / denom)) > eps)
				{
					oldmean = *mean;
				}
				else
				{
					goto less_than_eps;
				}
			}
		}
		return b_result = true;
	}
	else
	{
		return b_result = false;
	}

less_than_eps:
	return b_result = false;
}


bool  skew(double ** population, const int len_i, const int len_j, double * mean)
{
	bool b_result;
	if ((population == NULL) && (len_i <= 0) && (len_j <= 0) )
	{
		return b_result = false;
	}

	double average = 0.0L;
	double nom = 0.0L;
	double denom = 0.0L;
	
	for (int i = 0; i < len_i; ++i)
	{
		for (int j = 0; j < len_j; ++j)
		{
			average += population[i][j];
		}
	}

	average /= (len_i * len_j);

	for (int i = 0; i < len_i; ++i)
	{
		for (int j = 0; j < len_j; ++j)
		{
			nom += (population[i][j] - average) * (population[i][j] - average) * (population[i][j] - average) / (len_i * len_j);
			denom += (population[i][j] - average) * (population[i][j] - average) / (len_i * (len_j - 1));
			if (fabs(denom) > eps)
			{
				*mean = nom / pow(denom, 1.5);
			}
			else
			{
				goto less_than_eps;
			}
			
		}
	}
	return b_result = true;

less_than_eps:
	return b_result = false;

}

bool  arith_mean_wf(double ** population, const int len_i, const int len_j, double * mean)
{
	bool b_result = false;
	if ((population == NULL) && (len_i <= 0) && (len_j <= 0) && (fname == NULL))
	{
		return b_result;
	}

	else
	{
		
		for (int i = 0; i < len_i; ++i)
		{
			for (int j = 0; j < len_j; ++j)
			{
				*mean += population[i][j];
				
			}
		}
		*mean /= (len_i * len_j);
		return b_result = true;
	}
}

bool  geo_mean_wf(double ** population, const int len_i, const int len_j, double * mean)
{
	bool b_result = false;
	if ((population == NULL) && (len_i <= 0) && (len_j <= 0) && (fname == NULL))
	{
		return b_result;
	}
	else
	{
		
		for (int i = 0; i < len_i; ++i)
		{
			for (int j = 0; j < len_j; ++j)
			{
				*mean *= population[i][j];
				
			}
		}
		*mean = pow(*mean, 1.0 / (len_i * len_j));
		return b_result = true;
	}
}

bool  std_deviation_wf(double ** population, const int len_i, const int len_j, double * mean)
{
	bool b_result = false;
	if ((population == NULL) && (len_i <= 0) && (len_j <= 0) && (fname == NULL))
	{
		return b_result;
	}
	else
	{
		
		double average = 0.0L;
		for (int i = 0; i < len_i; ++i)
		{
			for (int j = 0; j < len_j; ++j)
			{
				average += population[i][j];
			}
		}

		average /= (len_i * len_j);

		for (int i = 0; i < len_i; ++i)
		{
			for (int j = 0; j < len_j; ++j)
			{
				population[i][j] = (population[i][j] - average) * (population[i][j] - average);
				*mean += population[i][j];
				

			}
		}
		is_population_modified = true;
		*mean = sqrt(*mean / (len_i * len_j));
		return b_result = true;
		
	}
}

bool  geo_median_wf(double ** population, const int len_i, const int len_j, double * mean)
{
	bool b_result = false;
	if ((population == NULL) && (len_i <= 0) && (len_j <= 0) && (fname == NULL))
	{
		return b_result;
	}
	FILE * fp;
	else
	{
		
		if (!is_population_modified)
		{
			double oldmean = 0.0L;
			double nom = 0.0L;
			double denom = 0.0L;
			for (int i = 0; i < len_i; ++i)
			{
				for (int j = 0; j < len_j; ++j)
				{
					nom += population[i][j] / (sqrt((population[i][j] - oldmean) * (population[i][j] - oldmean)));
					denom += 1.0 / (sqrt((population[i][j] - oldmean) * (population[i][j] - oldmean)));
					if ((*mean = fabs(nom / denom)) > eps)
					{
						oldmean = *mean;
						
					}
					else
					{
						goto less_than_eps;
					}
				}
			}
			return b_result = true;
		}
		else
		{
			
			return b_result = false;
		}
	less_than_eps:
		
		return b_result = false;
	}
}

bool  skew_wf(double population[][100] , const int len_i, const int len_j, double * mean,)
{
	bool b_result = false;
	if ((population == NULL) && (len_i <= 0) && (len_j <= 0) && (fname == NULL))
	{
		return b_result;
	}
	else
	{
		
		double average = 0.0L;
		double nom = 0.0L;
		double denom = 0.0L;

		for (int i = 0; i < len_i; ++i)
		{
			for (int j = 0; j < len_j; ++j)
			{
				average += population[i][j];
			}
		}

		average /= (len_i * len_j);
		for (int i = 0; i < len_i; ++i)
		{
			for (int j = 0; j < len_j; ++j)
			{
				    nom += (population[i][j] - average) * (population[i][j] - average) * (population[i][j] - average);
					denom += (population[i][j] - average) * (population[i][j] - average);
				
			}
		}
		nom /= (len_i * len_j);
		denom /= (len_i * (len_j - 1));
		if (fabs(denom) > eps)
		{
			*mean = nom / pow(denom, 1.5);
		}
		else
		{
			goto less_than_eps;
		}
		
		return b_result = true;

	less_than_eps:
		
		return b_result = false;
	}
}

void Test_skew_wf()
{
	unsigned long long ts0 = 0ULL; unsigned long long  ts1 = 0ULL;
	ts0 = _rdtsc();
	//Sleep(1);
	ts1 = _rdtsc();
	unsigned long long elapsed = ts1 - ts0;
	printf("TSC cycles = %llu\n", elapsed);
	const char *name = "sin_function_timing3.csv";
	const int max_outer_iters = 10;
	const int max_inner_iters = 100;
	double res = 0.0L;
	double mean_value = 0.0L;
	long long clock_ticks_start[max_outer_iters][max_inner_iters] = { 0LL };
	long long clock_ticks_stop[max_outer_iters][max_inner_iters] = { 0LL };
	double ns[max_outer_iters][max_inner_iters] = { 0.0L };
	unsigned int rand = 0U;
	unsigned int status = 0U;
#pragma noinline
	for (int i = 0; i < max_outer_iters; ++i)
	{
		for (int j = 0; j < max_inner_iters; ++j)
		{
			if ((status = _rdrand32_step(&rand)) == 0)
			{
				printf("RDRAND Failed\n");
				goto rdrand_error;
			}
			else
			{
				clock_ticks_start[i][j] = _rdtsc();
				res = sin((double)rand);
				clock_ticks_stop[i][j] = _rdtsc();
				ns[i][j] = (double)(clock_ticks_stop[i][j] - clock_ticks_start[i][j]) / (elapsed * 1.0e-6);
			}
		}
	}
	bool b_res = false;
	b_res = skew_wf(ns, max_outer_iters, max_inner_iters, &mean_value, name);
	printf("Skewennes = %.15f\n", mean_value);
	if (!b_res)
		printf("skew_wf() failed: %d\n", b_res);
	return;

rdrand_error:
	printf("RDRAND HW error \n");
	return;
}
