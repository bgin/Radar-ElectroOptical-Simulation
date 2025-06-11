

#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_malloc.h" // for demangle
#include "GMS_swap_avx_unroll8x.h"

/*
   icpc -o unit_test_dswap_avx_unroll8x -fp-model fast=2 -ftz -ggdb -ipo  -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_swap_avx_unroll8x.h GMS_swap_avx_unroll8x.cpp unit_test_dswap_avx_unroll8x.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_swap_avx_unroll8x.h GMS_swap_avx_unroll8x.cpp unit_test_dswap_avx_unroll8x.cpp
                        
*/

void unit_test_dswap_u_ymm4r8_unroll8x();

void unit_test_dswap_u_ymm4r8_unroll8x()
{
    constexpr std::size_t buf_len{4096ull};
    double buf_x[buf_len];
    double buf_y[buf_len];
    std::clock_t seed_x{0ull};
    std::clock_t seed_y{0ull};
    std::clock_t seed{0ull};
    FILE * fp{NULL};
    const char * fname{"OUTPUT-DATA_dswap_u_ymm4r8_unroll8x.csv"};
    int32_t which{-1};
    seed = std::clock();
    auto rand{std::bind(std::uniform_int_distribution<int32_t>(0,4),
                         std::mt19937(seed))};
    which = rand();
    switch (which)
    {
        case 0 : 
        {
             std::uniform_real_distribution<double> uniform_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(uniform_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             fp = fopen(fname,"w+");
             if(!fp)
             {
                    std::perror("fopen failed to open a file -- TERMINATING!!");
                    std::terminate();
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::uniform_real_distribution<double>(0.0,1.0),
                                   std::mt19937(seed_x))};
             std::generate((double*)&buf_x[0],(double*)&buf_x[0]+buf_len,rand_x);
             seed_y = std::clock();
             auto rand_y{std::bind(std::uniform_real_distribution<double>(0.0,1.0),
                                    std::mt19937(seed_y))};
             std::generate((double*)&buf_y[0],(double*)&buf_y[0]+buf_len,rand_y);
             printf("[UNIT-TEST]: -- End of data initialization\n");
             printf("[UNIT-TEST]: -- START: dswap_u_ymm4r8_unroll8x\n");
             gms::math::dswap_u_ymm4r8_unroll8x(buf_len,&buf_x[0],1ull,&buf_y[0],1ull);
             for(std::size_t __i{0ull}; __i != buf_len; ++__i)
             {
                   fprintf(fp,"[buf_x]=%.16f, [buf_y]=%.16f\n",buf_x[__i],buf_y[__i]);
             }
             printf("[UNIT-TEST]: -- END: dswap_u_ymm4r8_unroll8x\n");   
             fclose(fp);  
        }                    
        break;
        case 1 : 
        {
             std::normal_distribution<double> normal_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(normal_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(normal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             fp = fopen(fname,"w+");
             if(!fp)
             {
                    std::perror("fopen failed to open a file -- TERMINATING!!");
                    std::terminate();
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::normal_distribution<double>(0.0,1.0),
                                   std::mt19937(seed_x))};
             std::generate((double*)&buf_x[0],(double*)&buf_x[0]+buf_len,rand_x);
             seed_y = std::clock();
             auto rand_y{std::bind(std::normal_distribution<double>(0.0,1.0),
                                    std::mt19937(seed_y))};
             std::generate((double*)&buf_y[0],(double*)&buf_y[0]+buf_len,rand_y);
             printf("[UNIT-TEST]: -- End of data initialization\n");
             printf("[UNIT-TEST]: -- START: dswap_u_ymm4r8_unroll8x\n");
             gms::math::dswap_u_ymm4r8_unroll8x(buf_len,&buf_x[0],1ull,&buf_y[0],1ull);
             for(std::size_t __i{0ull}; __i != buf_len; ++__i)
             {
                   fprintf(fp,"[buf_x]=%.16f, [buf_y]=%.16f\n",buf_x[__i],buf_y[__i]);
             }
             printf("[UNIT-TEST]: -- END: dswap_u_ymm4r8_unroll8x\n");   
             fclose(fp);  
        }
        break;
        case 2 : 
        {
             std::lognormal_distribution<double> lognormal_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(lognormal_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognormal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             fp = fopen(fname,"w+");
             if(!fp)
             {
                    std::perror("fopen failed to open a file -- TERMINATING!!");
                    std::terminate();
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::lognormal_distribution<double>(1.55,0.245),
                                   std::mt19937(seed_x))};
             std::generate((double*)&buf_x[0],(double*)&buf_x[0]+buf_len,rand_x);
             seed_y = std::clock();
             auto rand_y{std::bind(std::lognormal_distribution<double>(0.99f,1.110f),
                                    std::mt19937(seed_y))};
             std::generate((double*)&buf_y[0],(double*)&buf_y[0]+buf_len,rand_y);
             printf("[UNIT-TEST]: -- End of data initialization\n");
             printf("[UNIT-TEST]: -- START: dswap_u_ymm4r8_unroll8x\n");
             gms::math::dswap_u_ymm4r8_unroll8x(buf_len,&buf_x[0],1ull,&buf_y[0],1ull);
             for(std::size_t __i{0ull}; __i != buf_len; ++__i)
             {
                   fprintf(fp,"[buf_x]=%.16f, [buf_y]=%.16f\n",buf_x[__i],buf_y[__i]);
             }
             printf("[UNIT-TEST]: -- END: dswap_u_ymm4r8_unroll8x\n");   
             fclose(fp);  
        }
        break;
        case 3 : 
        {
             std::fisher_f_distribution<double> fisher_f_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(fisher_f_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(fisher_f_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             fp = fopen(fname,"w+");
             if(!fp)
             {
                    std::perror("fopen failed to open a file -- TERMINATING!!");
                    std::terminate();
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::fisher_f_distribution<double>(1.55,5.245),
                                   std::mt19937(seed_x))};
             std::generate((double*)&buf_x[0],(double*)&buf_x[0]+buf_len,rand_x);
             seed_y = std::clock();
             auto rand_y{std::bind(std::fisher_f_distribution<double>(5.99,1.110),
                                    std::mt19937(seed_y))};
             std::generate((double*)&buf_y[0],(double*)&buf_y[0]+buf_len,rand_y);
             printf("[UNIT-TEST]: -- End of data initialization\n");
             printf("[UNIT-TEST]: -- START: dswap_u_ymm4r8_unroll8x\n");
             gms::math::dswap_u_ymm4r8_unroll8x(buf_len,&buf_x[0],1ull,&buf_y[0],1ull);
             for(std::size_t __i{0ull}; __i != buf_len; ++__i)
             {
                   fprintf(fp,"[buf_x]=%.16f, [buf_y]=%.16f\n",buf_x[__i],buf_y[__i]);
             }
             printf("[UNIT-TEST]: -- END: dswap_u_ymm4r8_unroll8x\n");   
             fclose(fp);  
        }
        break;
        case 4 : 
        {
             std::cauchy_distribution<double> cauchy_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(cauchy_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(cauchy_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             fp = fopen(fname,"w+");
             if(!fp)
             {
                    std::perror("fopen failed to open a file -- TERMINATING!!");
                    std::terminate();
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::cauchy_distribution<double>(-0.55,3.245),
                                   std::mt19937(seed_x))};
             std::generate((double*)&buf_x[0],(double*)&buf_x[0]+buf_len,rand_x);
             seed_y = std::clock();
             auto rand_y{std::bind(std::cauchy_distribution<double>(-4.99,1.110),
                                    std::mt19937(seed_y))};
             std::generate((double*)&buf_y[0],(double*)&buf_y[0]+buf_len,rand_y);
             printf("[UNIT-TEST]: -- End of data initialization\n");
             printf("[UNIT-TEST]: -- START: dswap_u_ymm4r8_unroll8x\n");
             gms::math::dswap_u_ymm4r8_unroll8x(buf_len,&buf_x[0],1ull,&buf_y[0],1ull);
             for(std::size_t __i{0ull}; __i != buf_len; ++__i)
             {
                   fprintf(fp,"[buf_x]=%.16f, [buf_y]=%.16f\n",buf_x[__i],buf_y[__i]);
             }
             printf("[UNIT-TEST]: -- END: dswap_u_ymm4r8_unroll8x\n");   
             fclose(fp);  
        }
        break;
        default : 
             printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
             std::terminate();
    }
}

__attribute__((hot))
void unit_test_dswap_a_ymm4r8_unroll8x();

void unit_test_dswap_a_ymm4r8_unroll8x()
{
    constexpr std::size_t buf_len{4096ull};
    double __attribute__((aligned(32))) buf_x[buf_len];
    double __attribute__((aligned(32))) buf_y[buf_len];
    std::clock_t seed_x{0ull};
    std::clock_t seed_y{0ull};
    std::clock_t seed{0ull};
    FILE * fp{NULL};
    const char * fname{"OUTPUT-DATA_dswap_a_ymm4r8_unroll8x.csv"};
    int32_t which{-1};
    seed = std::clock();
    auto rand{std::bind(std::uniform_int_distribution<int32_t>(0,4),
                         std::mt19937(seed))};
    which = rand();
    switch (which)
    {
        case 0 : 
        {
             std::uniform_real_distribution<double> uniform_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(uniform_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             fp = fopen(fname,"w+");
             if(!fp)
             {
                    std::perror("fopen failed to open a file -- TERMINATING!!");
                    std::terminate();
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::uniform_real_distribution<double>(0.0,1.0),
                                   std::mt19937(seed_x))};
             std::generate((double*)&buf_x[0],(double*)&buf_x[0]+buf_len,rand_x);
             seed_y = std::clock();
             auto rand_y{std::bind(std::uniform_real_distribution<double>(0.0f,1.0f),
                                    std::mt19937(seed_y))};
             std::generate((double*)&buf_y[0],(double*)&buf_y[0]+buf_len,rand_y);
             printf("[UNIT-TEST]: -- End of data initialization\n");
             printf("[UNIT-TEST]: -- START: dswap_a_ymm4r8_unroll8x\n");
             gms::math::dswap_a_ymm4r8_unroll8x(buf_len,&buf_x[0],1ull,&buf_y[0],1ull);
             for(std::size_t __i{0ull}; __i != buf_len; ++__i)
             {
                   fprintf(fp,"[buf_x]=%.16f, [buf_y]=%.16f\n",buf_x[__i],buf_y[__i]);
             }
             printf("[UNIT-TEST]: -- END: dswap_a_ymm4r8_unroll8x\n");   
             fclose(fp);  
        }                    
        break;
        case 1 : 
        {
             std::normal_distribution<double> normal_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(normal_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(normal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             fp = fopen(fname,"w+");
             if(!fp)
             {
                    std::perror("fopen failed to open a file -- TERMINATING!!");
                    std::terminate();
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::normal_distribution<double>(0.0,1.0),
                                   std::mt19937(seed_x))};
             std::generate((double*)&buf_x[0],(double*)&buf_x[0]+buf_len,rand_x);
             seed_y = std::clock();
             auto rand_y{std::bind(std::normal_distribution<double>(0.0f,1.0f),
                                    std::mt19937(seed_y))};
             std::generate((double*)&buf_y[0],(double*)&buf_y[0]+buf_len,rand_y);
             printf("[UNIT-TEST]: -- End of data initialization\n");
             printf("[UNIT-TEST]: -- START: sswap_a_ymm4r8_unroll8x\n");
             gms::math::dswap_a_ymm4r8_unroll8x(buf_len,&buf_x[0],1ull,&buf_y[0],1ull);
             for(std::size_t __i{0ull}; __i != buf_len; ++__i)
             {
                   fprintf(fp,"[buf_x]=%.16f, [buf_y]=%.16f\n",buf_x[__i],buf_y[__i]);
             }
             printf("[UNIT-TEST]: -- END: dswap_a_ymm4r8_unroll8x\n");   
             fclose(fp);  
        }
        break;
        case 2 : 
        {
             std::lognormal_distribution<double> lognormal_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(lognormal_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognormal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             fp = fopen(fname,"w+");
             if(!fp)
             {
                    std::perror("fopen failed to open a file -- TERMINATING!!");
                    std::terminate();
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::lognormal_distribution<double>(1.55,0.245),
                                   std::mt19937(seed_x))};
             std::generate((double*)&buf_x[0],(double*)&buf_x[0]+buf_len,rand_x);
             seed_y = std::clock();
             auto rand_y{std::bind(std::lognormal_distribution<double>(0.99,1.110),
                                    std::mt19937(seed_y))};
             std::generate((double*)&buf_y[0],(double*)&buf_y[0]+buf_len,rand_y);
             printf("[UNIT-TEST]: -- End of data initialization\n");
             printf("[UNIT-TEST]: -- START: dswap_a_ymm4r8_unroll8x\n");
             gms::math::dswap_a_ymm4r8_unroll8x(buf_len,&buf_x[0],1ull,&buf_y[0],1ull);
             for(std::size_t __i{0ull}; __i != buf_len; ++__i)
             {
                   fprintf(fp,"[buf_x]=%.16f, [buf_y]=%.16f\n",buf_x[__i],buf_y[__i]);
             }
             printf("[UNIT-TEST]: -- END: dswap_a_ymm4r8_unroll8x\n");   
             fclose(fp);  
        }
        break;
        case 3 : 
        {
             std::fisher_f_distribution<double> fisher_f_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(fisher_f_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(fisher_f_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             fp = fopen(fname,"w+");
             if(!fp)
             {
                    std::perror("fopen failed to open a file -- TERMINATING!!");
                    std::terminate();
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::fisher_f_distribution<double>(1.55,5.245),
                                   std::mt19937(seed_x))};
             std::generate((double*)&buf_x[0],(double*)&buf_x[0]+buf_len,rand_x);
             seed_y = std::clock();
             auto rand_y{std::bind(std::fisher_f_distribution<double>(5.99,1.110),
                                    std::mt19937(seed_y))};
             std::generate((double*)&buf_y[0],(double*)&buf_y[0]+buf_len,rand_y);
             printf("[UNIT-TEST]: -- End of data initialization\n");
             printf("[UNIT-TEST]: -- START: dswap_a_ymm4r8_unroll8x\n");
             gms::math::dswap_a_ymm4r8_unroll8x(buf_len,&buf_x[0],1ull,&buf_y[0],1ull);
             for(std::size_t __i{0ull}; __i != buf_len; ++__i)
             {
                   fprintf(fp,"[buf_x]=%.16f, [buf_y]=%.16f\n",buf_x[__i],buf_y[__i]);
             }
             printf("[UNIT-TEST]: -- END: dswap_a_ymm4r8_unroll8x\n");   
             fclose(fp);  
        }
        break;
        case 4 : 
        {
             std::cauchy_distribution<double> cauchy_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(cauchy_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(cauchy_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             fp = fopen(fname,"w+");
             if(!fp)
             {
                    std::perror("fopen failed to open a file -- TERMINATING!!");
                    std::terminate();
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::cauchy_distribution<double>(-0.55,3.245),
                                   std::mt19937(seed_x))};
             std::generate((double*)&buf_x[0],(double*)&buf_x[0]+buf_len,rand_x);
             seed_y = std::clock();
             auto rand_y{std::bind(std::cauchy_distribution<double>(-4.99,1.110),
                                    std::mt19937(seed_y))};
             std::generate((double*)&buf_y[0],(double*)&buf_y[0]+buf_len,rand_y);
             printf("[UNIT-TEST]: -- End of data initialization\n");
             printf("[UNIT-TEST]: -- START: dswap_a_ymm4r8_unroll8x\n");
             gms::math::dswap_a_ymm4r8_unroll8x(buf_len,&buf_x[0],1ull,&buf_y[0],1ull);
             for(std::size_t __i{0ull}; __i != buf_len; ++__i)
             {
                   fprintf(fp,"[buf_x]=%.16f, [buf_y]=%.16f\n",buf_x[__i],buf_y[__i]);
             }
             printf("[UNIT-TEST]: -- END: sswap_a_ymm4r8_unroll8x\n");   
             fclose(fp);  
        }
        break;
        default : 
             printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
             std::terminate();
    }
}

int main()
{   
    unit_test_dswap_u_ymm4r8_unroll8x();
    unit_test_dswap_a_ymm4r8_unroll8x();
    return 0;
}