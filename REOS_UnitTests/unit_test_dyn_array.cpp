#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"

/*
   icpc -o unit_test_dyn_array -fp-model fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_dyn_array.h unit_test_dyn_array.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -falign-functions=32   GMS_config.h GMS_malloc.h GMS_dyn_array.h unit_test_dyn_array.cpp

*/


void unit_test_darray_c4_t_2nd_Ctor();

void unit_test_darray_c4_t_2nd_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{18ull*4096ull};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_darray_c4_t_2nd_Ctor.csv"};
       std::clock_t seeds[2] = {0ULL};
       std::uniform_real_distribution<float> distros[2]{};
       FILE * fp{NULL};
       darray_c4_t test_2nd_Ctor = darray_c4_t(nx);
       char * ctor_name{gms::common::demangle(typeid(test_2nd_Ctor).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(test_2nd_Ctor).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       test_2nd_Ctor.info_size_alignment();
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");
       for(int __i{0}; __i != 2; ++__i)
       {
           seeds[__i] = std::clock();
           distros[__i]  = std::uniform_real_distribution<float>(0.0f,1.0f);
       }
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       for(std::size_t __i{0ULL}; __i != test_2nd_Ctor.mnx; ++__i)
       {
           const std::complex<float> rc{distros[0].operator()(rdev1),distros[1].operator()(rdev2)};
           test_2nd_Ctor.m_data[__i] = rc;
       }
       
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_data to file=%s\n",fname);
       for(std::size_t __i{0}; __i != test_2nd_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_2nd_Ctor.m_data[__i].real(),
                                        test_2nd_Ctor.m_data[__i].imag());
       }
     
      fclose(fp);
      printf("[UNIT-TEST #1]: darray_c4_t_2nd_Ctor() -- Ended correctly!!\n");
}

void unit_test_darray_c4_t_3rd_Ctor();

void unit_test_darray_c4_t_3rd_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{14ull*4096ull};
      
       constexpr int32_t prot{PROT_READ | PROT_WRITE};
       constexpr int32_t flags{MAP_ANONYMOUS | MAP_PRIVATE};
       constexpr int32_t fd{-1};
       constexpr off_t offset{0};
       constexpr int32_t fsize{0};
       const char * fname{"UNIT_TEST_Output_Random_Data_darray_c4_t_3rd_Ctor.csv"};
       std::clock_t seeds[2] = {0ULL};
       std::uniform_real_distribution<float> distros[2]{};
       FILE * fp{NULL};
       int32_t status;
       darray_c4_t test_3rd_Ctor = darray_c4_t(nx,prot,flags,fd,offset,fsize);
       char * ctor_name{gms::common::demangle(typeid(test_3rd_Ctor).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(test_3rd_Ctor).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       test_3rd_Ctor.info_size_alignment();
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");                                
       for(int32_t __i{0}; __i != 2; ++__i)
       {
           seeds[__i]   = std::clock();
           distros[__i] = std::uniform_real_distribution<float>(0.0f,1.0f); 
       }  
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       for(std::size_t __i{0ULL}; __i != test_3rd_Ctor.mnx; ++__i)
       {
           test_3rd_Ctor.m_data[__i] = std::complex<float>{distros[0].operator()(rdev1),
                                                           distros[1].operator()(rdev2)};
       }
       
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_data to file=%s \n",fname);
       for(std::size_t __i{0}; __i != test_3rd_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_3rd_Ctor.m_data[__i].real(),
                                        test_3rd_Ctor.m_data[__i].imag());
       }
       
       fclose(fp);
       printf("[UNIT-TEST #2]: darray_c4_t_3rd_Ctor() -- Ended correctly!!\n");
}


void unit_test_darray_c4_t_4th_Ctor();

void unit_test_darray_c4_t_4th_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{11ull*4096ull};
      
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_darray_c4_t_4th_Ctor.csv"};
       std::clock_t seeds[2] = {0ULL};
       std::uniform_real_distribution<float> distros[2]{};
       FILE * fp{NULL};
       std::vector<std::complex<float>> data = std::vector<std::complex<float>>(nx);
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");                                
       for(int32_t __i{0}; __i != 2; ++__i)
       {
           seeds[__i]   = std::clock();
           distros[__i] = std::uniform_real_distribution<float>(0.0f,1.0f); 
       }  
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       for(std::size_t __i{0ULL}; __i != nx; ++__i)
       {
           data[__i] = std::complex<float>{distros[0].operator()(rdev1),
                                          distros[1].operator()(rdev2)};
       }
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       darray_c4_t test_4th_Ctor = darray_c4_t(data);
       char * ctor_name{gms::common::demangle(typeid(test_4th_Ctor).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(test_4th_Ctor).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       test_4th_Ctor.info_size_alignment();
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST #3]: Dumping: m_data to file=%s \n",fname);
       for(std::size_t __i{0}; __i != test_4th_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_4th_Ctor.m_data[__i].real(),
                                        test_4th_Ctor.m_data[__i].imag());
       }
       
       fclose(fp);
       printf("[UNIT-TEST #3]: darray_c4_t_4th_Ctor() -- Ended correctly!!\n");
}


void unit_test_darray_c4_t_5th_Ctor();

void unit_test_darray_c4_t_5th_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{11ull*4096ull};
      
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_darray_c4_t_5th_Ctor.csv"};
       std::clock_t seeds[2] = {0ULL};
       std::uniform_real_distribution<float> distros[2]{};
       FILE * fp{NULL};
       std::valarray<std::complex<float>> data = std::valarray<std::complex<float>>(nx);
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");                                
       for(int32_t __i{0}; __i != 2; ++__i)
       {
           seeds[__i]   = std::clock();
           distros[__i] = std::uniform_real_distribution<float>(0.0f,1.0f); 
       }  
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       for(std::size_t __i{0ULL}; __i != nx; ++__i)
       {
           data[__i] = std::complex<float>{distros[0].operator()(rdev1),
                                          distros[1].operator()(rdev2)};
       }
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       darray_c4_t test_5th_Ctor = darray_c4_t(data);
       char * ctor_name{gms::common::demangle(typeid(test_5th_Ctor).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(test_5th_Ctor).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       test_5th_Ctor.info_size_alignment();
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST #4]: Dumping: m_data to file=%s \n",fname);
       for(std::size_t __i{0}; __i != test_5th_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_5th_Ctor.m_data[__i].real(),
                                        test_5th_Ctor.m_data[__i].imag());
       }
       
       fclose(fp);
       printf("[UNIT-TEST #4]: darray_c4_t_5th_Ctor() -- Ended correctly!!\n");
}



int main()
{
   unit_test_darray_c4_t_2nd_Ctor();
   unit_test_darray_c4_t_3rd_Ctor();
   unit_test_darray_c4_t_4th_Ctor();
   unit_test_darray_c4_t_5th_Ctor();
   return 0;
}
