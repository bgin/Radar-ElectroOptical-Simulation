#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"

/*
   icpc -o unit_test_dyn_array_r8 -fp-model fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_dyn_array.h unit_test_dyn_array_r8.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -falign-functions=32   GMS_config.h GMS_malloc.h GMS_dyn_array.h unit_test_dyn_array_r8.cpp

*/


void unit_test_darray_r8_t_2nd_Ctor();

void unit_test_darray_r8_t_2nd_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{18ull*4096ull};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_darray_r8_t_2nd_Ctor.csv"};
       std::clock_t seed{0ULL};
       std::uniform_real_distribution<double> distro{};
       FILE * fp{NULL};
       darray_r8_t test_2nd_Ctor = darray_r8_t(nx);
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
       seed = std::clock();
       distro  = std::uniform_real_distribution<double>(0.0,1.0);
       auto rdev1{std::mt19937(seed)};
       for(std::size_t __i{0ULL}; __i != test_2nd_Ctor.mnx; ++__i)
       {
           const double rc{distro.operator()(rdev1)};
           test_2nd_Ctor.m_data[__i] = rc;
       }
       
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_data to file=%s\n",fname);
       for(std::size_t __i{0}; __i != test_2nd_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f]\n", test_2nd_Ctor.m_data[__i]);
                   
       }
     
      fclose(fp);
      printf("[UNIT-TEST #1]: darray_r8_t_2nd_Ctor() -- Ended correctly!!\n");
}

void unit_test_darray_r8_t_3rd_Ctor();

void unit_test_darray_r8_t_3rd_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{14ull*4096ull};
      
       constexpr int32_t prot{PROT_READ | PROT_WRITE};
       constexpr int32_t flags{MAP_ANONYMOUS | MAP_PRIVATE};
       constexpr int32_t fd{-1};
       constexpr off_t offset{0};
       constexpr int32_t fsize{0};
       const char * fname{"UNIT_TEST_Output_Random_Data_darray_r8_t_3rd_Ctor.csv"};
       std::clock_t seed{0ULL};
       std::uniform_real_distribution<double> distro{};
       FILE * fp{NULL};
       int32_t status;
       darray_r8_t test_3rd_Ctor = darray_r8_t(nx,prot,flags,fd,offset,fsize);
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
       seed   = std::clock();
       distro = std::uniform_real_distribution<double>(0.0,1.0); 
       auto rdev1{std::mt19937(seed)};
       for(std::size_t __i{0ULL}; __i != test_3rd_Ctor.mnx; ++__i)
       {
           double rc{distro.operator()(rdev1)};
           test_3rd_Ctor.m_data[__i] = rc; 
                                      
       }
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_data to file=%s \n",fname);
       for(std::size_t __i{0}; __i != test_3rd_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f]\n", test_3rd_Ctor.m_data[__i]);
                                        
       }
       
       fclose(fp);
       printf("[UNIT-TEST #2]: darray_r8_t_3rd_Ctor() -- Ended correctly!!\n");
}


void unit_test_darray_r8_t_4th_Ctor();

void unit_test_darray_r8_t_4th_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{11ull*4096ull};
      
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_darray_r8_t_4th_Ctor.csv"};
       std::clock_t seed{0ULL};
       std::uniform_real_distribution<double> distro{};
       FILE * fp{NULL};
       std::vector<double> data = std::vector<double>(nx);
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");                                
       seed = std::clock();
       distro = std::uniform_real_distribution<double>(0.0,1.0); 
       auto rdev1{std::mt19937(seed)};
       for(std::size_t __i{0ULL}; __i != nx; ++__i)
       {
           const double rc{distro.operator()(rdev1)};
           data[__i] = rc; 
                                        
       }
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       darray_r8_t test_4th_Ctor = darray_r8_t(data);
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
           fprintf(fp,"[%.7f]\n", test_4th_Ctor.m_data[__i]);
                                        
       }
       
       fclose(fp);
       printf("[UNIT-TEST #3]: darray_r8_t_4th_Ctor() -- Ended correctly!!\n");
}


void unit_test_darray_r8_t_5th_Ctor();

void unit_test_darray_r8_t_5th_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{11ull*4096ull};
      
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_darray_r8_t_5th_Ctor.csv"};
       std::clock_t seed{0ULL};
       std::uniform_real_distribution<double> distro{};
       FILE * fp{NULL};
       std::valarray<double> data = std::valarray<double>(nx);
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");                                
       seed   = std::clock();
       distro = std::uniform_real_distribution<double>(0.0,1.0); 
       auto rdev1{std::mt19937(seed)};
       for(std::size_t __i{0ULL}; __i != nx; ++__i)
       {
           const double rc{distro.operator()(rdev1)};
           data[__i] = rc;
       }
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       darray_r8_t test_5th_Ctor = darray_r8_t(data);
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
           fprintf(fp,"[%.7f]\n", test_5th_Ctor.m_data[__i]);
                                       
       }
       
       fclose(fp);
       printf("[UNIT-TEST #4]: darray_r8_t_5th_Ctor() -- Ended correctly!!\n");
}

void unit_test_darray_r8_t_move_assign();

void unit_test_darray_r8_t_move_assign()
{
       using namespace gms;
       constexpr std::size_t nx{11ull*4096ull};
      
       int32_t status;
       const char * fname1{"UNIT_TEST_Output_Random_Data_1_darray_r8_t_move_assign.csv"};
       const char * fname2{"UNIT_TEST_Output_Random_Data_2_darray_r8_t_move_assign.csv"};
       std::clock_t seed{0ULL};
       std::uniform_real_distribution<double> distro{};
       FILE * fp1{NULL};
       FILE * fp2{NULL};
       std::valarray<double> data1 = std::valarray<double>(nx);
       std::valarray<double> data2 = std::valarray<double>(nx);
       printf("[UNIT-TEST #5]: -- Start of random data generation and initialization.\n");                                
       seed   = std::clock();
       distro = std::uniform_real_distribution<double>(0.0,1.0); 
       auto rdev1{std::mt19937(seed)};
       for(std::size_t __i{0ULL}; __i != nx; ++__i)
       {
           const double rc1{distro.operator()(rdev1)};
           data1[__i] = rc1;
           const double rc2{distro.operator()(rdev1)};
           data2[__i] = rc2;
       }
       printf("[UNIT-TEST #5]: End of random data generation and initialization.\n");
       darray_r8_t test_5th_Ctor_1 = darray_r8_t(data1);
       darray_r8_t test_5th_Ctor_2 = darray_r8_t(data2);
       char * ctor_name{gms::common::demangle(typeid(test_5th_Ctor_1).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST #5]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST #5]: Instantiation of object Constructor of type: %s\n\n", typeid(test_5th_Ctor_1).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       test_5th_Ctor_1.info_size_alignment();
       test_5th_Ctor_2.info_size_alignment();
       fp1 = fopen(fname1,"w+");
       if(!fp1)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       fp2 = fopen(fname2,"w+");
       if(!fp2)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       // test move assignment 
       test_5th_Ctor_2 = std::move(test_5th_Ctor_1);
       printf("[UNIT-TEST #5]: Dumping: m_data to file=%s, file=%s  \n",fname1,fname2);
       for(std::size_t __i{0}; __i != test_5th_Ctor_1.mnx; ++__i)
       {
           fprintf(fp1,"[%.7f]\n", test_5th_Ctor_1.m_data[__i]);
           
                                       
       }
       for(std::size_t __i{0}; __i != test_5th_Ctor_2.mnx; ++__i)
       {
           
           fprintf(fp2,"[%.7f]\n", test_5th_Ctor_2.m_data[__i]);
                                       
       }

       
       fclose(fp1);
       fclose(fp2);
       printf("[UNIT-TEST #5]: darray_r8_t_move_assign() -- Ended correctly!!\n");
}



int main()
{
   unit_test_darray_r8_t_2nd_Ctor();
   unit_test_darray_r8_t_3rd_Ctor();
   unit_test_darray_r8_t_4th_Ctor();
   unit_test_darray_r8_t_5th_Ctor();
   unit_test_darray_r8_t_move_assign();
   return 0;
}
