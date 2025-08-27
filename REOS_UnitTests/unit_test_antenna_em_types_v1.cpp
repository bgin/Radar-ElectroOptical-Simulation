#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <complex>
#include "GMS_malloc.h"
#include "GMS_antenna_em_types_v1.h"

/*
   icpc -o unit_test_antenna_em_types_v1 -fp-model -std=c++17 fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_antenna_em_types_v1.h.h unit_test_antenna_em_types_v1.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32   GMS_config.h GMS_malloc.h GMS_antenna_em_types_v1.h.h unit_test_antenna_em_types_v1.cpp

*/

void unit_test_Exyz_t_ctors();

void unit_test_Exyz_t_ctors()
{
       using namespace gms::radiolocation;
       
       constexpr std::size_t nx{4096};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_Exyz_t_Ctor.csv"};
       std::clock_t seeds[6] = {0ULL};
       std::uniform_real_distribution<float> distros[6]{};
       FILE * fp{NULL};
       Exyz_t<std::complex<float>> E_field = Exyz_t<std::complex<float>>(nx);
       char * ctor_name{gms::common::demangle(typeid(E_field).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(E_field).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");
       for(int __i{0}; __i != 6; ++__i)
       {
           seeds[__i] = std::clock();
           distros[__i]  = std::uniform_real_distribution<float>(0.0f,1.0f);
       }
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       auto rdev3{std::mt19937(seeds[2])};
       auto rdev4{std::mt19937(seeds[3])};
       auto rdev5{std::mt19937(seeds[4])};
       auto rdev6{std::mt19937(seeds[5])};
       for(std::size_t __i{0ULL}; __i != E_field.m_npts; ++__i)
       {
           const std::complex<float> rc1{distros[0].operator()(rdev1),distros[1].operator()(rdev2)};
           E_field.m_ex[__i] = rc1;
           const std::complex<float> rc2{distros[2].operator()(rdev3),distros[3].operator()(rdev4)};
           E_field.m_ey[__i] = rc2;
           const std::complex<float> r3{distros[4].operator()(rdev5),distros[5].operator()(rdev6)};
           E_field.m_ez[__i] = r3;
       }
       
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_ex \n");
       for(std::size_t __i{0}; __i != E_field.m_npts; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", E_field.m_ex[__i].real(),
                                        E_field.m_ex[__i].imag());
       }
      printf("[UNIT-TEST]: Dumping: m_ey \n");
      for(std::size_t __i{0}; __i != E_field.m_npts; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", E_field.m_ey[__i].real(),
                                        E_field.m_ey[__i].imag());
      }
      printf("[UNIT-TEST]: Dumping: m_Ez \n");
      for(std::size_t __i{0}; __i != E_field.m_npts; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", E_field.m_ez[__i].real(),
                                        E_field.m_ez[__i].imag());
      }
      fclose(fp);
}


void unit_test_Exyz_t_operators();

void unit_test_Exyz_t_operators()
{
     using namespace gms::radiolocation;
       constexpr std::size_t nx{4096};
       constexpr std::size_t nx2{8192};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_Exyz_t_operators.csv"};
       std::clock_t seeds[6] = {0ULL};
       std::uniform_real_distribution<float> distros[6]{};
       FILE * fp{NULL};
       Exyz_t<std::complex<float>> E_field = Exyz_t<std::complex<float>>(nx);
       Exyz_t<std::complex<float>> E_field_cp = Exyz_t<std::complex<float>>(nx2);
       char * ctor_name{gms::common::demangle(typeid(E_field).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(E_field).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");
       for(int __i{0}; __i != 6; ++__i)
       {
           seeds[__i] = std::clock();
           distros[__i]  = std::uniform_real_distribution<float>(0.0f,1.0f);
       }
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       auto rdev3{std::mt19937(seeds[2])};
       auto rdev4{std::mt19937(seeds[3])};
       auto rdev5{std::mt19937(seeds[4])};
       auto rdev6{std::mt19937(seeds[5])};
       for(std::size_t __i{0ULL}; __i != E_field.m_npts; ++__i)
       {
           const std::complex<float> rc1{distros[0].operator()(rdev1),distros[1].operator()(rdev2)};
           E_field.m_ex[__i] = rc1;
           const std::complex<float> rc2{distros[2].operator()(rdev3),distros[3].operator()(rdev4)};
           E_field.m_ey[__i] = rc2;
           const std::complex<float> r3{distros[4].operator()(rdev5),distros[5].operator()(rdev6)};
           E_field.m_ez[__i] = r3;
       }

       for(std::size_t __i{0ULL}; __i != E_field_cp.m_npts; ++__i)
       {
           const std::complex<float> rc1{distros[0].operator()(rdev1),distros[1].operator()(rdev2)};
           E_field_cp.m_ex[__i] = rc1;
           const std::complex<float> rc2{distros[2].operator()(rdev3),distros[3].operator()(rdev4)};
           E_field_cp.m_ey[__i] = rc2;
           const std::complex<float> r3{distros[4].operator()(rdev5),distros[5].operator()(rdev6)};
           E_field_cp.m_ez[__i] = r3;
       }

       E_field_cp.operator=(std::move(E_field));
       
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_ex \n");
       for(std::size_t __i{0}; __i != E_field_cp.m_npts; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", E_field_cp.m_ex[__i].real(),
                                        E_field_cp.m_ex[__i].imag());
       }
      printf("[UNIT-TEST]: Dumping: m_ey \n");
      for(std::size_t __i{0}; __i != E_field_cp.m_npts; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", E_field_cp.m_ey[__i].real(),
                                        E_field_cp.m_ey[__i].imag());
      }
      printf("[UNIT-TEST]: Dumping: m_Ez \n");
      for(std::size_t __i{0}; __i != E_field_cp.m_npts; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", E_field_cp.m_ez[__i].real(),
                                        E_field_cp.m_ez[__i].imag());
      }
      fclose(fp);
}


void unit_test_Hxyz_t_ctors();

void unit_test_Hxyz_t_ctors()
{
       using namespace gms::radiolocation;
       
       constexpr std::size_t nx{4096};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_Hxyz_t_Ctor.csv"};
       std::clock_t seeds[6] = {0ULL};
       std::uniform_real_distribution<float> distros[6]{};
       FILE * fp{NULL};
       Hxyz_t<std::complex<float>> H_field = Hxyz_t<std::complex<float>>(nx);
       char * ctor_name{gms::common::demangle(typeid(H_field).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(H_field).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");
       for(int __i{0}; __i != 6; ++__i)
       {
           seeds[__i] = std::clock();
           distros[__i]  = std::uniform_real_distribution<float>(0.0f,1.0f);
       }
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       auto rdev3{std::mt19937(seeds[2])};
       auto rdev4{std::mt19937(seeds[3])};
       auto rdev5{std::mt19937(seeds[4])};
       auto rdev6{std::mt19937(seeds[5])};
       for(std::size_t __i{0ULL}; __i != H_field.m_npts; ++__i)
       {
           const std::complex<float> rc1{distros[0].operator()(rdev1),distros[1].operator()(rdev2)};
           H_field.m_hx[__i] = rc1;
           const std::complex<float> rc2{distros[2].operator()(rdev3),distros[3].operator()(rdev4)};
           H_field.m_hy[__i] = rc2;
           const std::complex<float> r3{distros[4].operator()(rdev5),distros[5].operator()(rdev6)};
           H_field.m_hz[__i] = r3;
       }
       
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_hx \n");
       for(std::size_t __i{0}; __i != H_field.m_npts; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", H_field.m_hx[__i].real(),
                                        H_field.m_hx[__i].imag());
       }
      printf("[UNIT-TEST]: Dumping: m_hy \n");
      for(std::size_t __i{0}; __i != H_field.m_npts; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", H_field.m_hy[__i].real(),
                                        H_field.m_hy[__i].imag());
      }
      printf("[UNIT-TEST]: Dumping: m_hz \n");
      for(std::size_t __i{0}; __i != H_field.m_npts; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", H_field.m_hz[__i].real(),
                                        H_field.m_hz[__i].imag());
      }
      fclose(fp);
}


void unit_test_Hxyz_t_operators();

void unit_test_Hxyz_t_operators()
{
     using namespace gms::radiolocation;
       constexpr std::size_t nx{4096};
       constexpr std::size_t nx2{8192};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_Hxyz_t_operators.csv"};
       std::clock_t seeds[6] = {0ULL};
       std::uniform_real_distribution<float> distros[6]{};
       FILE * fp{NULL};
       Hxyz_t<std::complex<float>> H_field = Hxyz_t<std::complex<float>>(nx);
       Hxyz_t<std::complex<float>> H_field_cp = Hxyz_t<std::complex<float>>(nx2);
       char * ctor_name{gms::common::demangle(typeid(H_field).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(H_field).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");
       for(int __i{0}; __i != 6; ++__i)
       {
           seeds[__i] = std::clock();
           distros[__i]  = std::uniform_real_distribution<float>(0.0f,1.0f);
       }
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       auto rdev3{std::mt19937(seeds[2])};
       auto rdev4{std::mt19937(seeds[3])};
       auto rdev5{std::mt19937(seeds[4])};
       auto rdev6{std::mt19937(seeds[5])};
       for(std::size_t __i{0ULL}; __i != H_field.m_npts; ++__i)
       {
           const std::complex<float> rc1{distros[0].operator()(rdev1),distros[1].operator()(rdev2)};
           H_field.m_hx[__i] = rc1;
           const std::complex<float> rc2{distros[2].operator()(rdev3),distros[3].operator()(rdev4)};
           H_field.m_hy[__i] = rc2;
           const std::complex<float> r3{distros[4].operator()(rdev5),distros[5].operator()(rdev6)};
           H_field.m_hz[__i] = r3;
       }

       for(std::size_t __i{0ULL}; __i != H_field_cp.m_npts; ++__i)
       {
           const std::complex<float> rc1{distros[0].operator()(rdev1),distros[1].operator()(rdev2)};
           H_field_cp.m_hx[__i] = rc1;
           const std::complex<float> rc2{distros[2].operator()(rdev3),distros[3].operator()(rdev4)};
           H_field_cp.m_hy[__i] = rc2;
           const std::complex<float> r3{distros[4].operator()(rdev5),distros[5].operator()(rdev6)};
           H_field_cp.m_hz[__i] = r3;
       }

       H_field_cp.operator=(std::move(H_field));
       
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_hx \n");
       for(std::size_t __i{0}; __i != H_field_cp.m_npts; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", H_field_cp.m_hx[__i].real(),
                                        H_field_cp.m_hx[__i].imag());
       }
      printf("[UNIT-TEST]: Dumping: m_hy \n");
      for(std::size_t __i{0}; __i != H_field_cp.m_npts; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", H_field_cp.m_hy[__i].real(),
                                        H_field_cp.m_hy[__i].imag());
      }
      printf("[UNIT-TEST]: Dumping: m_hz \n");
      for(std::size_t __i{0}; __i != H_field_cp.m_npts; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", H_field_cp.m_hz[__i].real(),
                                        H_field_cp.m_hz[__i].imag());
      }
      fclose(fp);
}



int main()
{
    unit_test_Exyz_t_ctors();
    unit_test_Exyz_t_operators();
    unit_test_Hxyz_t_ctors();
    unit_test_Hxyz_t_operators();
    return 0;
}