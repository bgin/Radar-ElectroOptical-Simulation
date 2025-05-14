
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_containers.hpp"

static const char *ctor_names[] = {"DC3D_c4_t_2nd_Ctor,DC3D_c4_t_3rd_Ctor,DC3D_c4_t_4th_Ctor,DC3D_c4_t_5th_Ctor"};

void unit_test_DC3D_c4_t_2nd_Ctor();

void unit_test_DC3D_c4_t_2nd_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{4096};
       constexpr std::size_t ny{1024};
       constexpr std::size_t nz{512};
       const char * fname{"UNIT_TEST_Output_Random_Data_DC3D_c4_t_2nd_Ctor.csv"};
       std::clock_t seeds[6] = {0ULL};
       std::uniform_real_distribution<float> distros[6]{};
       FILE * fp{NULL};
       printf("UNIT-TEST]: Started test of: %s()\n", ctor_names[0]);
       DC3D_c4_t test_2nd_Ctor = DC3D_c4_t(nx,ny,nz);
       printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(test_2nd_Ctor).name());
       test_2nd_Ctor.info_size_and_alignment();
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
       for(std::size_t __i{0ULL}; __i != test_2nd_Ctor.mnx; ++__i)
       {
           const std::complex<float> rc{distros[0].operator()(rdev1),distros[1].operator()(rdev2)};
           test_2nd_Ctor.m_Ex[__i] = rc;
       }
       auto rdev3{std::mt19937(seeds[2])};
       auto rdev4{std::mt19937(seeds[3])};
       for(std::size_t __i{0ULL}; __i != test_2nd_Ctor.mny; ++__i)
       {
           const std::complex<float> rc{distros[2].operator()(rdev3),distros[3].operator()(rdev4)};
           test_2nd_Ctor.m_Ey[__i] = rc;
       }
       auto rdev5{std::mt19937(seeds[4])};
       auto rdev6{std::mt19937(seeds[5])}
       for(std::size_t __i{0ULL}; __i != test_2nd_Ctor.mnz; ++__i)
       {
           const std::complex<float> rc{distros[4].operator()(rdev5),distros[5].operator()(rdev6)};
           test_2nd_Ctor.m_Ez[__i] = rc;
       }
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_Ex \n");
       for(std::size_t __i{0}; __i != test_2nd_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_2nd_Ctor.m_Ex[__i].real(),
                                        test_2nd_Ctor.m_Ex[__i].imag());
       }
      printf("[UNIT-TEST]: Dumping: m_Ey \n");
      for(std::size_t __i{0}; __i != test_2nd_Ctor.mny; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", test_2nd_Ctor.m_Ey[__i].real(),
                                        test_2nd_Ctor.m_Ey[__i].imag());
      }
      printf("[UNIT-TEST]: Dumping: m_Ez \n");
      for(std::size_t __i{0}; __i != test_2nd_Ctor.mnz; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", test_2nd_Ctor.m_Ez[__i].real(),
                                        test_2nd_Ctor.m_Ez[__i].imag());
      }
      fclose(fp);
}










int main()
{   
    unit_test_DC3D_c4_t_2nd_Ctor();
    return 0;
}