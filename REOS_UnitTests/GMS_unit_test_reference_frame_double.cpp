
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_reference_frame_double.h"

static const char *fun_names[] = {"Args_7_ReferenceFrame_double_t_Ctor", "Args_11_ReferenceFrame_double_t_Ctor"};

void unit_test_Args_7_ReferenceFrame_double_t_Ctor();

void unit_test_Args_7_ReferenceFrame_double_t_Ctor()
{
     using namespace gms::fdm;
     constexpr std::size_t nx{4096ULL};
     constexpr std::size_t ny{4096ULL};
     constexpr std::size_t nz{4096ULL};
     constexpr double orig_x{0.0};
     constexpr double orig_y{0.0};
     constexpr double orig_z{0.0};
     constexpr double dt{0.01};
     const char * fname{"UNIT_TEST_Output_Random_Data_Args_7_ReferenceFrame_double_t_Ctor.csv"};
     std::clock_t seed_FI{0ULL};
     std::clock_t seed_dFI{0ULL};
     std::clock_t seed_ddFI{0ULL};
     FILE * fp = NULL; 
     printf("[Unit-Test]: Started test of:  %s()\n", fun_names[0]);
     ReferenceFrame_double_t testReferenceFrame = ReferenceFrame_double_t(nx,ny,nz,orig_x,orig_y,orig_z,dt);
     printf("[Unit-Test]: Instantiation of object Constructor of type: %s\n\n", typeid(testReferenceFrame).name());
     testReferenceFrame.info_size_and_alignment();
     fp = fopen(fname,"w+");
     if(!fp)
     {
        std::perror("fopen failed to open a file -- TERMINATING!!");
        std::exit(EXIT_FAILURE);
     }
     printf("[Unit-Test]: -- Start of random data generation and initialization.\n");
     seed_FI   = std::clock();
     auto srand_FI{std::bind(std::uniform_real_distribution<double>{0.0,static_cast<float>(testReferenceFrame.mnx)},
                                                                  std::mt19937(seed_FI))};
     seed_dFI  = std::clock();
     auto srand_dFI{std::bind(std::uniform_real_distribution<double>{0.0,static_cast<float>(testReferenceFrame.mny)},
                                                                   std::mt19937(seed_dFI))};
     seed_ddFI = std::clock();
     auto srand_ddFI{std::bind(std::uniform_real_distribution<double>{0.0,static_cast<float>(testReferenceFrame.mnz)},
                                                                   std::mt19937(seed_ddFI))};
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
          const double rFI{srand_FI()};
          testReferenceFrame.mFI_x[__i] = rFI;
          testReferenceFrame.mFI_y[__i] = rFI;
          testReferenceFrame.mFI_z[__i] = rFI;
     }
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
          const double rdFI{srand_dFI()};
          testReferenceFrame.mdFI_x[__i] = rdFI;
          testReferenceFrame.mdFI_y[__i] = rdFI;
          testReferenceFrame.mdFI_z[__i] = rdFI;
     }
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
          const double rddFI{srand_ddFI()};
          testReferenceFrame.mddFI_x[__i] = rddFI;
          testReferenceFrame.mddFI_y[__i] = rddFI;
          testReferenceFrame.mddFI_z[__i] = rddFI;
     }
     printf("[Unit-Test]: End of random data generation and initialization.\n");
     printf("[Unit-Test]: Dumping: FI_x \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
        fprintf(fp, "%.16f\n", testReferenceFrame.mFI_x[__i]);
     }
     printf("[Unit-Test]: Dumping: FI_y \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
         fprintf(fp, "%.16f\n", testReferenceFrame.mFI_y[__i]);
     }
     printf("Unit-Test]: Dumping: FI_z \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
          fprintf(fp, "%.16f\n", testReferenceFrame.mFI_z[__i]); 
     }
     printf("[Unit-Test]: Dumping: dFI_x \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
        fprintf(fp, "%.16f\n", testReferenceFrame.mdFI_x[__i]);
     }
     printf("[Unit-Test]: Dumping: dFI_y \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
         fprintf(fp, "%.16f\n", testReferenceFrame.mdFI_y[__i]);
     }
     printf("Unit-Test]: Dumping: dFI_z \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
         fprintf(fp, "%.16f\n", testReferenceFrame.mdFI_z[__i]); 
     }
     printf("[Unit-Test]: Dumping: ddFI_x \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
        fprintf(fp, "%.16f\n", testReferenceFrame.mddFI_x[__i]);
     }
     printf("[Unit-Test]: Dumping: ddFI_y \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
         fprintf(fp, "%.16f\n", testReferenceFrame.mddFI_y[__i]);
     }
     printf("Unit-Test]: Dumping: ddFI_z \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
          fprintf(fp, "%.16f\n", testReferenceFrame.mddFI_z[__i]); 
     }
     printf("[Unit-Test]: Reached Normal End of Test: %s\n", fun_names[0]);
     fclose(fp);
}

void unit_test_Args_11_ReferenceFrame_double_t_Ctor();

void unit_test_Args_11_ReferenceFrame_double_t_Ctor()
{
     using namespace gms::fdm;
     constexpr std::size_t nx{4096ULL};
     constexpr std::size_t ny{4096ULL};
     constexpr std::size_t nz{4096ULL};
     constexpr double orig_x{0.0};
     constexpr double orig_y{0.0};
     constexpr double orig_z{0.0};
     constexpr double dt{0.01};
     constexpr int32_t prot{PROT_READ | PROT_WRITE};
     constexpr int32_t flags{MAP_ANONYMOUS | MAP_PRIVATE};
     constexpr int32_t fd{-1};
     constexpr off_t offset{0};
     constexpr int32_t fsize{0};
     const char * fname{"UNIT_TEST_Output_Random_Data_Args_11_ReferenceFrame_double_t_Ctor.csv"};
     std::clock_t seed_FI{0ULL};
     std::clock_t seed_dFI{0ULL};
     std::clock_t seed_ddFI{0ULL};
     FILE * fp = NULL; 
     printf("[Unit-Test]: Started test of:  %s()\n", fun_names[1]);
     ReferenceFrame_double_t testReferenceFrame = ReferenceFrame_double_t(nx,ny,nz,orig_x,orig_y,orig_z,dt,
                                                                          prot,flags,fd,offset,fsize);
     printf("[Unit-Test]: Instantiation of object Constructor of type: %s\n\n", typeid(testReferenceFrame).name());
     testReferenceFrame.info_size_and_alignment();
     fp = fopen(fname,"w+");
     if(!fp)
     {
        std::perror("fopen failed to open a file -- TERMINATING!!");
        std::exit(EXIT_FAILURE);
     }
     printf("[Unit-Test]: -- Start of random data generation and initialization.\n");
     seed_FI   = std::clock();
     auto srand_FI{std::bind(std::uniform_real_distribution<double>{0.0,static_cast<float>(testReferenceFrame.mnx)},
                                                                  std::mt19937(seed_FI))};
     seed_dFI  = std::clock();
     auto srand_dFI{std::bind(std::uniform_real_distribution<double>{0.0,static_cast<float>(testReferenceFrame.mny)},
                                                                   std::mt19937(seed_dFI))};
     seed_ddFI = std::clock();
     auto srand_ddFI{std::bind(std::uniform_real_distribution<double>{0.0,static_cast<float>(testReferenceFrame.mnz)},
                                                                   std::mt19937(seed_ddFI))};
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
          const double rFI{srand_FI()};
          testReferenceFrame.mFI_x[__i] = rFI;
          testReferenceFrame.mFI_y[__i] = rFI;
          testReferenceFrame.mFI_z[__i] = rFI;
     }
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
          const double rdFI{srand_dFI()};
          testReferenceFrame.mdFI_x[__i] = rdFI;
          testReferenceFrame.mdFI_y[__i] = rdFI;
          testReferenceFrame.mdFI_z[__i] = rdFI;
     }
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
          const double rddFI{srand_ddFI()};
          testReferenceFrame.mddFI_x[__i] = rddFI;
          testReferenceFrame.mddFI_y[__i] = rddFI;
          testReferenceFrame.mddFI_z[__i] = rddFI;
     }
     printf("[Unit-Test]: End of random data generation and initialization.\n");
     printf("[Unit-Test]: Dumping: FI_x \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
        fprintf(fp, "%.16f\n", testReferenceFrame.mFI_x[__i]);
     }
     printf("[Unit-Test]: Dumping: FI_y \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
         fprintf(fp, "%.16f\n", testReferenceFrame.mFI_y[__i]);
     }
     printf("Unit-Test]: Dumping: FI_z \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
          fprintf(fp, "%.16f\n", testReferenceFrame.mFI_z[__i]); 
     }
     printf("[Unit-Test]: Dumping: dFI_x \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
        fprintf(fp, "%.16f\n", testReferenceFrame.mdFI_x[__i]);
     }
     printf("[Unit-Test]: Dumping: dFI_y \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
         fprintf(fp, "%.16f\n", testReferenceFrame.mdFI_y[__i]);
     }
     printf("Unit-Test]: Dumping: dFI_z \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
         fprintf(fp, "%.16f\n", testReferenceFrame.mdFI_z[__i]); 
     }
     printf("[Unit-Test]: Dumping: ddFI_x \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
        fprintf(fp, "%.16f\n", testReferenceFrame.mddFI_x[__i]);
     }
     printf("[Unit-Test]: Dumping: ddFI_y \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
         fprintf(fp, "%.16f\n", testReferenceFrame.mddFI_y[__i]);
     }
     printf("Unit-Test]: Dumping: ddFI_z \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
          fprintf(fp, "%.16f\n", testReferenceFrame.mddFI_z[__i]); 
     }
     printf("[Unit-Test]: Reached Normal End of Test: %s\n", fun_names[0]);
     fclose(fp);                                                                   
     
}


int main()
{
    unit_test_Args_7_ReferenceFrame_double_t_Ctor();
    unit_test_Args_11_ReferenceFrame_double_t_Ctor();
    return 0;
}


