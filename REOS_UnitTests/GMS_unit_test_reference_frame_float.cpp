
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_reference_frame_float.h"

static const char *fun_names[] = {"Args_7_ReferenceFrame_float_t", "ARGS_11_ReferenceFrame_float_t",
                                  "ReferenceFrame_float_t(ReferenceFrame_float_t &&)", "operator=(ReferenceFrame_float_t &&)"};


void unit_test_Args_7_ReferenceFrame_float_t_Ctor()
{
    using namespace gms::fdm;
     constexpr std::size_t nx{4096ULL};
     constexpr std::size_t ny{4096ULL};
     constexpr std::size_t nz{4096ULL};
     constexpr float orig_x{0.0f};
     constexpr float orig_y{0.0};
     constexpr float orig_z{0.0};
     constexpr float dt{0.01f};
     const char * fname{"UNIT_TEST_Output_Random_Data_Args_7_ReferenceFrame_float_t_Ctor.csv"};
     std::clock_t seed_FI{0ULL};
     std::clock_t seed_dFI{0ULL};
     std::clock_t seed_ddFI{0ULL};
     FILE * fp = NULL; 
     printf("[Unit-Test]: Start %s()\n", fun_names[0]);
     ReferenceFrame_float_t testReferenceFrame = ReferenceFrame_float_t(nx,ny,nz,orig_x,orig_y,orig_z,dt);
     printf("[Unit-Test]: Instantiation of object Constructor of type: %s\n", typeid(testReferenceFrame).name());
     printf("[Unit-Test]: End %s()\n", fun_names[0]);
     fp = fopen(fname,"w+");
     if(!fp)
     {
        std::perror("fopen failed to open a file -- TERMINATING!!");
        std::exit(EXIT_FAILURE);
     }
     printf("[Unit-Test]: -- Start of random data generation and initialization.\n");
     seed_FI   = std::clock();
     auto srand_FI{std::bind(std::uniform_real_distribution<float>{0.0f,static_cast<float>(testReferenceFrame.mnx)},
                                                                  std::mt19937(seed_FI))};
     seed_dFI  = std::clock();
     auto srand_dFI{std::bind(std::uniform_real_distribution<float>{0.0f,static_cast<float>(testReferenceFrame.mny)},
                                                                   std::mt19937(seed_dFI))};
     seed_ddFI = std::clock();
     auto srand_ddFI{std::bind(std::uniform_real_distribution<float>{0.0f,static_cast<float>(testReferenceFrame.mnz)},
                                                                   std::mt19937(seed_ddFI))};
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
          const float rFI{srand_FI()};
          testReferenceFrame.mFI_x[__i] = rFI;
          testReferenceFrame.mFI_y[__i] = rFI;
          testReferenceFrame.mFI_z[__i] = rFI;
     }
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
          const float rdFI{srand_dFI()};
          testReferenceFrame.mdFI_x[__i] = rdFI;
          testReferenceFrame.mdFI_y[__i] = rdFI;
          testReferenceFrame.mdFI_z[__i] = rdFI;
     }
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
          const float rddFI{srand_ddFI()};
          testReferenceFrame.mddFI_x[__i] = rddFI;
          testReferenceFrame.mddFI_y[__i] = rddFI;
          testReferenceFrame.mddFI_z[__i] = rddFI;
     }
     printf("[Unit-Test]: End of random data generation and initialization.\n");
     printf("[Unit-Test]: Dumping: FI_x \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
        fprintf(fp, "%.7f\n", testReferenceFrame.mFI_x[__i]);
     }
     printf("[Unit-Test]: Dumping: FI_y \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
         fprintf(fp, "%.7f\n", testReferenceFrame.mFI_y[__i]);
     }
     printf("Unit-Test]: Dumping: FI_z \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
          fprintf(fp, "%.7f\n", testReferenceFrame.mFI_z[__i]); 
     }
     printf("[Unit-Test]: Dumping: dFI_x \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
        fprintf(fp, "%.7f\n", testReferenceFrame.mdFI_x[__i]);
     }
     printf("[Unit-Test]: Dumping: dFI_y \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
         fprintf(fp, "%.7f\n", testReferenceFrame.mdFI_y[__i]);
     }
     printf("Unit-Test]: Dumping: dFI_z \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
         fprintf(fp, "%.7f\n", testReferenceFrame.mdFI_z[__i]); 
     }
     printf("[Unit-Test]: Dumping: ddFI_x \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnx; ++__i)
     {
        fprintf(fp, "%.7f\n", testReferenceFrame.mddFI_x[__i]);
     }
     printf("[Unit-Test]: Dumping: ddFI_y \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mny; ++__i)
     {
         fprintf(fp, "%.7f\n", testReferenceFrame.mddFI_y[__i]);
     }
     printf("Unit-Test]: Dumping: ddFI_z \n");
     for(std::size_t __i = 0ULL; __i != testReferenceFrame.mnz; ++__i)
     {
          fprintf(fp, "%.7f\n", testReferenceFrame.mddFI_z[__i]); 
     }
     printf("[Unit-Test]: Reached Normal End of Test: %s\n", fun_names[0]);
     fclose(fp);
}


int main()
{
    unit_test_Args_7_ReferenceFrame_float_t_Ctor();
    return 0;
}