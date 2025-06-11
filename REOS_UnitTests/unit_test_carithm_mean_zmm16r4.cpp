

#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_malloc.h" // for demangle
#include "GMS_carithm_mean_zmm16r4.h"
/*
icpc -o unit_test_carithm_mean_zmm16r4 -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_carithm_mean_zmm16r4.h GMS_carithm_mean_zmm16r4.cpp unit_test_carithm_mean_zmm16r4.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_carithm_mean_zmm16r4.h GMS_carithm_mean_zmm16r4.cpp unit_test_carithm_mean_zmm16r4.cpp 
                        
*/
__attribute__((hot))
void unit_test_carithm_mean_zmm16r4_u();

void unit_test_carithm_mean_zmm16r4_u()
{
     
     constexpr std::size_t nelems{16384ull};
     float buf_re[nelems];
     float buf_im[nelems];
     std::clock_t seedre{0ull};
     std::clock_t seedim{0ull};
     std::clock_t seed{};
     float mean_re{};
     float mean_im{};
     int32_t      which{-1};
     seed = std::clock();
     auto rand{std::bind(std::uniform_int_distribution<int32_t>{0,4},
                          std::mt19937(seed))};
     which = rand();
     switch(which)
     {
         case 0 : 
         {
            std::uniform_real_distribution<float> uni_distro{};
            char * __restrict__ distro_name{NULL};
            int32_t status{};
            distro_name = gms::common::demangle(typeid(uni_distro).name(),status);
            if(status==0 && distro_name != NULL)
            {
               printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", distro_name);
               gms::common::gms_mm_free(distro_name);
            }
            else
            {
              printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(uni_distro).name());
              if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            printf("UNIT-TEST]: -- Start of data initialization\n");
            seedre = std::clock();
            seedim = std::clock();
            auto srand_re{std::bind(std::uniform_real_distribution<float>{0.0f,1.0f},
                                 std::mt19937(seedre))};
            auto srand_im{std::bind(std::uniform_real_distribution<float>{0.0f,1.0f},
                                 std::mt19937(seedim))};
            std::generate((float*)&buf_re[0],(float*)&buf_re[0]+nelems,srand_re);
            std::generate((float*)&buf_im[0],(float*)&buf_im[0]+nelems,srand_im);
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cmean_arithm_u10x_zmm16r4_u\n");
            gms::math::cmean_arithm_u10x_zmm16r4_u(&buf_re[0],&buf_im[0],
                                                   &mean_re,&mean_im,nelems);
            printf("[UNIT-TEST]: -- END: cmean_arithm_u10x_zmm16r4_u\n");
            printf("mean_re=%.7f, mean_im=%.7f\n",mean_re,mean_im);
         }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro{};
            char * __restrict__ distro_name{NULL};
            int32_t status{};
            distro_name = gms::common::demangle(typeid(normal_distro).name(),status);
            if(status==0 && distro_name != NULL)
            {
               printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", distro_name);
               gms::common::gms_mm_free(distro_name);
            }
            else
            {
              printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(normal_distro).name());
              if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            printf("UNIT-TEST]: -- Start of data initialization\n");
            seedre = std::clock();
            seedim = std::clock();
            auto srand_re{std::bind(std::normal_distribution<float>{0.0f,1.0f},
                                 std::mt19937(seedre))};
            auto srand_im{std::bind(std::normal_distribution<float>{0.0f,1.0f},
                                 std::mt19937(seedim))};
            std::generate((float*)&buf_re[0],(float*)&buf_re[0]+nelems,srand_re);
            std::generate((float*)&buf_im[0],(float*)&buf_im[0]+nelems,srand_im);
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cmean_arithm_u10x_zmm16r4_u\n");
            gms::math::cmean_arithm_u10x_zmm16r4_u(&buf_re[0],&buf_im[0],
                                                   &mean_re,&mean_im,nelems);
            printf("[UNIT-TEST]: -- END: cmean_arithm_u10x_zmm16r4_u\n");
            printf("mean_re=%.7f, mean_im=%.7f\n",mean_re,mean_im);
        }
        break;
        case 2 : 
        {
            std::lognormal_distribution<float> lognormal_distro{};
            char * __restrict__ distro_name{NULL};
            int32_t status{};
            distro_name = gms::common::demangle(typeid(lognormal_distro).name(),status);
            if(status==0 && distro_name != NULL)
            {
               printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", distro_name);
               gms::common::gms_mm_free(distro_name);
            }
            else
            {
              printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognormal_distro).name());
              if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            printf("UNIT-TEST]: -- Start of data initialization\n");
            seedre = std::clock();
            seedim = std::clock();
            auto srand_re{std::bind(std::lognormal_distribution<float>{1.6f,0.25f},
                                 std::mt19937(seedre))};
            auto srand_im{std::bind(std::lognormal_distribution<float>{1.6f,0.25f},
                                 std::mt19937(seedim))};
            std::generate((float*)&buf_re[0],(float*)&buf_re[0]+nelems,srand_re);
            std::generate((float*)&buf_im[0],(float*)&buf_im[0]+nelems,srand_im);
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cmean_arithm_u10x_zmm16r4_u\n");
            gms::math::cmean_arithm_u10x_zmm16r4_u(&buf_re[0],&buf_im[0],
                                                   &mean_re,&mean_im,nelems);
            printf("[UNIT-TEST]: -- END: cmean_arithm_u10x_zmm16r4_u\n");
            printf("mean_re=%.7f, mean_im=%.7f\n",mean_re,mean_im);
        }
        break;
        case 3 : 
        {
            std::fisher_f_distribution<float> fisher_f_distro{};
            char * __restrict__ distro_name{NULL};
            int32_t status{};
            distro_name = gms::common::demangle(typeid(fisher_f_distro).name(),status);
            if(status==0 && distro_name != NULL)
            {
               printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", distro_name);
               gms::common::gms_mm_free(distro_name);
            }
            else
            {
              printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(fisher_f_distro).name());
              if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            printf("UNIT-TEST]: -- Start of data initialization\n");
            seedre = std::clock();
            seedim = std::clock();
            auto srand_re{std::bind(std::fisher_f_distribution<float>{1.0f,5.0f},
                                 std::mt19937(seedre))};
            auto srand_im{std::bind(std::fisher_f_distribution<float>{1.0f,5.0f},
                                 std::mt19937(seedim))};
            std::generate((float*)&buf_re[0],(float*)&buf_re[0]+nelems,srand_re);
            std::generate((float*)&buf_im[0],(float*)&buf_im[0]+nelems,srand_im);
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cmean_arithm_u10x_zmm16r4_u\n");
            gms::math::cmean_arithm_u10x_zmm16r4_u(&buf_re[0],&buf_im[0],
                                                   &mean_re,&mean_im,nelems);
            printf("[UNIT-TEST]: -- END: cmean_arithm_u10x_zmm16r4_u\n");
            printf("mean_re=%.7f, mean_im=%.7f\n",mean_re,mean_im);
        }
        break;
        case 4 : 
        {
            std::cauchy_distribution<float> cauchy_distro{};
            char * __restrict__ distro_name{NULL};
            int32_t status{};
            distro_name = gms::common::demangle(typeid(cauchy_distro).name(),status);
            if(status==0 && distro_name != NULL)
            {
               printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", distro_name);
               gms::common::gms_mm_free(distro_name);
            }
            else
            {
              printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(cauchy_distro).name());
              if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            printf("UNIT-TEST]: -- Start of data initialization\n");
            seedre = std::clock();
            seedim = std::clock();
            auto srand_re{std::bind(std::cauchy_distribution<float>{-1.0f,5.0f},
                                 std::mt19937(seedre))};
            auto srand_im{std::bind(std::cauchy_distribution<float>{-1.0f,5.0f},
                                 std::mt19937(seedim))};
            std::generate((float*)&buf_re[0],(float*)&buf_re[0]+nelems,srand_re);
            std::generate((float*)&buf_im[0],(float*)&buf_im[0]+nelems,srand_im);
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cmean_arithm_u10x_zmm16r4_u\n");
            gms::math::cmean_arithm_u10x_zmm16r4_u(&buf_re[0],&buf_im[0],
                                                   &mean_re,&mean_im,nelems);
            printf("[UNIT-TEST]: -- END: cmean_arithm_u10x_zmm16r4_u\n");
            printf("mean_re=%.7f, mean_im=%.7f\n",mean_re,mean_im);
        }
        break;
        default : 
            printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
            std::terminate();
     }
}

__attribute__((hot))
void unit_test_carithm_mean_zmm16r4_a();

void unit_test_carithm_mean_zmm16r4_a()
{
     constexpr std::size_t nelems{16384ull};
     __attribute__((aligned(64)))
     float buf_re[nelems];
     __attribute__((aligned(64)))
     float buf_im[nelems];
     std::clock_t seedre{0ull};
     std::clock_t seedim{0ull};
     std::clock_t seed{};
     float mean_re{};
     float mean_im{};
     int32_t      which{-1};
     seed = std::clock();
     auto rand{std::bind(std::uniform_int_distribution<int32_t>{0,4},
                          std::mt19937(seed))};
     which = rand();
     switch(which)
     {
         case 0 : 
         {
            std::uniform_real_distribution<float> uni_distro{};
            char * __restrict__ distro_name{NULL};
            int32_t status{};
            distro_name = gms::common::demangle(typeid(uni_distro).name(),status);
            if(status==0 && distro_name != NULL)
            {
               printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", distro_name);
               gms::common::gms_mm_free(distro_name);
            }
            else
            {
              printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(uni_distro).name());
              if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            printf("UNIT-TEST]: -- Start of data initialization\n");
            seedre = std::clock();
            seedim = std::clock();
            auto srand_re{std::bind(std::uniform_real_distribution<float>{0.0f,1.0f},
                                 std::mt19937(seedre))};
            auto srand_im{std::bind(std::uniform_real_distribution<float>{0.0f,1.0f},
                                 std::mt19937(seedim))};
            std::generate((float*)&buf_re[0],(float*)&buf_re[0]+nelems,srand_re);
            std::generate((float*)&buf_im[0],(float*)&buf_im[0]+nelems,srand_im);
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cmean_arithm_u10x_zmm16r4_a\n");
            gms::math::cmean_arithm_u10x_zmm16r4_a(&buf_re[0],&buf_im[0],
                                                   &mean_re,&mean_im,nelems);
            printf("[UNIT-TEST]: -- END: cmean_arithm_u10x_zmm16r4_a\n");
            printf("mean_re=%.7f, mean_im=%.7f\n",mean_re,mean_im);
         }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro{};
            char * __restrict__ distro_name{NULL};
            int32_t status{};
            distro_name = gms::common::demangle(typeid(normal_distro).name(),status);
            if(status==0 && distro_name != NULL)
            {
               printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", distro_name);
               gms::common::gms_mm_free(distro_name);
            }
            else
            {
              printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(normal_distro).name());
              if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            printf("UNIT-TEST]: -- Start of data initialization\n");
            seedre = std::clock();
            seedim = std::clock();
            auto srand_re{std::bind(std::normal_distribution<float>{0.0f,1.0f},
                                 std::mt19937(seedre))};
            auto srand_im{std::bind(std::normal_distribution<float>{0.0f,1.0f},
                                 std::mt19937(seedim))};
            std::generate((float*)&buf_re[0],(float*)&buf_re[0]+nelems,srand_re);
            std::generate((float*)&buf_im[0],(float*)&buf_im[0]+nelems,srand_im);
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cmean_arithm_u10x_zmm16r4_a\n");
            gms::math::cmean_arithm_u10x_zmm16r4_a(&buf_re[0],&buf_im[0],
                                                   &mean_re,&mean_im,nelems);
            printf("[UNIT-TEST]: -- END: cmean_arithm_u10x_zmm16r4_a\n");
            printf("mean_re=%.7f, mean_im=%.7f\n",mean_re,mean_im);
        }
        break;
        case 2 : 
        {
            std::lognormal_distribution<float> lognormal_distro{};
            char * __restrict__ distro_name{NULL};
            int32_t status{};
            distro_name = gms::common::demangle(typeid(lognormal_distro).name(),status);
            if(status==0 && distro_name != NULL)
            {
               printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", distro_name);
               gms::common::gms_mm_free(distro_name);
            }
            else
            {
              printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognormal_distro).name());
              if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            printf("UNIT-TEST]: -- Start of data initialization\n");
            seedre = std::clock();
            seedim = std::clock();
            auto srand_re{std::bind(std::lognormal_distribution<float>{1.6f,0.25f},
                                 std::mt19937(seedre))};
            auto srand_im{std::bind(std::lognormal_distribution<float>{1.6f,0.25f},
                                 std::mt19937(seedim))};
            std::generate((float*)&buf_re[0],(float*)&buf_re[0]+nelems,srand_re);
            std::generate((float*)&buf_im[0],(float*)&buf_im[0]+nelems,srand_im);
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cmean_arithm_u10x_zmm16r4_a\n");
            gms::math::cmean_arithm_u10x_zmm16r4_a(&buf_re[0],&buf_im[0],
                                                   &mean_re,&mean_im,nelems);
            printf("[UNIT-TEST]: -- END: cmean_arithm_u10x_zmm16r4_a\n");
            printf("mean_re=%.7f, mean_im=%.7f\n",mean_re,mean_im);
        }
        break;
        case 3 : 
        {
            std::fisher_f_distribution<float> fisher_f_distro{};
            char * __restrict__ distro_name{NULL};
            int32_t status{};
            distro_name = gms::common::demangle(typeid(fisher_f_distro).name(),status);
            if(status==0 && distro_name != NULL)
            {
               printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", distro_name);
               gms::common::gms_mm_free(distro_name);
            }
            else
            {
              printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(fisher_f_distro).name());
              if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            printf("UNIT-TEST]: -- Start of data initialization\n");
            seedre = std::clock();
            seedim = std::clock();
            auto srand_re{std::bind(std::fisher_f_distribution<float>{1.0f,5.0f},
                                 std::mt19937(seedre))};
            auto srand_im{std::bind(std::fisher_f_distribution<float>{1.0f,5.0f},
                                 std::mt19937(seedim))};
            std::generate((float*)&buf_re[0],(float*)&buf_re[0]+nelems,srand_re);
            std::generate((float*)&buf_im[0],(float*)&buf_im[0]+nelems,srand_im);
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cmean_arithm_u10x_zmm16r4_a\n");
            gms::math::cmean_arithm_u10x_zmm16r4_a(&buf_re[0],&buf_im[0],
                                                   &mean_re,&mean_im,nelems);
            printf("[UNIT-TEST]: -- END: cmean_arithm_u10x_zmm16r4_a\n");
            printf("mean_re=%.7f, mean_im=%.7f\n",mean_re,mean_im);
        }
        break;
        case 4 : 
        {
            std::cauchy_distribution<float> cauchy_distro{};
            char * __restrict__ distro_name{NULL};
            int32_t status{};
            distro_name = gms::common::demangle(typeid(cauchy_distro).name(),status);
            if(status==0 && distro_name != NULL)
            {
               printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", distro_name);
               gms::common::gms_mm_free(distro_name);
            }
            else
            {
              printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(cauchy_distro).name());
              if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            printf("UNIT-TEST]: -- Start of data initialization\n");
            seedre = std::clock();
            seedim = std::clock();
            auto srand_re{std::bind(std::cauchy_distribution<float>{-1.0f,5.0f},
                                 std::mt19937(seedre))};
            auto srand_im{std::bind(std::cauchy_distribution<float>{-1.0f,5.0f},
                                 std::mt19937(seedim))};
            std::generate((float*)&buf_re[0],(float*)&buf_re[0]+nelems,srand_re);
            std::generate((float*)&buf_im[0],(float*)&buf_im[0]+nelems,srand_im);
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cmean_arithm_u10x_zmm16r4_a\n");
            gms::math::cmean_arithm_u10x_zmm16r4_a(&buf_re[0],&buf_im[0],
                                                   &mean_re,&mean_im,nelems);
            printf("[UNIT-TEST]: -- END: cmean_arithm_u10x_zmm16r4_a\n");
            printf("mean_re=%.7f, mean_im=%.7f\n",mean_re,mean_im);
        }
        break;
        default : 
            printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
            std::terminate();
     }
}

int main()
{
    unit_test_carithm_mean_zmm16r4_u();
    unit_test_carithm_mean_zmm16r4_a();
    return 0;
}