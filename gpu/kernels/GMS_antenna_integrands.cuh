
#ifndef __GMS_ANTENNA_INTEGRANDS_CUH__
#define __GMS_ANTENNA_INTEGRANDS_CUH__



#include <cuda/std/complex>

using namespace cuda::std;

// Integrand: (2-42)
    
      struct Intf242 {
      
          __device__ complex<float>
          operator()(const complex<float> ikz,
                     const float th,
                     const float z,
                     const float k,
                     const float m,
                     const float I0) {
             
             constexpr float pi = 3.14159265358979323846264f;
             complex<float> earg;
             complex<float> cexp;
             complex<float> value;
             float cost,mpi2,kz,sarg;
             mpi2 = m*pi*0.5f;
             cost = cosf(th);
             kz   = k*z;
             sarg = I0*sinf(kz+mpi2);
             earg = ikz*cost;
             cexp = exp(earg);
             value= sarg*cexp;
             return (value);     
         }
    };
    
    
// Integrand: (2-69)
   
   

      
   

   
   













#endif /*__GMS_ANTENNA_INTEGRANDS_CUH__*/
