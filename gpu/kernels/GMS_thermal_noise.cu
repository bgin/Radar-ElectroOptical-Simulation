

#include <cuda_runtime.h>
#include <cuda.h>
#include "GMS_cpu_config.cuh"
#include "GMS_thermal_noise.cuh"
#include "GMS_jamming_common.cuh"
#if (SOFTWARE_PREFETCH) == 1
#include "GMS_cuda_soft_prefetch.cuh"
#endif

__global__ void
therm_noise_range_kernel(  const float Frdr,
                           const float Kth,
                           const float rho,
                           const float tf,
                           const float tr,
	                   const float th,
		           const float * __restrict d_Pt,
		           const float gamm,
			   const float * __restrict d_w,
			   const float * __restrict d_h,
			   const float Ln,
			   const float * __restrict d_Ts,
			   const float sig,
			   const float F,
			   const float Fp,
		           const float Flens,
			   const float Dx,
			   const float * __restrict d_Lt,
			   const float * __restrict d_La,
			   float * __restrict d_Rm,
			   const int nth) {

     int tid = blockDim.x*blockIdx.x+threadIdx.x;
     if(tid < nth) {
#if (SOFTWARE_PREFETCH) == 1
      __prefetch_global_l1(d_Pt);
      __prefetch_global_l1(d_h);
      __prefetch_global_l1(d_w);
      __prefetch_global_l1(d_Ts);
      __prefetch_global_l1(d_Lt);
      __prefetch_global_l1(d_La);
#endif
        const float dc  = duty_cycle(rho,tr);
        const float Pav = radar_avg_pow(d_Pt[tid],dc);
        const float ag  = radar_ant_gain(azimuth_bw(Kth,gamm,
                                                    d_h[tid]),
                                         elevation_bw(Kth,gamm,
                                                    d_w[tid]),Ln);
        const float N0  = noise_density(d_Ts[tid]);
        const float den = 1984.4017075391884912304967f*N0*Dx*d_Lt[tid]*d_La[tid];
        const float t1  = gamm*gamm;
        const float t2  = Pav*tf;
        const float t3  = ag*ag;
        const float t4  = sig*Frdr*Frdr*Fp*Fp;
        const float t5  = F*F*F*F*Flens*Flens;
        const float num = t1*t2*t3*t4*t5;
        const float rat = num/den;
        Rm[tid]         = powf(rat,0.25f);
     }
}


__global__ void
tropo_range_loss_kernel(   const float Frdr,
                           const float Kth,
                           const float rho,
                           const float tf,
                           const float tr,
	                   const float th,
		           const float * __restrict d_Pt,
			   const float * __restrict d_Rmj,
			   const float gamm,
			   const float * __restrict d_w,
			   const float * __restrict d_h,
			   const float Ln,
			   const float * __restrict d_Ts,
			   const float sig,
			   const float F,
			   const float Fp,
		           const float Flens,
			   const float Dx,
			   const float * __restrict d_Lt,
			   const float * __restrict d_La,
			   float * __restrict d_Rm,
			   float * __restrict d_La1,
			   const int nth) {

        
        int threadsBlock = 256;
        int blocksGrid  = (nth + threadsBlock - 1) / threadsBlock;
	therm_noise_range_kernel<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,
	                                                      d_Pt,gamm.d_w,d_h,Ln,
							      d_Ts,sig,F,Fp,Flens,
							      Dx,d_Lt,d_La,d_Rm,nth);
	int tid = blockDim.x*blockIdx.x+threadIdx.x;
        if(tid < nth) {
#if (SOFTWARE_PREFETCH) == 1
	   __prefetch_global_l1(d_La);
	   __prefetch_global_l1(d_Rmj);
	   __prefetch_global_l1(d_Rm);
#endif
           const float La = d_La[tid];
	   const float Rmj= d_Rmj[tid];
	   const float Rm = d_Rm[tid];
	   d_La1[tid]     = La*(Rmj/Rm);
        }
}


 void
 therm_noise_range_cuda(           const float Frdr,
                                   const float Kth,
                                   const float rho,
                                   const float tf,
                                   const float tr,
				   const float th,
				   const float * __restrict d_Pt,
				   const float gamm,
				   const float * __restrict d_w,
				   const float * __restrict d_h,
				   const float Ln,
				   const float * __restrict d_Ts,
				   const float sig,
				   const float F,
				   const float Fp,
				   const float Flens,
				   const float Dx,
				   const float * __restrict d_Lt,
				   const float * __restrict d_La,
				   float * __restrict d_Rm,
				   const int nth) {
				  
          int threadsBlock = 32;
          int blocksGrid  = (nth + threadsBlock - 1) / threadsBlock;
          therm_noise_range_kernel<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,
	                                                        d_Pt,gamm.d_w,d_h,Ln,
							        d_Ts,sig,F,Fp,Flens,
							        Dx,d_Lt,d_La,d_Rm,nth);


}

 void
 tropo_range_loss_cuda(            const float Frdr,
                                   const float Kth,
                                   const float rho,
                                   const float tf
                                   const float tr,
				   const float th,
				   const float * __restrict d_Pt,
				   const float * __restrict d_Rmj,
				   const float gamm,
				   const float * __restrict d_w,
				   const float * __restrict d_h,
				   const float Ln,
				   const float * __restrict d_Ts,
				   const float sig,
				   const float F,
				   const float Fp,
				   const float Flens,
				   const float Dx,
				   const float * __restrict d_Lt,
				   const float * __restrict d_La,
				   float * __restrict d_Rm,
				   float * __restrict d_La1,
				   const int nth) {

          int threadsBlock = 32;
          int blocksGrid  = (nth + threadsBlock - 1) / threadsBlock;
	  tropo_range_loss_kernel<<<blocksGrid,threadsBlock>>>( Frdr,Kth,rho,tf,tr,d_Pt,d_Rmj,
	                                                        gamm,d_w,d_h,Ln,d_Ts,sig,F,Fp,
								Flens,Dx,d_Lt,d_La,d_Rm,d_La1);
}
