



#include "GMS_jamming_kernels.cuh"
#include "GMS_jamming_common.cuh"





__global__ void
therm_noise_range_kernel1(  const float Frdr,
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
			   const uint32_t n_threads) {

     uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
     if(tid < n_threads) {
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
therm_noise_range_kernel2(  const float Frdr,
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
			   const uint32_t n) {

     uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
     uint32_t stride = blockDim.x*gridDim.x;
     for(uint32_t i = tid; i < n; i += stride) {
         const float dc  = duty_cycle(rho,tr);
         const float Pav = radar_avg_pow(d_Pt[i],dc);
         const float ag  = radar_ant_gain(azimuth_bw(Kth,gamm,
                                                    d_h[i]),
                                         elevation_bw(Kth,gamm,
                                                    d_w[i]),Ln);
         const float N0  = noise_density(d_Ts[i]);
         const float den = 1984.4017075391884912304967f*N0*Dx*d_Lt[i]*d_La[i];
         const float t1  = gamm*gamm;
         const float t2  = Pav*tf;
         const float t3  = ag*ag;
         const float t4  = sig*Frdr*Frdr*Fp*Fp;
         const float t5  = F*F*F*F*Flens*Flens;
         const float num = t1*t2*t3*t4*t5;
         const float rat = num/den;
         Rm[i]         = powf(rat,0.25f);
     }
}



__global__ void
tropo_range_loss_kernel1(   const float Frdr,
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
			    const uint n_threads,
                            const int32_t type,
                            const uint32_t n) {

        
      if(type==1) {
         uint threadsBlock = 32;
         uint blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
	 therm_noise_range_kernel1<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,
	                                                      d_Pt,gamm.d_w,d_h,Ln,
							      d_Ts,sig,F,Fp,Flens,
							      Dx,d_Lt,d_La,d_Rm,n_threads);
      }
      else if(type==2) {
         uint threadsBlock = 256;
         uint blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
	 therm_noise_range_kernel2<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,
	                                                      d_Pt,gamm.d_w,d_h,Ln,
							      d_Ts,sig,F,Fp,Flens,
							      Dx,d_Lt,d_La,d_Rm,n);
      }
	uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
        if(tid < n_threads) {
           const float La = d_La[tid];
	   const float Rmj= d_Rmj[tid];
	   const float Rm = d_Rm[tid];
	   d_La1[tid]     = La*(Rmj/Rm);
        }
}


__global__ void
tropo_range_loss_kernel2(   const float Frdr,
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
			    const uint n_threads,
                            const int32_t type,
                            const uint32_t n) {

        
      if(type==1) {
         uint threadsBlock = 32;
         uint blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
	 therm_noise_range_kernel1<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,
	                                                      d_Pt,gamm.d_w,d_h,Ln,
							      d_Ts,sig,F,Fp,Flens,
							      Dx,d_Lt,d_La,d_Rm,n_threads);
      }
      else if(type==2) {
         uint threadsBlock = 256;
         uint blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
	 therm_noise_range_kernel2<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,
	                                                      d_Pt,gamm.d_w,d_h,Ln,
							      d_Ts,sig,F,Fp,Flens,
							      Dx,d_Lt,d_La,d_Rm,n);
      }
      uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
      uint32_t stride = blockDim.x*gridDim.x;
      for(uint32_t i = tid; i < n; i += stride) {
           const float La = d_La[i];
	   const float Rmj= d_Rmj[i];
	   const float Rm = d_Rm[i];
	   d_La1[i]     = La*(Rmj/Rm);
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
				   const uint32_t n_threads,
                                   const uint32_t type,
                                   const uint32_t n) {
	  
          if(type==1) {			  
             uint32_t threadsBlock = 32;
             uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
             therm_noise_range_kernel1<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,
	                                                            d_Pt,gamm.d_w,d_h,Ln,
							            d_Ts,sig,F,Fp,Flens,
							            Dx,d_Lt,d_La,d_Rm,n_threads);
           }
           else if(type==2) {
             uint32_t threadsBlock = 256;
             uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
             therm_noise_range_kernel2<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,
	                                                            d_Pt,gamm.d_w,d_h,Ln,
							            d_Ts,sig,F,Fp,Flens,
							            Dx,d_Lt,d_La,d_Rm,n);
           }


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
				   const uint32_t n_threads,
                                   const uint32_t type,
                                   const uint32_t n) {
         
         if(type==1) {
            uint32_t threadsBlock = 32;
            uint32_t blocksGrid   = (n_threads + threadsBlock - 1) / threadsBlock;
	    tropo_range_loss_kernel1<<<blocksGrid,threadsBlock>>>( Frdr,Kth,rho,tf,tr,d_Pt,d_Rmj,
	                                                           gamm,d_w,d_h,Ln,d_Ts,sig,F,Fp,
								   Flens,Dx,d_Lt,d_La,d_Rm,d_La1,
                                                                   n_threads,type,n);
         }
         else if(type==2) {
            uint32_t threadsBlock = 256;
            uint32_t blocksGrid   = (n + threadsBlock - 1) / threadsBlock;
	    tropo_range_loss_kernel2<<<blocksGrid,threadsBlock>>>( Frdr,Kth,rho,tf,tr,d_Pt,d_Rmj,
	                                                           gamm,d_w,d_h,Ln,d_Ts,sig,F,Fp,
								   Flens,Dx,d_Lt,d_La,d_Rm,d_La1,
                                                                   n_threads,type,n);
         }
}
