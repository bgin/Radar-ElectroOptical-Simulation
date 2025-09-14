



#include "GMS_radar_jamming_kernels.cuh"
#include "GMS_radar_jamming_common.cuh"


static __device__
float therm_noise_range(   const float Frdr,
                           const float Kth,
                           const float rho,
                           const float tf,
                           const float tr,
	                   const float th,
		           const float d_Pt,
		           const float gamm,
			   const float d_w,
			   const float d_h,
			   const float Ln,
			   const float d_Ts,
			   const float sig,
			   const float F,
			   const float Fp,
		           const float Flens,
			   const float Dx,
			   const float d_Lt,
			   const float d_La) {

        const float dc  = duty_cycle(rho,tr);
        const float Pav = radar_avg_pow(d_Pt,dc);
        const float ag  = radar_ant_gain(azimuth_bw(Kth,gamm,
                                                    d_h),
                                         elevation_bw(Kth,gamm,
                                                    d_w),Ln);
        const float N0  = noise_density(d_Ts);
        const float den = 1984.4017075391884912304967f*N0*Dx*d_Lt*d_La;
        const float t1  = gamm*gamm;
        const float t2  = Pav*tf;
        const float t3  = ag*ag;
        const float t4  = sig*Frdr*Frdr*Fp*Fp;
        const float t5  = F*F*F*F*Flens*Flens;
        const float num = t1*t2*t3*t4*t5;
        const float rat = num/den;
        return (powf(rat,0.25f));
}


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


void therm_noise_range1_cuda( const float Frdr,
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

         uint threadsBlock = 32;
         uint blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
         therm_noise_range_kernel1<<<blocksGrid,threadsBlock>>>( Frdr,Kth,rho,tf,tr,th,d_Pt,gamm,
                                                                 d_w,d_h,Ln,d_Ts,sig,F,Fp,Flens,Dx,
                                                                 d_Lt,d_La,d_Rm,n_threads);
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

void therm_noise_range2_cuda( const float Frdr,
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

         uint threadsBlock = 256;
         uint blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
         therm_noise_range_kernel2<<<blocksGrid,threadsBlock>>>( Frdr,Kth,rho,tf,tr,th,d_Pt,gamm,
                                                                 d_w,d_h,Ln,d_Ts,sig,F,Fp,Flens,Dx,
                                                                 d_Lt,d_La,d_Rm,n);


static __device__
float tropo_range_loss(     const float Frdr,
                            const float Kth,
                            const float rho,
                            const float tf,
                            const float tr,
	                    const float th,
		            const float Pt,
			    const float Rmj,
			    const float gamm,
			    const float w,
			    const float h,
			    const float Ln,
			    const float Ts,
			    const float sig,
			    const float F,
			    const float Fp,
		            const float Flens,
			    const float Dx,
			    const float Lt,
			    const float La,
			    float Rm) {

           
           const float Rm = therm_noise_range(Frdr,Kth,rho,tf,tr,th,
                                              Pt,Rmj,gamm,
                                              w,h,Ln,Ts,
                                              sig,F,FP,Flens,Dx,Lt,La,Rm);
           return (La*(Rmj/Rm));                              
} 


__global__
void jammer_req_temp_kernel1( const float Frdr,
                              const float Kth,
                              const float rho,
                              const float tf,
                              const float tr,
                              const float th,
                              const float * __restrict__ d_Pt,
                              const float gamm,
                              const float * __restrict__ d_w,
                              const float * __restrict__ d_h,
                              const float Ln,
                              const float * __restrict__ d_Ts,
                              const float sig,
                              const float F,
                              const float Fp,
                              const float Flens,
                              const float Dx,
                              const float * __restrict__ d_Lt,
                              const float * __restrict__ d_Rm,
                              const float Rmj,
                              const float La,
                              const float Flen,
                              float * __restrict__ rt,
                              const uint32_t n_threads) {

    uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
    
    if(tid < n_threads) {
       Ts               = d_Ts[tid];
       const float Rm   = therm_noise_range(Frdr,Kth,rho,tf,tr,th,
                                            d_Pt[tid],Rmj,gamm,
                                            d_w[tid],d_h[tid],Ln,d_Ts[tid],
                                            sig,F,FP,Flens,Dx,d_Lt[tid],La,d_Rm[tid]);
                                            
       const float La1  = tropo_range_loss(Frdr,Kth,rho,tf,tr,th,d_Pt[tid],
                                           Rmj,gamm,d_w[tid],d_h[tid],Ln,d_Ts[tid],
                                           sig,F,Fp,Flens,Dx,d_Lt[tid],La);
       const float Fln1 = sqrtf(Flen);
       const float lrat = La/La1;
       const float mrat = Flen/Fln1;
       const float rrat = Rm/Rmj;
       const float rrat4= (rrat*rrat*rrat*rrat)-1.0f;
       rt[tid]          = Ts*lrat*(mrat*mrat)*rrat4;
    }
}

__global__
void jammer_req_temp_kernel2( const float Frdr,
                              const float Kth,
                              const float rho,
                              const float tf,
                              const float tr,
                              const float th,
                              const float * __restrict__ d_Pt,
                              const float gamm,
                              const float * __restrict__ d_w,
                              const float * __restrict__ d_h,
                              const float Ln,
                              const float * __restrict__ d_Ts,
                              const float sig,
                              const float F,
                              const float Fp,
                              const float Flens,
                              const float Dx,
                              const float * __restrict__ d_Lt,
                              const float * __restrict__ d_Rm,
                              const float Rmj,
                              const float La,
                              const float Flen,
                              float * __restrict__ rt,
                              const uint32_t n) {

   
    uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
    uint32_t stride = blockDim.x*gridDim.x;
    for(uint32_t i = tid; i < n; i += stride) {
       Ts               = d_Ts[tid];
       const float Rm   = therm_noise_range(Frdr,Kth,rho,tf,tr,th,
                                            d_Pt[i],Rmj,gamm,
                                            d_w[i],d_h[i],Ln,d_Ts[i],
                                            sig,F,FP,Flens,Dx,d_Lt[i],La,d_Rm[i]);
                                            
       const float La1  = tropo_range_loss(Frdr,Kth,rho,tf,tr,th,d_Pt[i],
                                           Rmj,gamm,d_w[i],d_h[i],Ln,d_Ts[i],
                                           sig,F,Fp,Flens,Dx,d_Lt[i],La);
       const float Fln1 = sqrtf(Flen);
       const float lrat = La/La1;
       const float mrat = Flen/Fln1;
       const float rrat = Rm/Rmj;
       const float rrat4= (rrat*rrat*rrat*rrat)-1.0f;
       rt[i]          = Ts*lrat*(mrat*mrat)*rrat4;
    }
}

	 
void jammer_req_temp_cuda(    const float Frdr,
                              const float Kth,
                              const float rho,
                              const float tf,
                              const float tr,
                              const float th,
                              const float * __restrict__ d_Pt,
                              const float gamm,
                              const float * __restrict__ d_w,
                              const float * __restrict__ d_h,
                              const float Ln,
                              const float * __restrict__ d_Ts,
                              const float sig,
                              const float F,
                              const float Fp,
                              const float Flens,
                              const float Dx,
                              const float * __restrict__ d_Lt,
                              const float * __restrict__ d_Rm,
                              const float Rmj,
                              const float La,
                              const float Flen,
                              float * __restrict__ rt,
                              const uint32_t n_threads,
                              const uint32_t type,
                              const uint32_t n) {

     if(type==1U) {
         uint threadsBlock = 32;
         uint blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
         jammer_req_temp_kernel1<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,d_Pt,
                                                              gamm,d_w,d_h,Ln,d_Ts,sig,F,
                                                              Fp,Flens,Dx,d_Lt,d_Rm,Rmj,
                                                              La,Flen,rt,n_threads);
     }
     else if(type==2U) {
         uint threadsBlock = 256;
         uint blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
         jammer_req_temp_kernel1<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,d_Pt,
                                                              gamm,d_w,d_h,Ln,d_Ts,sig,F,
                                                              Fp,Flens,Dx,d_Lt,d_Rm,Rmj,
                                                              La,Flen,rt,n); 
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

        
     /* if(type==1) {
         uint threadsBlock = 32;
         uint blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
	 therm_noise_range_kernel1<<<blocksGrid,threadsBlock>>>(Frdr,Kth,rho,tf,tr,th,
	                                                      d_Pt,gamm,d_w,d_h,Ln,
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
      }*/
	uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
        if(tid < n_threads) {
           const float La = d_La[tid];
	   const float Rmj= d_Rmj[tid];
           const float Rm = therm_noise_range(Frdr,Kth,rho,tf,tr,th,
                                              d_Pt[tid],Rmj,gamm,
                                              d_w[tid],d_h[tid],Ln,d_Ts[tid],
                                              sig,F,FP,Flens,Dx,d_Lt[tid],La,
                                              d_Rm[tid]);
           
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

        
      /*if(type==1) {
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
      }*/
      uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
      uint32_t stride = blockDim.x*gridDim.x;
      for(uint32_t i = tid; i < n; i += stride) {
           const float La = d_La[i];
	   const float Rmj= d_Rmj[i];
           const float Rm = therm_noise_range(Frdr,Kth,rho,tf,tr,th,
                                              d_Pt[i],Rmj,gamm,
                                              d_w[tid],d_h[i],Ln,d_Ts[i],
                                              sig,F,FP,Flens,Dx,d_Lt[i],La,
                                              d_Rm[i]);
           
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

 // Effective radiated power of Jammer (W)
/*
                            float sig;  //m, RSC of target
			    float Pj;   //W, jammer power
			    float Gj;   //dB, jammer antenna gain
			    float Qj;   //dB, jammer noise quality
			    float Flenj;//dB, jammer lens factor
			    float Rj;   //km, jammer range
			    float Bj;   //Mhz,jammer noise BW
			    float Ltj;  //dB, jammer transmit loss
			    float Fpj;  //dB, jammer polarization
			    float Rmj;  //km, jammer screening range
			    float Fj;   //dB, jammer pattern factor of propagation
			    float Laj;  //dB, jammer troposhperic loss
*/

__global__
void jammer_erp_kernel1(const float * __restrict__ Pj,
                        const float Gj,
                        const float Ltj,
                        float * __restrict__ erp,
                        const uint32_t n_threads) {

      uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
      if(tid < n_threads) {
         const float xPj = Pj[tid];
         erp[tid] = (xPj*Gj)/Ltj;
      }
}


__global__
void jammer_erp_kernel2(const float * __restrict__ Pj,
                        const float Gj,
                        const float Ltj,
                        float * __restrict__ erp,
                        const uint32_t n) {

     uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
     uint32_t stride = blockDim.x*gridDim.x;
     for(uint32_t i = tid; i < n; i += stride) {
          const float xPj = Pj[i];
          erp[i] = (xPj*Gj)/Ltj;
     }
}


void jammer_erp_cuda(const float * __restrict__ Pj,
                     const float Gj,
                     const float Ltj,
                     float * __restrict__ erp,
                     const uint32_t n_threads,
                     const uint32_t type,
                     const uint32_t n) {
    
     if(type==1) {
         uint32_t threadsBlock = 32;
         uint32_t blocksGrid   = (n_threads + threadsBlock - 1) / threadsBlock;
         jammer_erp_kernel1<<<blocksGrid,threadsBlock>>>(Pj,Gj,Ltj,erp,n_threads);
     }
     else if(type==2) {
         uint32_t threadsBlock = 256;
         uint32_t blocksGrid   = (n + threadsBlock - 1) / threadsBlock;
         jammer_erp_kernel2<<<blocksGrid,threadsBlock>>>(Pj,Gj,Ltj,erp,n);
     }
}


// Effective radiated Jammer noise power (W)
                        
__global__
void jammer_ernp_kernel1(const float * __restrict__ Pj,
                         const float Qj,
                         const float Gj,
                         const float Fpj,
                         const float Ltj,
                         float * __restrict__ ernp,
                         const uint32_t n_threads) {

      uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
      const float xQj = Qj;
      const float xGj = Gj;
      const float xFpj= Fpj;
      const float xLtj= Ltj;
      if(tid < n_threads) {
          const float xPj = Pj[tid];
          const float xFpj2 = xFpj*xFpj;
          ernp[tid] = (xQj*xPj*xGj*xFpj2)/xLtj;
      }
}


__global__
void jammer_ernp_kernel2(const float * __restrict__ Pj,
                         const float Qj,
                         const float Gj,
                         const float Fpj,
                         const float Ltj,
                         float * __restrict__ ernp,
                         const uint32_t n) {

     uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
     uint32_t stride = blockDim.x*gridDim.x; 
     const float xQj = Qj;
     const float xGj = Gj;
     const float xFpj= Fpj;
     const float xLtj= Ltj;
     for(uint32_t i = tid; i < n; i += stride) {
          const float xPj = Pj[i];
          const float xFpj2 = xFpj*xFpj;
          ernp[i] = (xQj*xPj*xGj*xFpj2)/xLtj;
     }
} 


void jammer_ernp_cuda(const float * __restrict__ Pj,
                      const float Qj,
                      const float Gj,
                      const float Fpj,
                      const float Ltj,
                      float * __restrict__ ernp,
                      const uint32_t n_threads,
                      const uint32_t type,
                      const uint32_t n) {

     if(type==1) {
         uint32_t threadsBlock = 32;
         uint32_t blocksGrid   = (n_threads + threadsBlock - 1) / threadsBlock;
         jammer_ernp_kernel1<<<blocksGrid,threadsBlock>>>(Pj,Qj,Gj,Fpj,Ltj,ernp,n_threads);
     }
     else if(type==2) {
         uint32_t threadsBlock = 256;
         uint32_t blocksGrid   = (n + threadsBlock - 1) / threadsBlock;
         jammer_ernp_kerne2<<<blocksGrid,threadsBlock>>>(Pj,Qj,Gj,Fpj,Ltj,ernp,n);
     }
}

 // Jamming spectral density (W/Hz)

__global__
void jammer_spectr_dens_kernel1(const float * __restrict__ gamma,
                                const float Kth,
                                const float Ln,
                                const float h,
                                const float w,
                                const float Qj,
                                const float * __restrict__ Pj,
                                const float Gj,
                                const float Fpj,
                                const float Flnsj,
                                const float Fj,
                                const float Rj,
                                const float Bj,
                                const float Ltj,
                                const float Laj,
                                float * __restrict__ sd,
                                const uint32_t n_threads) {
          
        const float j0   = 0.0f;
        const float PI42 =  157.9136704174297379013522f;  
        uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
        if(tid < n_threads) {
           const float xgamma = gamma[tid];
           const float xPj    = Pj[tid];
           const float xgm2   = xgamma*xgamma;
           const float xGr    = radar_ant_gain(azimuth_bw(Kth,xgamma,h),
                                               elevation_bw(Kth,xgamma,w),Ln);
           const float r0     = Rj*Bj*Ltj*Laj;
           const float t0     = Qj*xPj*Gj*Gr;
           const float t1     = xgm2*Fpj*Fpj*Flnsj*Flnsj;
           const float t2     = t0*t1*Fj*Fj;
           j0                 = t2/(PI42*r0);
           sd[tid]            = j0;
        } 
}


__global__
void jammer_spectr_dens_kernel2(const float * __restrict__ gamma,
                                const float Kth,
                                const float Ln,
                                const float h,
                                const float w,
                                const float Qj,
                                const float * __restrict__ Pj,
                                const float Gj,
                                const float Fpj,
                                const float Flnsj,
                                const float Fj,
                                const float Rj,
                                const float Bj,
                                const float Ltj,
                                const float Laj,
                                float * __restrict__ sd,
                                const uint32_t n) {

      const float j0   = 0.0f;
      const float PI42 =  157.9136704174297379013522f;
      uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
      uint32_t stride = blockDim.x*gridDim.x; 
      for(uint32_t i = tid; i < n; i += stride) {
           const float xgamma = gamma[i];
           const float xPj    = Pj[i];
           const float xgm2   = xgamma*xgamma;
           const float xGr    = radar_ant_gain(azimuth_bw(Kth,xgamma,h),
                                               elevation_bw(Kth,xgamma,w),Ln);
           const float r0     = Rj*Bj*Ltj*Laj;
           const float t0     = Qj*xPj*Gj*Gr;
           const float t1     = xgm2*Fpj*Fpj*Flnsj*Flnsj;
           const float t2     = t0*t1*Fj*Fj;
           j0                 = t2/(PI42*r0);
           sd[i]            = j0;
      }  
}


void jammer_spectr_dens_cuda(   const float * __restrict__ gamma,
                                const float Kth,
                                const float Ln,
                                const float h,
                                const float w,
                                const float Qj,
                                const float * __restrict__ Pj,
                                const float Gj,
                                const float Fpj,
                                const float Flnsj,
                                const float Fj,
                                const float Rj,
                                const float Bj,
                                const float Ltj,
                                const float Laj,
                                float * __restrict__ sd,
                                const uint32_t n_threads,
                                const uint32_t type,
                                const uint32_t n) {

     if(type==1) {
         uint32_t threadsBlock = 32;
         uint32_t blocksGrid   = (n_threads + threadsBlock - 1) / threadsBlock;
         jammer_spectr_dens_kernel1<<<blocksGrid,threadsBlock>>>(gamma,Kth,Ln,h,w,Qj,Pj,
                                                                 Gj,Fpj,Flnsj,Fj,Rj,Bj,
                                                                 Ltj,Laj,sd,n_threads);
     }
     else if(type==2) {
         uint32_t threadsBlock = 256;
         uint32_t blocksGrid   = (n + threadsBlock - 1) / threadsBlock;
         jammer_spectr_dens_kernel2<<<blocksGrid,threadsBlock>>>(gamma,Kth,Ln,h,w,Qj,Pj,
                                                                 Gj,Fpj,Flnsj,Fj,Rj,Bj,
                                                                 Ltj,Laj,sd,n);
     }
}

// Available jamming temperature single jammer scenario
// Remark: the implementation is identical to function coded above.

__global__
void single_jammer_temp_kernel1(const float * __restrict__ gamma,
                                const float Kth,
                                const float Ln,
                                const float h,
                                const float w,
                                const float Qj,
                                const float * __restrict__ Pj,
                                const float Gj,
                                const float Fpj,
                                const float Flnsj,
                                const float Fj,
                                const float Rj,
                                const float Bj,
                                const float Ltj,
                                const float Laj,
                                float * __restrict__ sd,
                                const uint32_t n_threads) {
          
        const float jt1   = 0.0f;
        const float PI42 =  157.9136704174297379013522f;  
        uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
        if(tid < n_threads) {
           const float xgamma = gamma[tid];
           const float xPj    = Pj[tid];
           const float xgm2   = xgamma*xgamma;
           const float xGr    = radar_ant_gain(azimuth_bw(Kth,xgamma,h),
                                               elevation_bw(Kth,xgamma,w),Ln);
           const float r0     = Rj*k_B*Bj*Ltj*Laj;
           const float t0     = Qj*xPj*Gj*Gr;
           const float t1     = xgm2*Fpj*Fpj*Flnsj*Flnsj;
           const float t2     = t0*t1*Fj*Fj;
           jt1                = t2/(PI42*r0);
           sd[tid]            = jt1;
        } 
}


__global__
void single_jammer_temp_kernel2(const float * __restrict__ gamma,
                                const float Kth,
                                const float Ln,
                                const float h,
                                const float w,
                                const float Qj,
                                const float * __restrict__ Pj,
                                const float Gj,
                                const float Fpj,
                                const float Flnsj,
                                const float Fj,
                                const float Rj,
                                const float Bj,
                                const float Ltj,
                                const float Laj,
                                float * __restrict__ sd,
                                const uint32_t n) {

      const float jt1  = 0.0f;
      const float PI42 =  157.9136704174297379013522f;
      uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
      uint32_t stride = blockDim.x*gridDim.x; 
      for(uint32_t i = tid; i < n; i += stride) {
           const float xgamma = gamma[i];
           const float xPj    = Pj[i];
           const float xgm2   = xgamma*xgamma;
           const float xGr    = radar_ant_gain(azimuth_bw(Kth,xgamma,h),
                                               elevation_bw(Kth,xgamma,w),Ln);
           const float r0     = Rj*k_B*Bj*Ltj*Laj;
           const float t0     = Qj*xPj*Gj*Gr;
           const float t1     = xgm2*Fpj*Fpj*Flnsj*Flnsj;
           const float t2     = t0*t1*Fj*Fj;
           jt1                = t2/(PI42*r0);
           sd[i]              = jt1;
      }  
}


void single_jammer_temp_cuda(   const float * __restrict__ gamma,
                                const float Kth,
                                const float Ln,
                                const float h,
                                const float w,
                                const float Qj,
                                const float * __restrict__ Pj,
                                const float Gj,
                                const float Fpj,
                                const float Flnsj,
                                const float Fj,
                                const float Rj,
                                const float Bj,
                                const float Ltj,
                                const float Laj,
                                float * __restrict__ sd,
                                const uint32_t n_threads,
                                const uint32_t type,
                                const uint32_t n) {

     if(type==1) {
         uint32_t threadsBlock = 32;
         uint32_t blocksGrid   = (n_threads + threadsBlock - 1) / threadsBlock;
         single_jammer_temp_kernel1<<<blocksGrid,threadsBlock>>>(gamma,Kth,Ln,h,w,Qj,Pj,
                                                                 Gj,Fpj,Flnsj,Fj,Rj,Bj,
                                                                 Ltj,Laj,sd,n_threads);
     }
     else if(type==2) {
         uint32_t threadsBlock = 256;
         uint32_t blocksGrid   = (n + threadsBlock - 1) / threadsBlock;
         single_jammer_temp_kernel2<<<blocksGrid,threadsBlock>>>(gamma,Kth,Ln,h,w,Qj,Pj,
                                                                 Gj,Fpj,Flnsj,Fj,Rj,Bj,
                                                                 Ltj,Laj,sd,n);
     }
}


//// Number of jammers range (km)
__global__
void n_jammers_range_kernel1(const float Kth,
                             const float gamma,
                             const float * __restrict__ h,
                             const float * __restrict__ w,
                             const float Ln,
                             const float rho,
                             const float tf,
                             const float tr,
                             const float * __restrict__ Pt,
                             const float Frdr,
                             const float Fp,
                             const float F,
                             const float Flen,
                             const float * __restrict__ sig,
                             const float Ts,
                             const float Dx,
                             const float Lt,
                             const float La,
                             float * __restrict Rnj,
                             const uint32_t n_threads) {
       
       const float c0 = 1984.4017075391884912304842f
       uint32_t tid   = blockDim.x*blockIdx.x+threadIdx.x;
       if(tid < n_threads) {
          const float xh     = h[tid];
          const float xw     = w[tid];
          const float xPt    = Pt[tid];
          const float xsig   = sig[tid];
          const float gamma2 = gamma*gamma;
          const float Gr     = radar_ant_gain(azimuth_bw(Kth,gamma2,xh),
                                              elevation_bw(Kth,gamma2,xw),Ln);
          const float Gr2    = Gr*Gr;
          const float dc     = duty_cycle(rho,tr);
          const float pav    = radar_avg_power(xPt,dc);
          const float N0     = noise_density(Ts);
          const float F4     = F*F*F*F;
          const float num    = pav*tf*Gr2*gamma*xsig*F4*Fp;
          const float den    = c0*k_B*Ts*Dx*Lt*La;
          const float ratio  = num/den;
          Rnj[tid]           = powf(ratio,0.25f);
       } 
}

__global__
void n_jammers_range_kernel2(const float Kth,
                             const float gamma,
                             const float * __restrict__ h,
                             const float * __restrict__ w,
                             const float Ln,
                             const float rho,
                             const float tf,
                             const float tr,
                             const float * __restrict__ Pt,
                             const float Frdr,
                             const float Fp,
                             const float F,
                             const float Flen,
                             const float * __restrict__ sig,
                             const float Ts,
                             const float Dx,
                             const float Lt,
                             const float La,
                             float * __restrict Rnj,
                             const uint32_t n) {
       
       const float c0 = 1984.4017075391884912304842f
       uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
       uint32_t stride = blockDim.x*gridDim.x; 
       for(uint32_t i = tid; i < n; i += stride) {   
           const float xh     = h[i];
           const float xw     = w[i];
           const float xPt    = Pt[i];
           const float xsig   = sig[i];
           const float gamma2 = gamma*gamma;
           const float Gr     = radar_ant_gain(azimuth_bw(Kth,gamma2,xh),
                                              elevation_bw(Kth,gamma2,xw),Ln);
           const float Gr2    = Gr*Gr;
           const float dc     = duty_cycle(rho,tr);
           const float pav    = radar_avg_power(xPt,dc);
           const float N0     = noise_density(Ts);
           const float F4     = F*F*F*F;
           const float num    = pav*tf*Gr2*gamma*xsig*F4*Fp;
           const float den    = c0*k_B*Ts*Dx*Lt*La;
           const float ratio  = num/den;
           Rnj[i]           = powf(ratio,0.25f);
       } 
} 


void n_jammers_range_cuda(   const float Kth,
                             const float gamma,
                             const float * __restrict__ h,
                             const float * __restrict__ w,
                             const float Ln,
                             const float rho,
                             const float tf,
                             const float tr,
                             const float * __restrict__ Pt,
                             const float Frdr,
                             const float Fp,
                             const float F,
                             const float Flen,
                             const float * __restrict__ sig,
                             const float Ts,
                             const float Dx,
                             const float Lt,
                             const float La,
                             float * __restrict Rnj,
                             const uint32_t n_threads,
                             const uint32_t type,
                             const uint32_t n) {
     
     if(type==1U) {
         uint32_t threadsBlock = 32;
         uint32_t blocksGrid   = (n_threads + threadsBlock - 1) / threadsBlock;
         n_jammers_range_kernel1<<<blocksGrid,threadsBlock>>>(Kth,gamma,h,w,Ln,rho,tf,tr,
                                                              Pt,Frdr,Fp,F,Flen,sig,Ts,Dx,
                                                              Lt,La,Rnj,n_threads);
     }
     else if(type==2U) {
          uint32_t threadsBlock = 256;
          uint32_t blocksGrid   = (n + threadsBlock - 1) / threadsBlock;
          n_jammers_range_kernel2<<<blocksGrid,threadsBlock>>>(Kth,gamma,h,w,Ln,rho,tf,tr,
                                                              Pt,Frdr,Fp,F,Flen,sig,Ts,Dx,
                                                              Lt,La,Rnj,n);
     }
} 



 
                                  
