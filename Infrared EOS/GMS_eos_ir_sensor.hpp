

#ifndef __GMS_EOS_IR_SENSOR_HPP__
#define __GMS_EOS_IR_SENSOR_HPP__ 050920220913

namespace file_version {

    const unsigned int GMS_EOS_IR_SENSOR_MAJOR = 1U;
    const unsigned int GMS_EOS_IR_SENSOR_MINOR = 0U;
    const unsigned int GMS_EOS_IR_SENSOR_MICRO = 0U;
    const unsigned int GMS_EOS_IR_SENSOR_FULLVER =
      1000U*GMS_EOS_IR_SENSOR_MAJOR+
      100U*GMS_EOS_IR_SENSOR_MINOR+
      10U*GMS_EOS_IR_SENSOR_MICRO;
    const char * const GMS_EOS_IR_SENSOR_CREATION_DATE = "05-09-2022 09:13 AM +00200 (MON 05 SEP 2022 GMT+2)";
    const char * const GMS_EOS_IR_SENSOR_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_EOS_IR_SENSOR_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_EOS_IR_SENSOR_DESCRIPTION   = "Based on book: Mathematical Theory of Electro-Optical Sensors (rus)."

}


/*
   Various characteristics of Electro-Optical Sensors   
 ! Based mainly on Based mainly on Miroshenko M.M book (rus):          
 ! "Mathematical Theory of Electro-Optical Sensors".
*/

#include <cstdint>
#include <omp.h>
#include <cmath>
#include "GMS_config.h"

namespace gms {

         namespace eos {


	      __ATTR_ALWAYS_INLINE__
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	      static inline
#pragma omp declare simd simdlen(16)
	      float param_gamma(const float phi) {

	            float gamma = 0.0f;
		    gamma       = 0.5f*phi*0.5f;
		    return (gamma);
	   }


              __ATTR_ALWAYS_INLINE__
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	      static inline
#pragma omp declare simd simdlen(8)
	      double param_gamma(const double phi) {

	            double gamma = 0.0;
		    gamma        = 0.5*phi*0.5;
		    return (gamma);
	   }


	   //! Formula 1, p.54
           //!Тогда длина перпендикуляра SN, опущенного из 
           //!светящейся точки на плоскость зеркала
	      __ATTR_ALWAYS_INLINE__
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	      static inline
#pragma omp declare simd simdlen(16)
	      float compute_SN(const float R,
	                       const float phi,
			       const float gamma,
			       const bool present) {

		     float SN = 0.0f;
		     if(present) {  //if true compute phi, otherwise compute gamma
                        SN = R*std::sin(phi);
		     }
		     else {
                        SN = R*std::cos(gamma);
		     }
		     return (SN);
	    }
	    

	      __ATTR_ALWAYS_INLINE__
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	      static inline
#pragma omp declare simd simdlen(8)
	      double compute_SN(const double R,
	                        const double phi,
			        const double gamma,
			        const bool   present) {

		     double SN = 0.0;
		     if(present) {  //if true compute phi, otherwise compute gamma
                        SN = R*std::sin(phi);
		     }
		     else {
                        SN = R*std::cos(gamma);
		     }
		     return (SN);
	    }


	    //! Formula 2, p. 54
            //! расстояние SM от светящейся точки до ее изображения
	      __ATTR_ALWAYS_INLINE__
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	      static inline
#pragma omp declare simd simdlen(16)
              float compute_SM(const float R,
	                       const float phi,
			       const float gamma,
			       const bool present) {

                  float SM = 0.0f;
		  SM       = 2.0f*compute_SN(R,phi,gamma,present);
		  return (SM);
	    }


	      __ATTR_ALWAYS_INLINE__
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	      static inline
#pragma omp declare simd simdlen(8)
              double compute_SM(const double R,
	                        const double phi,
			        const double gamma,
			        const bool present) {

                  double SM = 0.0;
		  SM        = 2.0*compute_SN(R,phi,gamma,present);
		  return (SM);
	    }


	   //!Сканирующее зеркало для обеспечения осмотра всего поля
           //!обзора ф необходимо повернуть на угол, обеспечивающий 
           //!совмещение края изображения источника излучения с отверстием 
           //!диафрагмы а, находящимся в центре поля. Для этого необходимо 
           //!повернуть изображение светящейся точки S на угол ф/2
           //! Formula 1, p. 56
	      __ATTR_ALWAYS_INLINE__
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	      static inline
#pragma omp declare simd simdlen(16)
              float ratio_FH(const float psi,
	                     const float phi) {

                    float FH = 0.0f;
		    float hpsi,hphi;
		    hpsi     = std::tan(0.5f*psi)
		    hphi     = std::tan(0.5f*phi);
		    FH       = hpsi/hphi;
		    return (FH);
	    }


	      __ATTR_ALWAYS_INLINE__
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	      static inline
#pragma omp declare simd simdlen(8)
              double ratio_FH(const double psi,
	                      const double phi) {
			      
                    double FH = 0.0;
		    double hpsi,hphi;
		    hpsi     = std::tan(0.5*psi)
		    hphi     = std::tan(0.5*phi);
		    FH       = hpsi/hphi;
		    return (FH);
	    }


	    //! следовательно, угол установки сканирующего зеркала
            //! Formula 4, p. 56
	     __ATTR_ALWAYS_INLINE__
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     static inline
#pragma omp declare simd simdlen(16)
	     float scan_mirror_ang(const float gam0,
	                           const float psi,
				   const float phi,
				   const int32_t dir) { // value '1' -- positive, value '2' negative direction

                float gamma,t0,t1;
		t1     = 0.5f*phi*0.5f;
		gamma  = 0.0f;
		if(dir==1) {
                   t0    = gam0+t1;
		   gamma = t0*ratio_FH(psi,phi); 
		}
		else if(dir==2) {
                   t0    = gam0-t1;
		   gamma = t0*ratio_FH(psi,phi);
		}
		return (gamma);
	    }


	     __ATTR_ALWAYS_INLINE__
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     static inline
#pragma omp declare simd simdlen(8)
	     double scan_mirror_ang(const double gam0,
	                            const double psi,
				    const double phi,
				    const int32_t dir) { // value '1' -- positive, value '2' negative direction

                double gamma,t0,t1;
		t1     = 0.5*phi*0.5;
		gamma  = 0.0;
		if(dir==1) {
                   t0    = gam0+t1;
		   gamma = t0*ratio_FH(psi,phi); 
		}
		else if(dir==2) {
                   t0    = gam0-t1;
		   gamma = t0*ratio_FH(psi,phi);
		}
		return (gamma);
	    }


	    //! Maximum size of (verical) diameter of scanning mirror.
            //! Formula 2, page. 56, part: 1.3
            //! Anax = [h tg (6/2) + do6/2] [2 cos y' + sin yf {tg (у' + 6/2) +
            //!+ tg(Y'-6/2)}].
	     __ATTR_ALWAYS_INLINE__
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     static inline
#pragma omp declare simd simdlen(16)
             float compute_Dmax(const float h,
	                        const float delta,
				const float d_ob,
				const float gamma) {
             
                  float Dmax = 0.0f;
		  float delta2,d_ob2,cosg,sing,t0,t1,t2,tant0,tant1;
		  float t3,t4,t5;
		  delta2 = 0.5f*delta;
                  d_ob2  = 0.5f*d_ob;
                  cosg   = std::cos(gamma);
                  if(delta<=gamma) {
                     t0  = h*delta+d_ob;
                     Dmax= t0/cosg;
                 }
                  t0     = gamma+delta2;
                  t1     = gamma-delta2;
                  sing   = std::sin(gamma);
                  tant1  = std::tan(t0);
                  tant2  = std::tan(t1);
                  t3     = h*std::tan(delta2)+d_ob2;
                  t4     = 2.0f*cosg+sing;
                  t5     = tant1+tant2;
                  Dmax   = t3*t4*t5;
		  return (Dmax);
	   }


	   
             __ATTR_ALWAYS_INLINE__
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     static inline
#pragma omp declare simd simdlen(8)
             double compute_Dmax(const double h,
	                         const double delta,
				 const double d_ob,
				 const double gamma) {
             
                  double Dmax = 0.0;
		  double delta2,d_ob2,cosg,sing,t0,t1,t2,tant0,tant1;
		  double t3,t4,t5;
		  delta2 = 0.5*delta;
                  d_ob2  = 0.5*d_ob;
                  cosg   = std::cos(gamma);
                  if(delta<=gamma) {
                     t0  = h*delta+d_ob;
                     Dmax= t0/cosg;
                 }
                  t0     = gamma+delta2;
                  t1     = gamma-delta2;
                  sing   = std::sin(gamma);
                  tant1  = std::tan(t0);
                  tant2  = std::tan(t1);
                  t3     = h*std::tan(delta2)+d_ob2;
                  t4     = 2.0*cosg+sing;
                  t5     = tant1+tant2;
                  Dmax   = t3*t4*t5;
		  return (Dmax);
	   }


	   //! Размер зеркала в направлении, перпендикулярном плоскости
           //! чертежа, приблизительно равен
           //! Formula 2, p. 58
             __ATTR_ALWAYS_INLINE__
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     static inline
#pragma omp declare simd simdlen(16)
             float compute_Dmin(const float h,
	                        const float delta,
				const float d_ob) {

                    float Dmin = 0.0f;
		    Dmin       = h*delta+d_ob;
		    return (Dmin);
	  }


	     __ATTR_ALWAYS_INLINE__
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     static inline
#pragma omp declare simd simdlen(8)
             double compute_Dmin(const double h,
	                         const double delta,
				 const double d_ob) {

                    double Dmin = 0.0;
		    Dmin       = h*delta+d_ob;
		    return (Dmin);
	  }



          //!Если зеркало осуществляет сканирование в пространстве
          //!изображений его размеры
          //! Formula 7, p. 58
             __ATTR_ALWAYS_INLINE__
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     static inline
#pragma omp declare simd simdlen(16)
             float Dmax_imag_scan(const float H,
	                          const float F,
				  const float B,
				  const float d_ob,
				  const float gamma,
				  const float psi,
				  const float phi,
				  const float d) {

                 float Dmax = 0.0f;
		 float t0,t1,t2,t3;
		 float cosg,sing,tanp1,tanp2,psi2,phi2;
		 psi2  = 0.5f*psi;
                 if(psi2<=gamma && B<=d) {
                     phi2 = 0.5f*phi;
                     t0   = (F+F)*std::tan(phi2);
                     t1   = (H/F)*d_ob;
                     t2   = std::sin(gamma);
                     Dmax = (t0+t1)*t2;
                 }
                 t0    = (H/F)*(d_ob-B)+B;
                 cosg  = std::cos(gamma);
                 tanp1 = gamma+psi2;
                 tanp2 = gamma-psi2;
                 sing  = std::sin(gamma);
                 t1    = 2.0f*cosg+sing;
                 t2    = std::tan(tanp1)+std::tan(tanp2);
                 t3    = 0.5f*t1*t2;
                 Dmax  = t0*t3;
		 return (Dmax);
	  }


	     __ATTR_ALWAYS_INLINE__
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     static inline
#pragma omp declare simd simdlen(8)
             double Dmax_imag_scan(const double H,
	                           const double F,
				   const double B,
				   const double d_ob,
				   const double gamma,
				   const double psi,
				   const double phi,
				   const double d) {

                 double Dmax = 0.0;
		 double t0,t1,t2,t3;
		 double cosg,sing,tanp1,tanp2,psi2,phi2;
		 psi2  = 0.5*psi;
                 if(psi2<=gamma && B<=d) {
                     phi2 = 0.5*phi;
                     t0   = (F+F)*std::tan(phi2);
                     t1   = (H/F)*d_ob;
                     t2   = std::sin(gamma);
                     Dmax = (t0+t1)*t2;
                 }
                 t0    = (H/F)*(d_ob-B)+B;
                 cosg  = std::cos(gamma);
                 tanp1 = gamma+psi2;
                 tanp2 = gamma-psi2;
                 sing  = std::sin(gamma);
                 t1    = 2.0*cosg+sing;
                 t2    = std::tan(tanp1)+std::tan(tanp2);
                 t3    = 0.5*t1*t2;
                 Dmax  = t0*t3;
		 return (Dmax);
	  }


	     __ATTR_ALWAYS_INLINE__
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     static inline
#pragma omp declare simd simdlen(16)
             float Dmin_imag_scan(const float H,
	                          const float F,
				  const float d_ob,
				  const float B) {

                   float Dmin = 0.0f;
		   float t0,t1;
		   t0   = H/F;
                   t1   = (d_ob-B)+B;
                   Dmin = t0*t1;
		   return (Dmin);
	   }


	     __ATTR_ALWAYS_INLINE__
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     static inline
#pragma omp declare simd simdlen(16)
             double Dmin_imag_scan(const double H,
	                           const double F,
				   const double d_ob,
				   const double B) {

                   double Dmin = 0.0;
		   double t0,t1;
		   t0   = H/F;
                   t1   = (d_ob-B)+B;
                   Dmin = t0*t1;
		   return (Dmin);
	   }
	  

    } //eos

}// gms


#endif /*__GMS_EOS_IR_SENSOR_HPP__*/
