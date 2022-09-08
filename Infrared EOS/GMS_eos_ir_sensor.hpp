

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
#include <limits>
#include "GMS_config.h"

namespace gms {

         namespace eos {



	      namespace {
                  // Taken from StackOverflow article.
		  // https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison
#pragma omp declare simd simdlen(16)
                  bool approximatelyEqual(const float a,
		                          const float b,
					  const float epsilon) {
			   const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return std::fabs(a - b) <=
                           ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
                   }
#pragma omp declare simd simdlen(8)
		   bool approximatelyEqual(const double a,
		                           const double b,
					   const double epsilon) {
			   const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return fabs(a - b) <=
                           ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
                   }
#pragma omp declare simd simdlen(16)
                  bool essentiallyEqual(const float a,
		                        const float b,
					const float epsilon) {
                           const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return fabs(a - b) <=
			   ((fabsa > fabsb ? fabsb : fabsa) * epsilon);
                   }
#pragma omp declare simd simdlen(8)
                   bool essentiallyEqual(const double a,
		                         const double b,
					 const double epsilon) {
                           const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return fabs(a - b) <=
			   ((fabsa > fabsb ? fabsb : fabsa) * epsilon);
                   }
#pragma omp declare simd simdlen(16)		   
                  bool definitelyGreaterThan(const float a,
		                             const float b,
					     const float epsilon) {
                           const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return (a - b) >
			   ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
		  }
#pragma omp declare simd simdlen(8)
		  bool definitelyGreaterThan(const double a,
		                             const double b,
					     const double epsilon) {
                           const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return (a - b) >
			   ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
		  }
#pragma omp declare simd simdlen(16)
                  bool definitelyLessThan(const float a,
		                          const float b,
					  const float epsilon) {
                           const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return (b - a) >
			   ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
                 }
#pragma omp declare simd simdlen(8)
		 bool definitelyLessThan( const double a,
		                          const double b,
					  const double epsilon) {
                           const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return (b - a) >
			   ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
                 }


	    }


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


	   //!величина расфокусировки
           //!Formula 1, p. 59
	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float defocus_cof(const float l2,
	                      const float alpha,
			      const float O,
			      const int32_t inf) {

                 float df = 0.0f;
		 float cos2a,icos;
		 cos2a = std::cos(alpha+alpha);
                 icos  = 1.0f/cos2a;
                 if(inf) 
                    df    = l2/(icos-1.0f)*O;
                 else
                    df    = l2/(icos-1.0f);
                 return (df);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double defocus_cof(const double l2,
	                       const double alpha,
			       const double O,
			       const int32_t inf) {

                 double df = 0.0;
		 double cos2a,icos;
		 cos2a = std::cos(alpha+alpha);
                 icos  = 1.0/cos2a;
                 if(inf) 
                    df    = l2/(icos-1.0)*O;
                 else
                    df    = l2/(icos-1.0);
                 return (df);
	  }


	  //! Диаметр кружка рассеяния р
          //! Formula 3, p.59
	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float circle_dispersion(const float d,
	                            const float l1,
				    const float l2,
				    const float alpha,
				    const float O,
				    const int32_t inf) {

                  float rho = 0.0f;
		  float t0,t1;
		  t0  = d/(l1+l2);
                  t1  = defocus_cof(l2,alpha,O,inf);
                  rho = t0*t1;
		  return (rho);
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double circle_dispersion(const double d,
	                             const double l1,
				     const double l2,
				     const double alpha,
				     const double O,
				     const int32_t inf) {

                  double rho = 0.0;
		  double t0,t1;
		  t0  = d/(l1+l2);
                  t1  = defocus_cof(l2,alpha,O,inf);
                  rho = t0*t1;
		  return (rho);
	   }


	   //!Formula 2, p. 59
	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float circ_dispers_diam(const float l1,
	                            const float l2,
				    const float alpha,
				    const float O,
				    const int32_t inf) {

                 float ratio = 0.0f;
		 float t0,t1;
		 t0    = l1+l2;
                 t1    = defocus_cos(l2,alpha,O,inf);
                 ratio = t1/t0;
		 return (ratio);
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double circ_dispers_diam(const double l1,
	                             const double l2,
				     const double alpha,
				     const double O,
				     const int32_t inf) {

                 double ratio = 0.0;
		 double t0,t1;
		 t0    = l1+l2;
                 t1    = defocus_cos(l2,alpha,O,inf);
                 ratio = t1/t0;
		 return (ratio);
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float defocus_small_ang(const float O,
	                            const float l2,
				    const float alpha) {

                float rho = 0.0f;
		float t0,t1,t2,alpha2;
		const float eps = std::numeric_limits<float>::epsilon();
		alpha2 = alpha+alpha;
		t0     = std::cos(alpha2);
		t1     = 1.0f-(alpha2*alpha2*0.5f);
		if(approximatelyEqual(t0,t1,eps)) {
                   t2  = l2*0.5f;
		   rho = O*t2*alpha2*alpha2;
		   return (rho);
		}
		else {
                   rho = std::numeric_limits<float>::quiet_NaN();
		   return (rho);
		}
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double defocus_small_ang(const double O,
	                             const double l2,
				     const double alpha) {

                double rho = 0.0;
		double t0,t1,t2,alpha2;
		const double eps = std::numeric_limits<double>::epsilon();
		alpha2 = alpha+alpha;
		t0     = std::cos(alpha2);
		t1     = 1.0-(alpha2*alpha2*0.5);
		if(approximatelyEqual(t0,t1,eps)) {
                   t2  = l2*0.5;
		   rho = O*t2*alpha2*alpha2;
		   return (rho);
		}
		else {
                   rho = std::numeric_limits<double>::quiet_NaN();
		   return (rho);
		}
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float traj_scan_dxdt(const float dx[2],
	                         const float dt[2]) {
 
                  float dxdt = 0.0f;
		  dxdt = (dx[1]-dx[0])/(dt[1]-dt[0]);
		  return (dxdt);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double traj_scan_dxdt(const double dx[2],
	                         const double dt[2]) {
 
                  double dxdt = 0.0;
		  dxdt = (dx[1]-dx[0])/(dt[1]-dt[0]);
		  return (dxdt);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float traj_scan_dydt(const float dy[2],
	                         const float dt[2]) {
 
                  float dydt = 0.0f;
		  dxdt = (dy[1]-dy[0])/(dt[1]-dt[0]);
		  return (dydt);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double traj_scan_dydt(const double dy[2],
	                          const double dt[2]) {
 
                  double dydt = 0.0;
		  dxdt = (dy[1]-dy[0])/(dt[1]-dt[0]);
		  return (dydt);
	  }


	  //! СКАНИРОВАНИЕ ЗЕРКАЛОМ, ВРАЩАЮЩИМСЯ
          //! ВОКРУГ ОСИ, НЕПЕРПЕНДИКУЛЯРНОЙ К НЕМУ
          //! Formula 1, p. 100
            __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
	    float fov_x_axis(const float H,
	                     const float delta,
			     const float gamma) {

               float ax = 0.0f;
	       float gamm2,tdel;
	       gamm2    = 0.5f*gamma;
	       tdel     = std::tan(delta);
	       ax       = H*tdel*std::cos(gamm2);
	       return (ax);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
	    double fov_x_axis(const double H,
	                      const double delta,
			      const double gamma) {

               double ax = 0.0;
	       double gamm2,tdel;
	       gamm2    = 0.5*gamma;
	       tdel     = std::tan(delta);
	       ax       = H*tdel*std::cos(gamm2);
	       return (ax);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
	    float fov_y_axis(const float H,
	                     const float delta,
			     const float gamma) {

               float ay = 0.0f;
	       float ax,t0;
	       t0       = 0.5f*gamma;
	       ax       = fov_x_axis(H,delta,gamma);
	       ay       = t0*ax;
	       return (ay);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
	    double fov_y_axis(const double H,
	                      const double delta,
			      const double gamma) {

               double ay = 0.0;
	       double ax,t0;
	       t0       = 0.5*gamma;
	       ax       = fov_x_axis(H,delta,gamma);
	       ay       = t0*ax;
	       return (ay);
	  }


	 //!Если рабочая зона сканирования ограничена углом G, то
         //!ширина захвата
         //!Formula 3, p. 100
            __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float scan_width(const float H,
	                     const float gamma,
			     const float theta) {

                float B = 0.0f;
		float gam2,th2,t0,t1;
		gam2  = 0.5f*gamma;
                th2   = 0.5f*theta;
                t0    = std::tan(gam2);
                t1    = std::sin(th2);
                B     = (H+H)*t0*t1;
		return (B);
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double scan_width(const double H,
	                      const double gamma,
			      const double theta) {

                double B = 0.0;
		double gam2,th2,t0,t1;
		gam2  = 0.5*gamma;
                th2   = 0.5*theta;
                t0    = std::tan(gam2);
                t1    = std::sin(th2);
                B     = (H+H)*t0*t1;
		return (B);
	   }


	   //!Плоскопараллельная пластинка, установленная за 
           //!объективом, изменяет ход лучей таким образом, что изображение
           //! светящейся точки отодвигается и его положение зависит от угла у
           //!между оптической осью и нормалью N к поверхности пластинки
           //! Formula 7,8 p. 106
	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float refract_shift(const float i1,
	                        const float delta,
				const float alfa,
				const float gamma,
				const float n) {

                 float l = 0.0f;
		 float ag,num,den,sin2,sag,t0,t1;
		 const float eps = std::numeric_limits<float>::epsilon();
		 if(approximatelyEqual(i1,ag,eps)) {
                     sag  = std::sin(ag);
                     t0   = delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0f-sag;
                     den  = n*n-sag;
                     t1   = 1.0f-std::sqrt(num/den);
                     l    = t0*t2;
		 }
		 else if(alfa==0.0f) {
                     sag  = std::sin(gamma);
                     t0   = -delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0f-sin2;
                     den  = n*n-sin2;
                     t1   = 1.0f-sqrt(num/den);
                     l    = t0*t1;
		 }
		 else {
                     sag  = std::sin(i1);
                     t0   = delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0f-sin2;
                     den  = n*n-sin2;
                     t1   = 1.0f-sqrt(num/den);
                     l    = t0*t1;
		 }
		 return (l);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double refract_shift(const double i1,
	                         const double delta,
				 const double alfa,
				 const double gamma,
				 const double n) {

                 double l = 0.0;
		 double ag,num,den,sin2,sag,t0,t1;
		 const double eps = std::numeric_limits<double>::epsilon();
		 if(approximatelyEqual(i1,ag,eps)) {
                     sag  = std::sin(ag);
                     t0   = delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0-sag;
                     den  = n*n-sag;
                     t1   = 1.0-std::sqrt(num/den);
                     l    = t0*t2;
		 }
		 else if(alfa==0.0) {
                     sag  = std::sin(gamma);
                     t0   = -delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0-sin2;
                     den  = n*n-sin2;
                     t1   = 1.0-sqrt(num/den);
                     l    = t0*t1;
		 }
		 else {
                     sag  = std::sin(i1);
                     t0   = delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0-sin2;
                     den  = n*n-sin2;
                     t1   = 1.0-sqrt(num/den);
                     l    = t0*t1;
		 }
		 return (l);
	  }


	  //!Formula 1, p. 108
	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            void project_xy_axis( const float l,
	                          const float alpha,
				  float &xl,
				  float &yl) {

                float absl = 0.0f;
		absl       = std::abs(l);
		xl         = std::cos(alpha);
		yl         = std::sin(alpha);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            void project_xy_axis( const double l,
	                          const double alpha,
				  double &xl,
				  double &yl) {

                double absl = 0.0;
		absl       = std::abs(l);
		xl         = std::cos(alpha);
		yl         = std::sin(alpha);
	  }


	 //!Величину смещения луча s вдоль перпендикуляра к 
         //!поверхности пластинки
         //! Formula 2, p. 108
            __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float s_shift(const float l,
	                  const float alpha,
		          const float gamma) {

              float s = 0.0f;
	      float ag,sag;
	      ag = alpha-gamma;
              sag= std::sin(ag);
              s  = l/sag;
	      return (s);
	 }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double s_shift(const double l,
	                   const double alpha,
		           const double gamma) {

              double s = 0.0;
	      double ag,sag;
	      ag = alpha-gamma;
              sag= std::sin(ag);
              s  = l/sag;
	      return (s);
	 }


	 //! Проекции s на оси координат равны
         //! Formula 4, p. 108
            __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            void  project_s_xy(const float s,
	                       const float gamma,
			       float &xs,
			       float &ys) {

                 xs = s*std::cos(gamma);
                 ys = s*std::sin(gamma);
	 }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            void  project_s_xy(const double s,
	                       const double gamma,
			       double &xs,
			       double &ys) {

                 xs = s*std::cos(gamma);
                 ys = s*std::sin(gamma);
	 }


      //! что расстояния от начала координат О до точек
      //! пересечения лучей, образующих с горизонталью угла ±а, с 
      //! перпендикуляром к пластинке
      //! Formula 1, p. 110
            __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float ray_intercept_pa(const float delta,
	                           const float alpha,
				   const float gamma,
				   const float n) {

               float sp = 0.0f;
	       float ag,num,den,n2,sag,sin2;
	       ag  = abs(alpha)-gamma;
               n2  = n*n;
               num = std::cos(ag);
               sag = std::sin(ag);
               sin2= std::sag*sag;
               den = std::sqrt(n2-sin2);
               sp  = std::delta*1.0f-(num/den);
	       return (sp);
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double ray_intercept_pa(const double delta,
	                            const double alpha,
				    const double gamma,
				    const double n) {

               double sp = 0.0;
	       double ag,num,den,n2,sag,sin2;
	       ag  = abs(alpha)-gamma;
               n2  = n*n;
               num = std::cos(ag);
               sag = std::sin(ag);
               sin2= std::sag*sag;
               den = std::sqrt(n2-sin2);
               sp  = std::delta*1.0-(num/den);
	       return (sp);
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float ray_intercept_na( const float delta,
	                            const float alpha,
				    const float gamma,
				    const float n) {

               float sn = 0.0f;
	       float ag,num,den,n2,sag,sin2;
	       ag  = abs(alpha)+gamma;
               n2  = n*n;
               num = std::cos(ag);
               sag = std::sin(ag);
               sin2= sag*sag;
               den = sqrt(n2-sin2);
               sn  = delta*1.0f-(num/den);
	       return (sn);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double ray_intercept_na( const double delta,
	                             const double alpha,
				     const double gamma,
				     const double n) {

               double sn = 0.0;
	       double ag,num,den,n2,sag,sin2;
	       ag  = abs(alpha)+gamma;
               n2  = n*n;
               num = std::cos(ag);
               sag = std::sin(ag);
               sin2= sag*sag;
               den = sqrt(n2-sin2);
               sn  = delta*1.0-(num/den);
	       return (sn);
	  }


	  //! Formula 3, p. 110
	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float ray_diff(const float delta,
	                   const float alpha,
			   const float gamma,
			   const float n,
			   const float u) {

              float ds = 0.0f;
	      float t0,t1,u2,u2g,su2,sg,t2,n2,t3,t4,t5;
	      n   = n*n;
              u2  = u*0.5f;
              u2g = u2-gamma;
              t2  = std::sin(u2g);
              su2 = t2*t2;
              if(n2>=su2){
                 t3 = (-2.0f*delta)/n;
                 t4 = std::sin(u2);
                 t5 = std::sin(gamma);
                 ds = t3*t4*t5;
	      }
              else {
                 t0  = ray_intercept_pa(delta,alpha,gamma,n)
                 t1  = ray_intercept_na(delta,alpha,gamma,n)
                 ds  = t0-t1;
             }
	     return (ds);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double ray_diff(const double delta,
	                    const double alpha,
			    const double gamma,
			    const double n,
			    const double u) {

              double ds = 0.0;
	      double t0,t1,u2,u2g,su2,sg,t2,n2,t3,t4,t5;
	      n   = n*n;
              u2  = u*0.5;
              u2g = u2-gamma;
              t2  = std::sin(u2g);
              su2 = t2*t2;
              if(n2>=su2){
                 t3 = (-2.0*delta)/n;
                 t4 = std::sin(u2);
                 t5 = std::sin(gamma);
                 ds = t3*t4*t5;
	      }
              else {
                 t0  = ray_intercept_pa(delta,alpha,gamma,n)
                 t1  = ray_intercept_na(delta,alpha,gamma,n)
                 ds  = t0-t1;
             }
	     return (ds);
	  }


	  //!Поле точек пересечения лучей, преломленных пластинкой,
          //!относительно оси Ох (рис. 87) имеет симметрию, поэтому 
          //!упростим обозначения и выполним расчет соответствующих 
          //!координат на основании
          //! Formula 6,7, p. 111
	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            void compute_dxdy(const float alpha,
	                      const float beta,
			      const float delta,
			      const float gamma,
			      const float n,
			      const float u,
			      float &dx,
			      float &dy) {

                 float ag,ds,t0,t1,t2;
		 ag  = alpha+gamma;
                 ds  = ray_diff(delta,alfa,gamma,n,u);
                 t0  = sin(ag);
                 t1  = 2.0f*std::sin(alpha);
                 t2  = 2.0f*std::cos(alpha);
                 dx  = t0/t1*ds;
                 dy  = t0/t2*ds;
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            void compute_dxdy(const double alpha,
	                      const double beta,
			      const double delta,
			      const double gamma,
			      const double n,
			      const double u,
			      double &dx,
			      double &dy) {

                 double ag,ds,t0,t1,t2;
		 ag  = alpha+gamma;
                 ds  = ray_diff(delta,alfa,gamma,n,u);
                 t0  = sin(ag);
                 t1  = 2.0*std::sin(alpha);
                 t2  = 2.0*std::cos(alpha);
                 dx  = t0/t1*ds;
                 dy  = t0/t2*ds;
	   }


	   //! Formula 7,8  p. 111
	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            void compute_xy(  const float alpha,
	                      const float beta,
			      const float delta,
			      const float gamma,
			      const float n,
			      const float u,
			      float &x,
			      float &y) {

                float sag,cag,pa,dx,dy,xs,ys;
		sag  = std::sin(gamma);
                cag  = std::cos(gamma);
                pa   = ray_intercept_pa(delta,alpha,gamma,n);
                xs   = pa*sag;
                ys   = pa*cag;
                compute_dxdy(alpha,beta,delta,gamma,n,u,dx,dy);
                x    = xs+dx;
                y    = ys+dx;
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            void compute_xy(  const double alpha,
	                      const double beta,
			      const double delta,
			      const double gamma,
			      const double n,
			      const double u,
			      double &x,
			      double &y) {

                double sag,cag,pa,dx,dy,xs,ys;
		sag  = std::sin(gamma);
                cag  = std::cos(gamma);
                pa   = ray_intercept_pa(delta,alpha,gamma,n);
                xs   = pa*sag;
                ys   = pa*cag;
                compute_dxdy(alpha,beta,delta,gamma,n,u,dx,dy);
                x    = xs+dx;
                y    = ys+dx;
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
	    void compute_xdyd(const float gamma,
	                      const float u,
			      const float n,
			      float &xd,
			      float &yd) {

                float cosg,sing,sin2s,sin2d,u2,u2gs,u2gd,n2,t0,t1,t2,t3,t4;
		cosg = std::cos(gamma);
                sing = std::sin(gamma);
                u2   = u*0.5f;
                n2   = n*n;
                u2gs = u2+gamma;
                ungd = u2-gamma;
                t0   = sin(u2gs);
                sin2s= t0*t0;
                t1   = std::sin(u2gd);
                sin2d= t1*t1;
                t2   = 1.0f/(4.0f*std::sin(u2));
                t3   = std::sqrt(n2-sin2s);
                t0   = sin2s/t3;
                t4   = std::sqrt(n2-sin2d);
                t1   = sin2d/t4;
                dx   = cosg-t2*(t0+t1);
                t2   = 1.0f/(4.0f*std::cos(u2));
                dy   = sing-t2*(t0-t1);
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
	    void compute_xdyd(const double gamma,
	                      const double u,
			      const double n,
			      double &xd,
			      double &yd) {

                double cosg,sing,sin2s,sin2d,u2,u2gs,u2gd,n2,t0,t1,t2,t3,t4;
		cosg = std::cos(gamma);
                sing = std::sin(gamma);
                u2   = u*0.5;
                n2   = n*n;
                u2gs = u2+gamma;
                ungd = u2-gamma;
                t0   = sin(u2gs);
                sin2s= t0*t0;
                t1   = std::sin(u2gd);
                sin2d= t1*t1;
                t2   = 1.0f/(4.0*std::sin(u2));
                t3   = std::sqrt(n2-sin2s);
                t0   = sin2s/t3;
                t4   = std::sqrt(n2-sin2d);
                t1   = sin2d/t4;
                dx   = cosg-t2*(t0+t1);
                t2   = 1.0/(4.0f*std::cos(u2));
                dy   = sing-t2*(t0-t1);
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            void paraxial_xdyd(const float gamma,
	                       const float alpha,
			       const float n,
			       float &xd,
			       float &yd) {

               float n2,cosg,sing,sin4g,sin2g,num,den,cos2g,t0,t1,n2ss;
	       n2    = n*n;
               cosg  = std::cos(gamma);
	       n2ss  = n2-sing*sing;
               cos2g = std::cos(gamma+gamma);
               sing  = std::sin(gamma);
               sin4g = sing*sing*sing*sing;
               num   = n2*cos2g+sin4g;
               den   = std::pow(n2ss,1.5f);
               xd    = cosg-num/den;
               t0    = sqrt(n2ss);
               t1    = 1.0f-cosg/t0;
               yd    = sing*t1;
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            void paraxial_xdyd(const double gamma,
	                       const double alpha,
			       const double n,
			       double &xd,
			       double &yd) {

               double n2,cosg,sing,sin4g,sin2g,num,den,cos2g,t0,t1,n2ss;
	       n2    = n*n;
               cosg  = std::cos(gamma);
	       n2ss  = n2-sing*sing;
               cos2g = std::cos(gamma+gamma);
               sing  = std::sin(gamma);
               sin4g = sing*sing*sing*sing;
               num   = n2*cos2g+sin4g;
               den   = std::pow(n2ss,1.5);
               xd    = cosg-num/den;
               t0    = sqrt(n2ss);
               t1    = 1.0-cosg/t0;
               yd    = sing*t1;
	  }


	  //!СКАНИРОВАНИЕ ВРАЩАЮЩИМИСЯ ОБЪЕКТИВАМИ
          //!Formula 1, p. 121
            __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            void fov_axay(const float H,
	                  const float delx,
			  const float dely,
			  const float phi,
			  float &ax,
			  float &ay) {

               float sec2,phi2,t0,t1,sec;
	       phi2  = 0.5f*phi;
               sec   = 1.0f/std::cos(phi2);
               sec2  = sec*sec;
               ax    = H*delx*sec2;
               ay    = H*dely*sec;
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            void fov_axay(const double H,
	                  const doube delx,
			  const double dely,
			  const double phi,
			  double &ax,
			  double &ay) {

               double sec2,phi2,t0,t1,sec;
	       phi2  = 0.5*phi;
               sec   = 1.0/std::cos(phi2);
               sec2  = sec*sec;
               ax    = H*delx*sec2;
               ay    = H*dely*sec;
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            void fov_dxdy(const float x,
	                  const float y,
			  const float F,
			  const float phi,
			  float &dx,
			  float &dy) {

              float d0x,d0y,phi2;
	      d0y   = y/F;
              phi2  = 0.5f*phi;
              d0x   = x/F;
              dy    = d0y;
              dx    = d0x*std::cos(phi2);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            void fov_dxdy(const double x,
	                  const double y,
			  const double F,
			  const double phi,
			  double &dx,
			  double &dy) {

              double d0x,d0y,phi2;
	      d0y   = y/F;
              phi2  = 0.5*phi;
              d0x   = x/F;
              dy    = d0y;
              dx    = d0x*std::cos(phi2);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            void volt_impulse_uxuy(const float u,
	                           const float om1,
				   const float om2,
				   const float t,
				   float &ux,
				   float &uy) {

                float om1t,om2t,t0,t1;
		om1t = om1*t;
                om2t = om2*t;
                t0   = std::sin(om1t)+std::sin(om2t);
                t1   = std::cos(om1t)+std::cos(om2t);
                ux   = u*t0;
                uy   = u*t1;
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            void volt_impulse_uxuy(const double u,
	                           const double om1,
				   const double om2,
				   const double t,
				   double &ux,
				   double &uy) {

                double om1t,om2t,t0,t1;
		om1t = om1*t;
                om2t = om2*t;
                t0   = std::sin(om1t)+std::sin(om2t);
                t1   = std::cos(om1t)+std::cos(om2t);
                ux   = u*t0;
                uy   = u*t1;
	   }


	//! Phase Modulation
        //! Formula 1, p. 143
        //! растрового анализатора со 
        //!скрещивающимися осями, выполненного в виде надетой на вращающийся 
        //!барабан тонкой пленки, прозрачность которой изменяется по 
        //!синусоидальному закону
	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(16)
            float raster_transparency(const float rho_avg,
	                              const float rho_max,
				      const float rho_min,
				      const float l,
				      const float L,
				      const float N) {

               float rho = 0.0f;
	       constexpr float twopi = 6.283185307179586476925286766559f;
	       float t0,t1,t2;
	       t0  = 0.5f*(rho_max-rho_min);
               t1  = L/N;
               t2  = std::sin(twopi*l*t1);
               rho = rho_avg+t0*t2;
	       return (rho);
	  }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
#pragma omp declare simd simdlen(8)
            double raster_transparency(const double rho_avg,
	                               const double rho_max,
				       const double rho_min,
				       const double l,
				       const double L,
				       const double N) {

               double rho = 0.0;
	       constexpr double twopi = 6.283185307179586476925286766559;
	       double t0,t1,t2;
	       t0  = 0.5*(rho_max-rho_min);
               t1  = L/N;
               t2  = std::sin(twopi*l*t1);
               rho = rho_avg+t0*t2;
	       return (rho);
	  }


#include "GMS_avint.hpp"
 

             //!СТРУКТУРА И СПЕКТР МОДУЛИРОВАННОГО ПОТОКА
             //!ИЗЛУЧЕНИЯ
             //!Formula 1, p. 178
             //! Ф(*) = Int rp(z,t)E(z,t) dsig
            __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
	    void raster_flux_integral_omp(const float * __restrict rhoE,
	                                  const float * __restrict absc,
					  const int32_t n,
					  const int32_t t,
					  const float * __restrict xlo,
					  const float * __restrict xup,
					  float * __restrict Phit,
					  float * __restrict ier) {

                 float ans_x   = 0.0f;
		 int32_t err_x = 0;
#pragma omp parallel for default(none) schedule(runtime) \
                 private(i,ans,err) shared(t,rhoE,absc,n,xlo,xup)
		 for(int32_t i = 0, i < t; ++i) {
                     avint(&rhoE[i*n],&absc[0],n,xlo,xup,ans,err);
		     Phit[i] = ans;
		     ier[i]  = err;
		 }
	   }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
	    void raster_flux_integral_omp(const double * __restrict rhoE, // points to memory of size: (0:n-1,t)
	                                  const double * __restrict absc,
					  const int32_t n,
					  const int32_t t,
					  const double * __restrict xlo,
					  const double * __restrict xup,
					  double * __restrict Phit, // results data size: 0:t-1
					  double * __restrict ier) {

                 double ans_x   = 0.0f;
		 int32_t err_x = 0;
#pragma omp parallel for default(none) schedule(runtime) \
                 private(i,ans,err) shared(t,rhoE,absc,n,xlo,xup)
		 for(int32_t i = 0, i < t; ++i) {
                     avint(&rhoE[i*n],&absc[0],n,xlo,xup,ans,err);
		     Phit[i] = ans;
		     ier[i]  = err;
		 }
	   }






    } //eos

}// gms


#endif /*__GMS_EOS_IR_SENSOR_HPP__*/
