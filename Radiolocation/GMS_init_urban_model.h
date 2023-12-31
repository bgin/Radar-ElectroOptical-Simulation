

#ifndef __GMS_INIT_URBAN_MODEL_H__
#define __GMS_INIT_URBAN_MODEL_H__ 311220231556


namespace file_info {

     const unsigned int GMS_INIT_URBAN_MODEL_MAJOR = 1;
     const unsigned int GMS_INIT_URBAN_MODEL_MINOR = 0;
     const unsigned int GMS_INIT_URBAN_MODEL_MICRO = 0;
     const unsigned int GMS_INIT_URBAN_MODEL_FULLVER =
       1000U*GMS_INIT_URBAN_MODEL_MAJOR+100U*GMS_INIT_URBAN_MODEL_MINOR+
       10U*GMS_INIT_URBAN_MODEL_MICRO;
     const char * const GMS_INIT_URBAN_MODEL_CREATION_DATE = "31-12-2023 15:56 +00200 (SUN 31 DEC 2023 GMT+2)";
     const char * const GMS_INIT_URBAN_MODEL_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_INIT_URBAN_MODEL_SYNOPSIS      = "Urban model data initialization routines.."

}

#include <cstdint>
#include "GMS_config.h"
#include "GMS_malloc.h"


namespace gms {

         namespace radiolocation {
       
       
                    /*
                          Initialize data type: BRCcount_t
                    */
                    __ATTR_COLD__
                    __ATTR_ALIGN__(32);
                    static
                    inline
                    void init_BRCcount(BRCcount_t & brc,
                                       const int32_t _nbpc,
                                       const int32_t _nbpr)  {
                         
                         using namespace gms::common;
                         brc.nbpc = _nbpc;
                         const std::size_t nbytes1 = sizeof(int32_t)*static_cast<std::size_t>(brc.nbpc);
                         brc.nbpr = _nbpr;
                         const std::size_t nbytes2 = sizeof(int32_t)*static_cast<std::size_t>(brc.nbpr);
                         brc.bpc = (int32_t*)gms_mm_malloc(nbytes,64ULL);
                         brc.bpr = (int32_t*)gms_mm_malloc(nbytes,64ULL);
                  }
                  
                  
                  /*
                         Destroy data type: BRCcount_t
                  */
                    __ATTR_COLD__
                    __ATTR_ALIGN__(32);
                    static
                    inline
                    void destroy_BRCcount(BRCcount_t & brc) {
                         using namespace gms::common;
                         if(__builtin_expect(brc.bpc!=NULL,1))
                              gms_mm_free(brc.bpc);
                         if(__builtin_expect(brc.bpr!=NULL,1))
                              gms_mm_free(brc.bpr); 
                    }
                                                       
                    /*
                          Initialize data type: BLatLondR1x_t
                    */  
                    __ATTR_COLD__
                    __ATTR_ALIGN__(32);
                    static
                    inline                                
                    template<typename T>
                    void init_BLatLondR1x(BLatLondR1x<T> &bll,
                                          const std::size_t _nblatd,
                                          const std::size_t _nblond) {
                          
                         bll.nblatd = _nblatd;
                         bll.nblond = _nblond;
                         bll.blatd  = DC1D<T>(bll.nblatd);
                         bll.blond  = DC1D<T>(bll.nblond);                      
                   }
                   
                   
                   /*
                        Initialize data type:  BLatLonrR1x_t 
                   */
                    __ATTR_COLD__
                    __ATTR_ALIGN__(32);
                    static
                    inline                                
                    template<typename T>
                    void init_BLatLondR1x(BLatLonrR1x_t<T> &bll,
                                          const std::size_t _nblatr,
                                          const std::size_t _nblonr) {
                         
                         bll.nblatr = _nblatr;
                         bll.nblonr = _nblonr;
                         bll.blatr  = DC1D<T>(bll.nblatr);
                         bll.blonr  = DC1D<T>(bll.nblonr);                      
                   }
                   
                   
                   /*
                         Initialize data type:  EllpbR1x_t
                   */
                    __ATTR_COLD__
                    __ATTR_ALIGN__(32);
                    static
                    inline                                
                    template<typename T>
                    void init_EllpbR1x(EllpbR1x_t &ell,
                                       const std::size_t _nellpb) {
                          
                          ell.nellpb = _nellpb;
                          ell.ellpb  = DC1D<T>(ell.nellpb);                   
                   }
                   
                   
                   /*
                         Initialize data type: PxybR1x_t
                   */
                    __ATTR_COLD__
                    __ATTR_ALIGN__(32);
                    static
                    inline                                
                    template<typename T>
                    void init_PxybR1x(PxybR1x_t<T> pxy,
                                      const std::size_t _nbpc,
                                      const std::size_t _npxb,
                                      const std::size_t _npyb) {
                         
                         pxy.nbpc = _nbpc;
                         pxy.npxb = _npxb;
                         pxy.npyb = _npyb;
                         pxy.pxb  = DC1D<T>(pxy.nbpc*pxy.npxb);
                         pxy.pyb  = DC1D<T>(pxy.nbpc*pxy.npyb);                  
                   }
            
     }

}

























#endif /*__GMS_INIT_URBAN_MODEL_H__*/
