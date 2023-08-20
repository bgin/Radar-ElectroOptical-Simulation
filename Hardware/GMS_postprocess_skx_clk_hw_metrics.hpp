

#ifndef __GMS_POSTPROCESS_SKX_CLK_HW_METRICS_HPP__
#define __GMS_POSTPROCESS_SKX_CLK_HW_METRICS_HPP__ 150820230813


/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

namespace file_version {

    const unsigned int GMS_POSTPROCESS_SKX_CLK_HW_METRICS_MAJOR = 1U;
    const unsigned int GMS_POSTPROCESS_SKX_CLK_HW_METRICS_MINOR = 0U;
    const unsigned int GMS_POSTPROCESS_SKX_CLK_HW_METRICS_MICRO = 0U;
    const unsigned int GMS_POSTPROCESS_SKX_CLK_HW_METRICS_FULLVER =
      1000U*GMS_POSTPROCESS_SKX_CLK_HW_METRICS_MAJOR+
      100U*GMS_POSTPROCESS_SKX_CLK_HW_METRICS_MINOR+
      10U*GMS_POSTPROCESS_SKX_CLK_HW_METRICS_MICRO;
    const char * const GMS_POSTPROCESS_SKX_CLK_HW_METRICS_CREATION_DATE = "15-08-2023 08:13 PM +00200 (TUE 15 AUG 2023 GMT+2)";
    const char * const GMS_POSTPROCESS_SKX_CLK_HW_METRICS_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_POSTPROCESS_SKX_CLK_HW_METRICS_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_POSTPROCESS_SKX_CLK_HW_METRICS_DESCRIPTION   = "Postprocessing (time-series analysis) of SKX and CLK HW metrics";

}

#include <cstdint>
#include "GMS_config.h"
#include "GMS_malloc.h"
#include "GMS_preprocess_skx_clk_hw_metrics.h"
#include "GMS_cpu_perf_time_series_analysis.h"
#include "GMS_pmc_events_descriptive_stats.hpp"


namespace gms {

  
        
                             
       
        
        /*
            For the two arguments metrics.
        */
        template<int32_t len, int32_t lagh>
        struct   SKX_HW_metric_2_t __ATTR_ALIGN__(64) {
        
                using fptr2 = void(*)(const double __restrict *,
                             const double __restrict *,
                             double       __restrict *,
                             const int32_t);  
                                      
                 __ATTR_ALIGN__(8) double * __restrict m_param1;
                 __ATTR_ALIGN__(8) double * __restrict m_param2;
                 __ATTR_ALIGN__(8) double * __restrict m_samples;
                 __ATTR_ALIGN__(8) fptr2               m_pfmetric;    
                 char                                  m_metric[64];       
                 
                 SKX_HW_metric_2_t() noexcept(true) {
                      
                      m_param1     = nullptr;
                      m_param2     = nullptr;
                      m_samples    = nullptr;
                      m_m_pfmetric = nullptr;
                      m_metric[64] = {};
               }
               
               SKX_HW_metric_2_t(  const fptr2 pfmetric,
                                   const char[64] metric) noexcept(false) {
                      
                      using namespace gms::common;
                      const std::size_t samp_len = (std::size_t)len;
                      const std::size_t align    = (std::size_t)64;
                      m_param1    = (double*)gms_mm_malloc(samp_len,align);
                      m_param2    = (double*)gms_mm_malloc(samp_len,align);
                      m_samples   = (double*)gms_mm_malloc(samp_len,align);                                                   
                      m_pfmetric  = pfmetric;
                      _mm512_storeu_si512((char*)&m_metric[0],
                            _mm512_loadu_si512((const char*)&metric[0]));                                               
                    
               }
               
               SKX_HW_metric_2_t(const SKX_HW_metric_2_t &) = delete;
               
               SKX_HW_metric_2_t(SKX_HW_metric_2_t &&) = delete
               
               ~SKX_HW_metric_2_t() {
                     
                       using namespace gms::common;
                       gms_mm_free(m_param1);
                       m_param1   = nullptr;
                       gms_mm_free(m_param2);
                       m_param2   = nullptr;
                       gms_mm_free(m_samples);
                       m_samples  = nullptr;
                       m_pfmetric = nullptr;
               }     
               
               
               SKX_HW_metric_2_t &
               operator=(const SKX_HW_metric_2_t &) = delete;
               
               SKX_HW_metric_2_t &
               operator=(SKX_HW_metric_2_t &&) = delete;  
               
               void compute_HW_metric() {
                    
                    m_pfmetric(m_param1,m_param2,m_samples,len);
               }     
               
               void analyze_metric_canarm(const char * __restrict fname,
                                          const char * __restrict data_type,
                                          const bool use_omp) {
                   
                    cpu_perf_time_series_canarm<len,lagh>(m_samples,fname,
                                                          data_type,use_omp);
                                                                            
             }  
             
             void analyze_metric_unimar(const char * __restrict fname,
                                        const char * __restrict data_type) {
                 
                    cpu_perf_time_series_unimar<len,lagh>(m_samples,
                                                          fname,data_type);                           
            } 
            
            void analyze_metric_unibar(const char * __restrict fname,
                                       const char * __restrict data_type) {
                                       
                    cpu_perf_time_series_unibar<len,lagh>(m_samples,
                                                          fname,data_type);                           
           }    
           
           void analyze_metric_exsar(const char * __restrict fname,
                                     const char * __restrict data_type) {
                                     
                    cpu_perf_time_series_exsar<len,lagh>(m_samples,
                                                         fname,data_type);                       
           }
           
           void analyze_metric_bispec(const char * __restrict fname,
                                      const char * __restrict data_type,
                                      const bool use_omp) {
                                      
                    cpu_perf_time_series_bispec<len,lagh>(m_samples,
                                                          fname,data_type,use_omp);                          
           }
           
           void analyze_metric_thirmo(const char * __restrict fname,
                                      const char * __restrict data_type) {
                                      
                    cpu_perf_time_series_thirmo<len,lagh>(m_samples,
                                                          fname,data_type);                
           }
           
           void analyze_metric_autocor(const char * __restrict fname,
                                       const char * __restrict data_type) {
                                       
                    cpu_perf_time_series_autocor<len,lagh>(m_samples,
                                                           fname,data_type);                            
          }
          
       
               
      };
      
      
     /////////////////////////////////////////////////////////////////////////
     
     /*
         For the three arguments metrics.
     */
     
        template<int32_t len, int32_t lagh>
        struct   SKX_HW_metric_3_t __ATTR_ALIGN__(64) {
        
                using fptr3 = void(*)(const double __restrict *,
                                      const double __restrict *,
                                      const double __restrict *,
                                      double       __restrict *,
                                      const int32_t);  
                                      
                 __ATTR_ALIGN__(8) double * __restrict m_param1;
                 __ATTR_ALIGN__(8) double * __restrict m_param2;
                 __ATTR_ALIGN__(8) double * __restrict m_param3;
                 __ATTR_ALIGN__(8) double * __restrict m_samples;
                 __ATTR_ALIGN__(8) fptr3               m_pfmetric;    
                 char                                  m_metric[64];       
                 
                 SKX_HW_metric_3_t() noexcept(true) {
                      
                      m_param1     = nullptr;
                      m_param2     = nullptr;
                      m_param3     = nullptr;
                      m_samples    = nullptr;
                      m_m_pfmetric = nullptr;
                      m_metric[64] = {};
               }
               
               SKX_HW_metric_3_t(  const fptr3 pfmetric,
                                   const char[64] metric) noexcept(false) {
                      
                      using namespace gms::common;
                      const std::size_t samp_len = (std::size_t)len;
                      const std::size_t align    = (std::size_t)64;
                      m_param1    = (double*)gms_mm_malloc(samp_len,align);
                      m_param2    = (double*)gms_mm_malloc(samp_len,align);
                      m_param3    = (double*)gms_mm_malloc(samp_len,align);
                      m_samples   = (double*)gms_mm_malloc(samp_len,align);                                                   
                      m_pfmetric  = pfmetric;
                      _mm512_storeu_si512((char*)&m_metric[0],
                            _mm512_loadu_si512((const char*)&metric[0]));                                               
                    
               }
               
               SKX_HW_metric_3_t(const SKX_HW_metric_3_t &) = delete;
               
               SKX_HW_metric_3_t(SKX_HW_metric_3_t &&) = delete
               
               ~SKX_HW_metric_3_t() {
                     
                       using namespace gms::common;
                       gms_mm_free(m_param1);
                       m_param1   = nullptr;
                       gms_mm_free(m_param2);
                       m_param2   = nullptr;
                       gms_mm_free(m_param3);
                       m_param3   = nullptr;
                       gms_mm_free(m_samples);
                       m_samples  = nullptr;
                       m_pfmetric = nullptr;
               }     
               
               
               SKX_HW_metric_3_t &
               operator=(const SKX_HW_metric_3_t &) = delete;
               
               SKX_HW_metric_3_t &
               operator=(SKX_HW_metric_3_t &&) = delete;  
               
               void compute_HW_metric() {
                    
                    m_pfmetric(m_param1,m_param2,
                               m_param3,m_samples,len);
               }     
               
               void analyze_metric_canarm(const char * __restrict fname,
                                          const char * __restrict data_type,
                                          const bool use_omp) {
                   
                    cpu_perf_time_series_canarm<len,lagh>(m_samples,fname,
                                                          data_type,use_omp);
                                                                            
             }  
             
             void analyze_metric_unimar(const char * __restrict fname,
                                        const char * __restrict data_type) {
                 
                    cpu_perf_time_series_unimar<len,lagh>(m_samples,
                                                          fname,data_type);                           
            } 
            
            void analyze_metric_unibar(const char * __restrict fname,
                                       const char * __restrict data_type) {
                                       
                    cpu_perf_time_series_unibar<len,lagh>(m_samples,
                                                          fname,data_type);                           
           }    
           
           void analyze_metric_exsar(const char * __restrict fname,
                                     const char * __restrict data_type) {
                                     
                    cpu_perf_time_series_exsar<len,lagh>(m_samples,
                                                         fname,data_type);                       
           }
           
           void analyze_metric_bispec(const char * __restrict fname,
                                      const char * __restrict data_type,
                                      const bool use_omp) {
                                      
                    cpu_perf_time_series_bispec<len,lagh>(m_samples,
                                                          fname,data_type,use_omp);                          
           }
           
           void analyze_metric_thirmo(const char * __restrict fname,
                                      const char * __restrict data_type) {
                                      
                    cpu_perf_time_series_thirmo<len,lagh>(m_samples,
                                                          fname,data_type);                
           }
           
           void analyze_metric_autocor(const char * __restrict fname,
                                       const char * __restrict data_type) {
                                       
                    cpu_perf_time_series_autocor<len,lagh>(m_samples,
                                                           fname,data_type);                            
          }
          
                      
      };
      
      
      
   
    
      
      
      
      

} // gms

















#endif /*__GMS_POSTPROCESS_SKX_CLK_HW_METRICS_HPP__*/
