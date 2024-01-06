

#ifndef __GMS_URBAN_ADT_STAT_H__
#define __GMS_URBAN_ADT_STAT_H__


namespace file_info {

     const unsigned int GMS_URBAN_ADT_STAT_MAJOR = 1;
     const unsigned int GMS_URBAN_ADT_STAT_MINOR = 0;
     const unsigned int GMS_URBAN_ADT_STAT_MICRO = 0;
     const unsigned int GMS_URBAN_ADT_STAT_FULLVER =
       1000U*GMS_URBAN_ADT_STAT_MAJOR+100U*GMS_URBAN_ADT_STAT_MINOR+
       10U*GMS_URBAN_ADT_STAT_MICRO;
     const char * const GMS_URBAN_ADT_STAT_CREATION_DATE = "23-12-2023 12:09 +00200 (SAT 23 DEC 2023 GMT+2)";
     const char * const GMS_URBAN_ADT_STAT_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_URBAN_ADT_STAT_SYNOPSIS      = "Abstract data types (static allocation) representing built-up area for diffraction and scattering calculations."

}

#include <cstdint>
#include <memory>
#include "GMS_config"


namespace gms {

       namespace radiolocation {
    
    
            // Building units per column
            // Building units per row
            template<int32_t nbpc,int32_t nbpr>
            struct BRCcount_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64) int32_t bpc[nbpc];
                   __ATTR_ALIGN__(64) int32_t bpr[nbpr];
#elif defined (__AVX__) || defined (__AVX2__)
                   __ATTR_ALIGN__(32) int32_t bpc[nbpc];
                   __ATTR_ALIGN__(32) int32_t bpr[nbpr];
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16) int32_t bpc[nbpc];
                   __ATTR_ALIGN__(16) int32_t bpr[nbpr];
#else
                   int32_t bpc[nbpc];
                   int32_t bpr[nbpr];
#endif
                   constexpr T * bpc_beg() { return (std::addressof(bpc[0]));} 
                   constexpr T * bpr_beg() { return (std::addressof(bpr[0]));}
                   constexpr int32_t bpc_size() { return (nbpc);}
                   constexpr int32_t bpr_size() { return (nbpr);}
            };
        
            // latitude   values (deg), per building
            // longtitude values (deg), per building
            // Number of latitude   values (deg), per building
            // Number of longtitude values (deg), per building
            template<typename T,int32_t nblatd,int32_t nblond>
            struct BLatLondR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64) T blatd[nblatd];
                   __ATTR_ALIGN__(64) T blond[nblond];
#elif defined (__AVX__) || defined (__AVX2__)
                   __ATTR_ALIGN__(32) T blatd[nblatd];
                   __ATTR_ALIGN__(32) T blond[nblond];
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16) T blatd[nblatd];
                   __ATTR_ALIGN__(16) T blond[nblond];
#else
                     T blatd[nblatd];
                     T blond[nblond];
#endif
                   constexpr T * blatd_beg() { return (std::addressof(blatd[0]));} 
                   constexpr T * blond_beg() { return (std::addressof(blond[0]));}
                   constexpr int32_t blatd_size() { return (nblatd);}
                   constexpr int32_t blond_size() { return (nblond);}
            };
            
            
            // latitude   values (rad), per building
            // longtitude values (rad), per building
            // Number of latitude   values (rad), per building
            // Number of longtitude values (rad), per building
            template<typename T,int32_t nblatr,int32_t nblonr>
            struct BLatLonrR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64) T blatr[nblatr];
                   __ATTR_ALIGN__(64) T blonr[nblonr];
#elif defined (__AVX__) || defined (__AVX2__)
                   __ATTR_ALIGN__(32) T blatr[nblatr];
                   __ATTR_ALIGN__(32) T blonr[nblonr];
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16) T blatr[nblatr];
                   __ATTR_ALIGN__(16) T blonr[nblonr];
#else
                   T blatr[nblatr];
                   T blonr[nblonr];
#endif
                   constexpr T * blatr_beg() { return (std::addressof(blatr[0]));} 
                   constexpr T * blonr_beg() { return (std::addressof(blonr[0]));}
                   constexpr int32_t blatd_size() { return (nblatr);}
                   constexpr int32_t blond_size() { return (nblonr);}
            };
            
            
            // ellipsoidal (radar waveform irradiating field) cells for building column
            // Number of ellipsoidal (radar waveform irradiating field) cells for building column
            template<typename T,int32_t nellpb>
            struct EllpbR1x_t {
#if defined (__AVX512F__)                    
                  __ATTR_ALIGN__(64) T ellpb[nellpb];
#elif defined (__AVX__) || defined (__AVX2__)   
                  __ATTR_ALIGN__(32) T ellpb[nellpb];  
#elif defined (__SSE__)
                  __ATTR_ALIGN__(16) T ellpb[nellpb];   
#else
                  T ellpb[nellpb];   
#endif      
                  constexpr T * ellpb_beg() { return (std::addressof(ellpb[0]));}
                  constexpr int32_t ellpb_size() { return (nellpb);}
            };
            
            
            // ! Parametric equation x (acos(t)) values (building)
            // ! 1st dimension building column, 2nd dimension 'x' parameter values
            // ! Parametric equation y (b(sin(t)) values (building)
            // 1st dimension building column, 2nd dimension 'y' parameter values
            // number of building columns
            template<typename T,int32_t nbpc,
                     int32_t npxb,int32_t npyb>
            struct PxybR1x_t {
#if defined (__AVX512F__)                    
                   __ATTR_ALIGN__(64) T pxb[nbpc*npxb];
                   __ATTR_ALIGN__(64) T pyb[nbpc*npyb];
#elif defined (__AVX__) || defined (__AVX2__)
                   __ATTR_ALIGN__(32) T pxb[nbpc*npxb];
                   __ATTR_ALIGN__(32) T pyb[nbpc*npyb];
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16) T pxb[nbpc*npxb];
                   __ATTR_ALIGN__(16) T pyb[nbpc*npyb];
#else
                    T pxb[nbpc*npxb];
                    T pyb[nbpc*npyb];
#endif
                   constexpr T * pxb_beg() { return (std::addressof(pxb[0]));}
                   constexpr T * pyb_beg() { return (std::addressof(pyb[0]));}
                   constexpr int32_t pxb_size() { return (nbpc*npxb);}
                   constexpr int32_t pyb_size() { return (nbpc*npyb);}
            };
            
            
            // Length, width and an area of every street
             // number of streets
            template<typename T,int32_t nstr>
            struct SLWAR1x_t {
#if defined (__AVX512F__)                    
                    __ATTR_ALIGN__(64) T lstr[nstr];
                    __ATTR_ALIGN__(64) T wstr[nstr];
                    __ATTR_ALIGN__(64) T astr[nstr];
#elif defined (__AVX__) || defined (__AVX2__)
                    __ATTR_ALIGN__(32) T lstr[nstr];
                    __ATTR_ALIGN__(32) T wstr[nstr];
                    __ATTR_ALIGN__(32) T astr[nstr]; 
#elif defined (__SSE__)
                    __ATTR_ALIGN__(16) T lstr[nstr];
                    __ATTR_ALIGN__(16) T wstr[nstr];
                    __ATTR_ALIGN__(16) T astr[nstr];
#else
                     T lstr[nstr];
                     T wstr[nstr];
                     T astr[nstr];
#endif
                    constexpr T * lstr_beg() { return (std::addressof(lstr[0]));}
                    constexpr T * wstr_beg() { return (std::addressof(wstr[0]));}
                    constexpr T * astr_beg() { return (std::addressof(astr[0]));}
                    constexpr int32_t xstr_size() { return (nstr);}
            };
            
            
            //   ! Moisture of every street (2D array)
            //   ! 2nd dimension humidity values (per street), 1st dimension street numbers
            //    Percent of moist to dry area of evey street at each cell
            // number of streets
            template<typename T,int32_t nstr,
                     int32_t nmstr,int32_t npmstr>
            struct MStrR1x_t {
#if defined (__AVX512F__)                    
                    __ATTR_ALIGN__(64) T mstr[nstr*nmstr];
                    __ATTR_ALIGN__(64) T pmstr[nstr*npmstr];
#elif defined (__AVX__) || defined (__AVX2__)
                    __ATTR_ALIGN__(32) T mstr[nstr*nmstr];
                    __ATTR_ALIGN__(32) T pmstr[nstr*npmstr];
#elif defined (__SSE__)
                    __ATTR_ALIGN__(16) T mstr[nstr*nmstr];
                    __ATTR_ALIGN__(16) T pmstr[nstr*npmstr];
#else
                     T mstr[nstr*nmstr];
                     T pmstr[nstr*npmstr];
#endif                    
                    constexpr T * mstr_beg()  { return (std::addressof(mstr[0]));}
                    constexpr T * pmstr_beg() { return (std::addressof(pmstr[0]));}
                    constexpr int32_t mstr_size()  { return (nstr*nmstr);}
                    constexpr int32_t pmstr_size() { return (nstr*npmstr);}
            };
            
            
            // !Coverage of every street (like: '1' for snow,'2' for mud, '3' for clay, ...etc)
            // number of streets
            template<int32_t nstr>
            struct CStrIx_t {
#if defined (__AVX512F__)                    
                 __ATTR_ALIGN__(64) int32_t cstr[nstr];
#elif defined (__AVX__) || defined (__AVX2__)
                 __ATTR_ALIGN__(32) int32_t cstr[nstr];
#elif defined (__SSE__)
                 __ATTR_ALIGN__(16) int32_t cstr[nstr];
#else
                 int32_t cstr[nstr];
#endif
                 constexpr int32_t * cstr_beg() { return (std::addressof(cstr[0]));}
                 constexpr int32_t cstr_size()  { return (nstr);}
            };
            
            
            // Percent of covered to non-covered portion of every street at each irradiated cell
            // Average thickness of each layer (cover) of every street at each irradiated cell
            // Thickness of cover along street (number of values) at each irradiated cell
             // number of streets
            template<typename T,int32_t nstr,int32_t npcstr,
                              int32_t natstr,int32_t ntcstr>
            struct CDStrR1x_t {
#if defined (__AVX512F__)                           
                    __ATTR_ALIGN__(64) T pcstr[nstr*npcstr];
                    __ATTR_ALIGN__(64) T atstr[nstr*natstr];
                    __ATTR_ALIGN__(64) T tcstr[nstr*ntcstr];
#elif defined (__AVX__) || defined (__AVX2__)
                    __ATTR_ALIGN__(32) T pcstr[nstr*npcstr];
                    __ATTR_ALIGN__(32) T atstr[nstr*natstr];
                    __ATTR_ALIGN__(32) T tcstr[nstr*ntcstr]; 
#elif defined (__SSE__)
                    __ATTR_ALIGN__(16) T pcstr[nstr*npcstr];
                    __ATTR_ALIGN__(16) T atstr[nstr*natstr];
                    __ATTR_ALIGN__(16) T tcstr[nstr*ntcstr]; 
#else
                    T pcstr[nstr*npcstr];
                    T atstr[nstr*natstr];
                    T tcstr[nstr*ntcstr];
#endif                  
                    constexpr T * pcstr_beg() { return (std::addressof(pcstr[0]));}
                    constexpr T * atstr_beg() { return (std::addressof(atstr[0]));}
                    constexpr T * tcstr_beg() { return (std::addressof(tcstr[0]));}
                    constexpr int32_t pcstr_size() { return (nstr*npcstr);}
                    constexpr int32_t atstr_size() { return (nstr*natstr);}
                    constexpr int32_t tcstr_size() { return (nstr*ntcstr);}
            };
            
            
            // Mu values for 'clean' street interpolated along the street length at each irradiated cell
            // Eps for 'clean' street street length interpolated at each irradiated cell
             // number of streets
            template<typename T,int32_t nstr,
                     int32_t nmustr,int32_t nepstr> 
            struct MEStr1C1x_t {       
#if defined (__AVX512F__)             
                    __ATTR_ALIGN__(64) T murstr[nstr*nmustr];
                    __ATTR_ALIGN__(64) T muistr[nstr*nmustr];
                    __ATTR_ALIGN__(64) T eprstr[nstr*nepstr];
                    __ATTR_ALIGN__(64) T epistr[nstr*nepstr];
#elif defined (__AVX__) || (__AVX2__)
                    __ATTR_ALIGN__(32) T murstr[nstr*nmustr];
                    __ATTR_ALIGN__(32) T muistr[nstr*nmustr];
                    __ATTR_ALIGN__(32) T eprstr[nstr*nepstr];
                    __ATTR_ALIGN__(32) T epistr[nstr*nepstr];
#elif defined (__SSE__) 
                    __ATTR_ALIGN__(16) T murstr[nstr*nmustr];
                    __ATTR_ALIGN__(16) T muistr[nstr*nmustr];
                    __ATTR_ALIGN__(16) T eprstr[nstr*nepstr];
                    __ATTR_ALIGN__(16) T epistr[nstr*nepstr];
#else
                     T murstr[nstr*nmustr];
                     T muistr[nstr*nmustr];
                     T eprstr[nstr*nepstr];
                     T epistr[nstr*nepstr];
#endif
                    constexpr T * murstr_beg() { return (std::addressof(murstr[0]));}
                    constexpr T * muistr_beg() { return (std::addressof(muistr[0]));}
                    constexpr T * eprstr_beg() { return (std::addressof(eprstr[0]));}
                    constexpr T * epistr_beg() { return (std::addressof(epistr[0]));}
                    constexpr int32_t murstr_size() { return (nstr*nmustr);}
                    constexpr int32_t muistr_size() { return (nstr*nmustr);}
                    constexpr int32_t eprstr_size() { return (nstr*nepstr);}
                    constexpr int32_t epistr_size() { return (nstr*nepstr);}
            };
            
            
            // Mu for covered (i.e. by mud,snow,clay, ..etc) 
            // street interpolated along the street length at each irradiated cell
            // Eps for covered (i.e. by mud,snow,clay, ..etc) 
            // street  length interpolated at each irradiated cell
             // number of streets
            template<typename T,int32_t nstr,
                     int32_t nmustr,int32_t nepstr>
            struct MEStr2C1x_t {
#if defined (__AVX512F__)                     
                    __ATTR_ALIGN__(64) T murstr[nstr*nmustr];
                    __ATTR_ALIGN__(64) T muistr[nstr*nmustr];
                    __ATTR_ALIGN__(64) T eprstr[nstr*nepstr];
                    __ATTR_ALIGN__(64) T epistr[nstr*nepstr];  
#elif defined (__AVX__) || defined (__AVX2__)
                    __ATTR_ALIGN__(32) T murstr[nstr*nmustr];
                    __ATTR_ALIGN__(32) T muistr[nstr*nmustr];
                    __ATTR_ALIGN__(32) T eprstr[nstr*nepstr];
                    __ATTR_ALIGN__(32) T epistr[nstr*nepstr];  
#elif defined (__SSE__)
                    __ATTR_ALIGN__(16) T murstr[nstr*nmustr];
                    __ATTR_ALIGN__(16) T muistr[nstr*nmustr];
                    __ATTR_ALIGN__(16) T eprstr[nstr*nepstr];
                    __ATTR_ALIGN__(16) T epistr[nstr*nepstr];  
#else
                     T murstr[nstr*nmustr];
                     T muistr[nstr*nmustr];
                     T eprstr[nstr*nepstr];
                     T epistr[nstr*nepstr];  
#endif
                    constexpr T * murstr_beg() { return (std::addressof(murstr[0]));}
                    constexpr T * muistr_beg() { return (std::addressof(muistr[0]));}
                    constexpr T * eprstr_beg() { return (std::addressof(eprstr[0]));}
                    constexpr T * epistr_beg() { return (std::addressof(epistr[0]));}
                    constexpr int32_t murstr_size() { return (nstr*nmustr);}
                    constexpr int32_t muistr_size() { return (nstr*nmustr);}
                    constexpr int32_t eprstr_size() { return (nstr*nepstr);}
                    constexpr int32_t epistr_size() { return (nstr*nepstr);}
            };
            
            
            // Street curvature parametric equation u-parameter
            // Street curvature parametric equation v-parameter
            // number of streets
            template<typename T,int32_t nstr,
                     int32_t nustr,int32_t nvstr>
            struct SCrvR1x_t {
#if defined (__AVX512F__)                     
                    __ATTR_ALIGN__(64) T ustr[nstr*nustr];
                    __ATTR_ALIGN__(64) T vstr[nstr*nvstr];
#elif defined (__AVX__) || defined (__AVX2__)
                    __ATTR_ALIGN__(32) T ustr[nstr*nustr];
                    __ATTR_ALIGN__(32) T vstr[nstr*nvstr];
#elif defined (__SSE__)
                    __ATTR_ALIGN__(16) T ustr[nstr*nustr];
                    __ATTR_ALIGN__(16) T vstr[nstr*nvstr];
#else
                     T ustr[nstr*nustr];
                     T vstr[nstr*nvstr];
#endif
                 
            };
            
            
            //  Street surface normal vectors x-components 
            //  along the street length at each irradiated cell
            //  Street surface normal vectors y-components 
            //  along the street length at each irradiated cell
            //  Street surface normal vectors z-components 
            //  along the street length at each irradiated cell
             // number of streets
            template<typename T,int32_t nstr,
                     int32_t nx,int32_t ny,int32_t,nz>
            struct SNrmR1x_t {
                    
                    __ATTR_ALIGN__(64) T nvx[nstr*nx];
                    __ATTR_ALIGN__(64) T nvy[nstr*ny];
                    __ATTR_ALIGN__(64) T nvz[nstr*nz]
                    
            };
            
            
            // latitude   values (deg), per street length (at irradiance point)
            // longtitude values (deg), per street length (at irradiance point)
           // number of streets
            template<typename T,int32_t nstr,
                     int32_t nlon,int32_t nlat>
            struct SIRCDR1x_t {
                     
                    __ATTR_ALIGN__(64) T irlon[nstr*nlon];
                    __ATTR_ALIGN__(64) T irlat[nstr*nlat];
                   
            } 
            
            
           // latitude   values (rad), per street length (at irradiance point)
           // longtitude values (rad), per street length (at irradiance point)
           // number of streets
           template<typename T,int32_t nstr,
                    int32_t nlon,int32_t nlat>
           struct SIRCRR1x_t {
                    
                    __ATTR_ALIGN__(64) T irlon[nstr*nlon];
                    __ATTR_ALIGN__(64) T irlat[nstr*nlat];
           };
           
           
           // latitude   values (deg), of the building area (at irradiance point)
           // longtitude values (deg), of the building area (at irradiance point)
            // number of buildings
            template<typename T,int32_t nbld,
                     int32_t nlon,int32_t nlat>
            struct BIRCDR1x_t {
                     
                    __ATTR_ALIGN__(64) T irlon[nbld*nlon]
                    __ATTR_ALIGN__(64) T irlat[nbld*nlat];
                    
            };
            
            
           // latitude   values (rad), of the building area (at irradiance point)
           // longtitude values (rad), of the building area (at irradiance point)
            // number of buildings
           template<typename T,int32_t nbld,
                    int32_t nlon,int32_t nlat>
           struct BIRCRR1x_t {
                    
                    __ATTR_ALIGN__(64) T irlon[nbld*nlon]
                    __ATTR_ALIGN__(64) T irlat[nbld*nlat];
                 
           }; 
           
           
           // Urban area height map (at single building resolution)
           template<typename T,int32_t nx,int32_t ny>
           struct UHMapR1x_t {
                  
                  __ATTR_ALIGN__(64) T hmap[nx*ny];
                 
           };
           
           
           // Urban area height map (at single building resolution) -- 1st derivative
           template<typename T,int32_t nx,int32_t ny>
           struct UHDxDyR1x_t {
                  
                  __ATTR_ALIGN__(64) T hdxdy[nx*ny];
           };
           
           
           // Urban area height map (at single building resolution) -- gradient x-component
           // Urban area height map (at single building resolution) -- gradient y-component
           template<typename T,int32_t nx,int32_t ny>
           struct UHGradR1x_t {
                  
                  __ATTR_ALIGN__(64) T uhgx[nx];
                  __ATTR_ALIGN__(64) T uhgy[ny];
           };
           
           
           // Smoothing and approximating curve for linearly-piecewise height function (x-coordinate)
           // Smoothing and approximating curve for linearly-piecewise height function (y-coordinate)
           template<typename T,int32_t nx,int32_t ny>
           struct XYSMBHR1x_t {
                   
                  __ATTR_ALIGN__(64) T xsmbh[nx];
                  __ATTR_ALIGN__(64) T ysmbh[ny];
           };
           
           
           // Empty space in-between of buildings (per single column) x number columns
           template<int32_t ncols,int32_t nval>
           struct ESBBI1x_t {
                  
                 __ATTR_ALIGN__(64) int32_t esbb[ncols*nval];
           };
           
           
           // An area values of in-between buildings empty spaces (per single column) x number columns
           template<typename T,int32_t ncols,int32_t nval> 
           struct AESBBR1x_t {
                  
                  __ATTR_ALIGN__(64) T aesbb[ncols*nval];
           };
           
           
           // An area values of each building (per single building column) x number columns
           template<typename T,int32_t ncols,int32_t nval>
           struct ABCR1x_t {
                  
                   __ATTR_ALIGN__(64) T  abc[ncols*nvals];
           };
           
           
           // Number of south-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct SWPCI1x_t {
                  
                   __ATTR_ALIGN__(64) int32_t swpc[ncols];
           };
           
           
           // Number of east-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct EWPCI1x_t {
                  
                 __ATTR_ALIGN__(64) int32_t ewpc[ncols];
           };
           
           
           // Number of west-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct WWPCI1x_t {
                   
                __ATTR_ALIGN__(64) int32_t  wwpc[ncols];
           };
           
           
           // Number of north-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct NWPCI1x_t {
                  
                __ATTR_ALIGN__(64) int32_t  nwpc[ncols];
           };
           
           
           // Number of building roofs per each column x number of columns
           template<int32_t ncols,int32_t nval>
           struct BRPCI1x_t {
                  
                 __ATTR_ALIGN__(64) int32_t brpc; 
           };
           
           
           //  An area of every building [flat] roof (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct BRAPCR1x_t {
                  
                  __ATTR_ALIGN__(64) T brapc[ncols*nval];
           };
           
           
           // Number of angled roof -- south facing roof wall (per each column) x number of columns
           template<int32_t ncols>
           struct SRWCI1x_t {
                  
                 __ATTR_ALIGN__(64) int32_t srwc[ncols];
           };
           
           
           // Number of angled roof -- east facing roof wall (per each column) x number of columns
           template<int32_t ncols>
           struct ERWCI1x_t {
                  
                  __ATTR_ALIGN__(64) int32_t  erwc[ncols];
           };
           
           
           // Number of angled roof -- west facing roof wall (per each column)  x number of columns
           template<int32_t ncols>
           struct WRWCI1x_t {
                  
                __ATTR_ALIGN__(64) int32_t  wrwc[ncols];
           };
           
           
           // Number angled roof -- north facing roof wall (per each column) x number of columns
           template<int32_t ncols>
           struct NRWCI1x_t {
                  
                __ATTR_ALIGN__(64) int32_t nrwc[ncols];
           };
           
           
           // An angled roof inclination (deg) -- south facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IDSRWR1x_t {
                   
                 __ATTR_ALIGN__(64) T idsrw[ncols*nval];
           };
           
           
           // An angled roof inclination (deg) -- east facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IDERWR1x_t {
                   
                  __ATTR_ALIGN__(64) T iderw[ncols*nval];
           };
           
           
           // An angled roof inclination (deg) -- west facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IDWRWR1x_t {
                
                  __ATTR_ALIGN__(64) T idwrw[ncols*nval];
           };
           
           
           // An angled roof inclination (rad) -- south facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IDNRWR1x_t {
                   
                  __ATTR_ALIGN__(64) T idnrw[ncols*nval]; 
           };
           
           
           //  An angled roof inclination (rad) -- south facing roof wall 
           //  (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IRSRWR1x_t {
                   
                  __ATTR_ALIGN__(64) T irsrw[ncols*nval]; 
           };
           
           
           //  An angled roof inclination (rad) -- east facing roof wall 
           //  (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IRERWR1x_t {
                    
                   __ATTR_ALIGN__(64) T irerw[ncols*nval];
           };
           
           
           //  An angled roof inclination (rad) -- north facing roof wall 
           //  (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IRNRWR1x_t {
                    
                   __ATTR_ALIGN__(64) T irnrw[ncols*nval]; 
           };
           
           
           // An angled roof inclination surface area -- south facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct ISRAR1x_t {
                   
                  __ATTR_ALIGN__(64) T isra[ncols*nval];  
           };
           
           
           // An angled roof inclination surface area -- west facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IWRAR1x_t {
                    
                   __ATTR_ALIGN__(64) T iwra[ncols*nval];  
           };
           
           
           // An angled roof inclination surface area -- east facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IERAR1x_t {
                 
                  __ATTR_ALIGN__(64) T iera[ncols*nval];   
           };
           
           
           // An angled roof inclination surface area -- north facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct INRAR1x_t {
                   
                  __ATTR_ALIGN__(64) T inra[ncols*nval];   
           };
           
           
           // South wall upper-facing edge inclination (rad) -- 
           // (per each column)  x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct SWUER1x_t {
                   
                 __ATTR_ALIGN__(64) T swue[ncols*nval];    
          };
          
          
          // East wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct EWUER1x_t {
                   
                 __ATTR_ALIGN__(64) T ewue[ncols*nval];    
          };
          
          
          // West wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct WWUER1x_t {
                  
                 __ATTR_ALIGN__(64) T wwue[ncols*nval];     
          };
          
          
          // North wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct NWUER1x_t {
                   
                 __ATTR_ALIGN__(64) T nwue[ncols*nval];       
          };
          
          
          // Shared right edge between the south wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct SEWER1x_t {
                   
                 __ATTR_ALIGN__(64) T sewe[ncols*nval];     
          };
          
          // Shared left edge between the south wall and west wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct SWWER1x_t {
                   
                 __ATTR_ALIGN__(64) T swwe[ncols*nval];   
          };
          
          
           // Shared right edge between the north wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct NWEER1x_t {
                    
                  __ATTR_ALIGN__(64) T nwee[ncols*nval];   
           };
           
           
           // Shared right edge between the north wall and west wall inclination (rad) 
           // ! -- (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct NWWER1x_t {
                   
                  __ATTR_ALIGN__(64) T nwee[ncols*nval];   
           };
           
           
           // Simple cell-based mesh
           template<typename T,int32_t nL>
           struct CellMeshR1x_t {
                  
                  // Coordinates (x,y,z) of the center
                  // of Lth cell.
                  __ATTR_ALIGN__(64) T cx[nL];
                  __ATTR_ALIGN__(64) T cy[nL];
                  __ATTR_ALIGN__(64) T cz[nL];
                  __ATTR_ALIGN__(64) T dv[nL];
                  // (X,Y,Z) dimensions of the Ith
                  // rectangular volume cell (this is needed for
                  // the numerical integration)
                  __ATTR_ALIGN__(64) T dx[nL];
                  __ATTR_ALIGN__(64) T dy[nL];
                  __ATTR_ALIGN__(64) T dz[nL];
                  // Number of divisions along the x,y,z
                  int32_t ndiv[3];
                  // Compute numerical integration.
                  bool nint;
                 
                  
           };
           
           
           // South walls surface area (for every building, per column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct SWSAR1x_t {
                   
                  __ATTR_ALIGN__(64) T swsa[ncols*nval];
           };
           
           
           // East walls surface area (for every building, per column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct EWSAR1x_t {
                   
                  __ATTR_ALIGN__(64) T ewsa[ncols*nval];
           };
           
           
           // West walls surface area (for every building, per column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct WWSAR1x_t {
                   
                  __ATTR_ALIGN__(64) T wwsa[ncols*nval];
           };
           
           
           // North walls surface area (for every building, per column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct NWSAR1x_t {
                    
                  __ATTR_ALIGN__(64) T wwsa[ncols*nval]; 
                  T * begin() { return (std::__adressof(wwsa[0]));}
                  constexpr int32_t size() { return (ncols*nval)};  
                 
           };
           
           
           // South walls moist/non moist logical (per column) x number of columns
           struct MNMSWB1x_t {
                    
                   int32_t ncols;
                   int32_t nval;
                   bool * __restrict mnmsw; 
           };
           
           
            // East walls moist/non moist logical (per column) x number of columns
           struct MNMEWB1x_t {
                    
                   int32_t ncols;
                   int32_t nval;
                   bool * __restrict mnmew; 
           };
           
           
             // West walls moist/non moist logical (per column) x number of columns
           struct MNMWWB1x_t {
                    
                   int32_t ncols;
                   int32_t nval;
                   bool * __restrict mnmww; 
           };
           
           
             // North walls moist/non moist logical (per column) x number of columns
           struct MNMNWB1x_t {
                    
                   int32_t ncols;
                   int32_t nval;
                   bool * __restrict mnmnw; 
           };
           
           
           // ! The values describing the ratio (percentage) of south wall 
                            // ! moisture to dryness (per each column) x number of columns
           template<typename T>
           struct MDSWRR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mdswr;
           };
           
           
           // ! The values describing the ratio (percentage) of east wall 
                            // ! moisture to dryness (per each column) x number of columns
           template<typename T>
           struct MDEWRR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mdewr;
           };
           
           
           //  The values describing the ratio (percentage) of west wall 
                             //! moisture to dryness (per each column) x number of columns
           template<typename T>
           struct MDWWRR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mdwwr;
           };  
           
           
           //  The values describing the ratio (percentage) of north wall 
                             //! moisture to dryness (per each column) x number of columns              
           template<typename T>
           struct MDNWRR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mdnwr;
           };   
           
           
           // The logical values of flat roof moistness (being either moist or dry) 
                              // ! (per column) x number of columns                
          struct MDRB1x_t {
                 
                 int32_t ncols;
                 int32_t nval;
                 bool * __restrict mdr;
          }; 
          
          
          // The values describing the ratio (percentage) of flat roof moisture to dryness 
                              // ! (per each column) x number of columns 
          template<typename T>
          struct MDRRR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mdrr; 
          };
          
          
          // The values describing the surface of moist part of the flat roof 
                               // ! (per each column) x number of columns
          template<typename T>
          struct MPFRR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mpfr; 
          };     
          
          
          // The values describing the surface of dry part of the flat roof 
                               // ! (per each column) x number of columns  
          template<typename T>
          struct DPFRR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   dpfr; 
          };   
          
          
          //  The values describing the surface of moist part of the south wall 
                                // ! (per each column) x number of columns   
          template<typename T>
          struct MPSWR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mpsw; 
          };  
          
          
          //  The values describing the surface of dry part of the south wall 
                               //  ! (per each column) x number of columns  
          template<typename T>
          struct DPSWR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   dpsw; 
          };   
          
          
          //  The values describing the surface of moist part of the east wall 
                                // ! (per each column) x number of columns
          template<typename T>
          struct MPEWR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mpew; 
          };   
          
         // The values describing the surface of dry part of the east wall 
                                // ! (per each column) x number of columns 
          template<typename T>
          struct DPEWR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   dpew; 
          };  
          
          
         // The values describing the surface of moist part of the west wall 
                                 //! (per each column) x number of columns
          template<typename T>
          struct MPWWR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mpww; 
          }; 
          
          
        //  The values describing the surface of dry part of the west wall 
                                 //! (per each column) x number of columns 
          template<typename T>
          struct DPWWR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   dpww; 
          };  
          
          
        // The values describing the surface of moist part of the north wall 
                                 //! (per each column) x number of columns
         template<typename T>
         struct MPNWR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mpnw; 
          }; 
          
          
         //  The values describing the surface of dry part of the north wall 
                                // ! (per each column) x number of columns
         template<typename T>
         struct DPNWR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   dpnw; 
          }; 
          
          
          // The values describing the surface of moist part of the angled south roof wall
                               // ! (per each column) x number of columns
          template<typename T>
          struct MPSARR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mpsar;
          };
          
          
          // The values describing the surface of dry part of the angled south roof wall
                               // ! (per each column) x number of columns
          template<typename T>
          struct DPSARR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   dpsar;
          };
          
          
          // The values describing the surface of moist part of the angled east roof wall
                               // ! (per each column) x number of columns 
          template<typename T>
          struct MPEARR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mpear;
          };
          
          
         // The values describing the surface of dry part of the angled east roof wall
                               // ! (per each column) x number of columns  
         template<typename T>
         struct DPEARR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   dpear;
          }; 
          
          
          // The values describing the surface of moist part of the angled west roof wall
                               // ! (per each column) x number of columns  
          template<typename T>
          struct MPWARR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mpwar;
          }; 
          
          
           // The values describing the surface of dry part of the angled west roof wall
                               // ! (per each column) x number of columns  
          template<typename T>
          struct DPWARR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   dpwar;
          }; 
          
          
           // The values describing the surface of moist part of the angled north roof wall
                               // ! (per each column) x number of columns  
          template<typename T>
          struct MPNARR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   mpnar;
          }; 
          
          
          // The values describing the surface of dry part of the angled north roof wall
                               // ! (per each column) x number of columns  
          template<typename T>
          struct DPNARR1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   dpnar;
          }; 
          
          
         // The values describing the complex permittivity of south walls
                               // ! (per each column) x number of columns  
         template<typename T>
         struct CESWC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cesw;
         };
         
         
         // The values describing the complex permeabillity of south walls
                               // ! (per each column) x number of columns  
         template<typename T>
         struct CMSWC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cmsw;
         };
         
         
         // The values describing the complex permittivity of west walls
                               // ! (per each column) x number of columns 
         template<typename T>
         struct CEWWC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   ceww;
         };
         
         
          // The values describing the complex permeability of west walls
                               // ! (per each column) x number of columns 
         template<typename T>
         struct CMWWC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cmww;
         };
         
         
          // The values describing the complex permittivity of east walls
                               // ! (per each column) x number of columns 
         template<typename T>
         struct CEEWC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   ceew;
         };
         
         
          // The values describing the complex permeability of east walls
                               // ! (per each column) x number of columns 
         template<typename T>
         struct CMEWC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cmew;
         };
         
         
          // The values describing the complex permittivity of north walls
                               // ! (per each column) x number of columns 
         template<typename T>
         struct CENWC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cenw;
         };
         
         
          // The values describing the complex permeability of north walls
                               // ! (per each column) x number of columns 
         template<typename T>
         struct CMNWC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cmnw;
         };
         
         
          // The values describing the complex permittivity of south angled roof
                               // ! (per each column) x number of columns  
         template<typename T>
         struct CESARC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cesar;
         };
         
         
         // The values describing the complex permeabillity of south angled roof
                               // ! (per each column) x number of columns  
         template<typename T>
         struct CMSARC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cmsar;
         };
         
         
         // The values describing the complex permittivity of east angled roof
                               // ! (per each column) x number of columns  
         template<typename T>
         struct CEEARC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   ceear;
         };
         
         
         // The values describing the complex permeabillity of east angled roof
                               // ! (per each column) x number of columns  
         template<typename T>
         struct CMEARC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cmear;
         };
         
         
         // The values describing the complex permittivity of west angled roof
                               // ! (per each column) x number of columns  
         template<typename T>
         struct CEWARC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cewar;
         };
         
         
         // The values describing the complex permeabillity of west angled roof
                               // ! (per each column) x number of columns  
         template<typename T>
         struct CMWARC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cmwar;
         };
         
         
          // The values describing the complex permittivity of north angled roof
                               // ! (per each column) x number of columns  
         template<typename T>
         struct CENARC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cenar;
         };
         
         
         // The values describing the complex permeabillity of north angled roof
                               // ! (per each column) x number of columns  
         template<typename T>
         struct CMNARC1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC2D<T>_t   cmnar;
         };
         
         
         // The components of south walls normal vector
                                    // ! (per each column) x number of columns  
         template<typename T>
         struct NVSWR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D<T>_t   nvx;
                  DC1D<T>_t   nvy;
                  DC1D<T>_t   nvz;
         };
         
         
        // The components of east walls normal vector
                                    // ! (per each column) x number of columns  
         template<typename T>
         struct NVEWR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D<T>_t   nvx;
                  DC1D<T>_t   nvy;
                  DC1D<T>_t   nvz;
         };
         
         
        // The components of west walls normal vector
                                    // ! (per each column) x number of columns  
         template<typename T>
         struct NVWWR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D<T>_t   nvx;
                  DC1D<T>_t   nvy;
                  DC1D<T>_t   nvz;
         };
         
         
        // The components of north walls normal vector
                                    // ! (per each column) x number of columns  
         template<typename T>
         struct NVNWR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D<T>_t   nvx;
                  DC1D<T>_t   nvy;
                  DC1D<T>_t   nvz;
         };
         
         
         // The components of each building normal vector
                                    // ! (per each column) x number of columns 
         template<typename T>
         struct NVBR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D<T>_t   nvx;
                  DC1D<T>_t   nvy;
                  DC1D<T>_t   nvz;  
         };
         
         
         // The components of each building flat roof normal vector
                                    // ! (per each column) x number of columns 
         template<typename T>
         struct NVFRR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D<T>_t   nvx;
                  DC1D<T>_t   nvy;
                  DC1D<T>_t   nvz;   
         };
         
         
          // The components of south angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<typename T>
         struct NVSARR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D<T>_t   nvx;
                  DC1D<T>_t   nvy;
                  DC1D<T>_t   nvz;
         };
         
         
        // The components of east angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<typename T>
         struct NVEARR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D<T>_t   nvx;
                  DC1D<T>_t   nvy;
                  DC1D<T>_t   nvz;
         };
         
         
        // The components of west angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<typename T>
         struct NVWARR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D<T>_t   nvx;
                  DC1D<T>_t   nvy;
                  DC1D<T>_t   nvz;
         };
         
         
        // The components of north angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<typename T>
         struct NVNARR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D<T>_t   nvx;
                  DC1D<T>_t   nvy;
                  DC1D<T>_t   nvz;
         };
         
         
         // The values of each south wall height and width
                        // ! (per each column) x number of columns  
         template<typename T>
         struct HWSWR1x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D<T>_t   hsw;
                DC1D<T>_t   wsw;
         };
         
         
          // The values of each east wall height and width
                        // ! (per each column) x number of columns  
         template<typename T>
         struct HWEWR1x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D<T>_t   hew;
                DC1D<T>_t   wew;
         };
         
         
         // The values of each west wall height and width
                        // ! (per each column) x number of columns  
         template<typename T>
         struct HWWWR1x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D<T>_t   hww;
                DC1D<T>_t   www;
         };
         
         
          // The values of each north wall height and width
                        // ! (per each column) x number of columns  
         template<typename T>
         struct HWNWR1x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D<T>_t   hnw;
                DC1D<T>_t   wnw;
         };
         
         
         // The values of each flat roof height and width
                        // ! (per each column) x number of columns  
         template<typename T>
         struct HWFRR1x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D<T>_t   hfr;
                DC1D<T>_t   wfr;
         };
         
         
        // The values of each south non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         template<typename T>
         struct HWSNFRR1x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D<T>_t   hsnfr;
                DC1D<T>_t   wsnfr;
         };
         
         
        // The values of each east non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         template<typename T>
         struct HWENFRR1x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D<T>_t   henfr;
                DC1D<T>_t   wenfr;
         };
         
         
        // The values of each west non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         template<typename T>
         struct HWWNFRR1x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D<T>_t   hwnfr;
                DC1D<T>_t   wwnfr;
         };
         
         
        // The values of each north non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         template<typename T>
         struct HWNNFRR1x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D<T>_t   hnnfr;
                DC1D<T>_t   wnnfr;
         };
         
         
         // Any kind of metallic structure fixed on the roof, e.g., an wire antenna, cylindrical object (ventillation)
         // or parabollic antenna, or similar type (per each column) x number of columns.
         struct BRMSB1x_t {
                
                int32_t ncols;
                int32_t nval;
                bool * __restrict want; // a wire antennae
                bool * __restrict pant; // a parabollic antennae
                bool * __restrict yant; // yagi type antennae
                bool * __restrict lpda; // log-periodic dipole array
                bool * __restrict cant; // cell phone sector bars antennae
                bool * __restrict cylo; // any kind of cylindrical (ventillation) object
         };
         
         
         // The number of ventillation objects per single building roof
         // for every building column.
         struct NVOI1x_t {
                
                int32_t ncols;
                int32_t nval;
                int32_t * __restrict nvo;
         };
         
         
         // The number of wire antennae per single building roof
         // for every building column.
         struct NWAI1x_t {
                
                int32_t ncols;
                int32_t nval;
                int32_t * __restrict nwa;
         };
         
         
          // The number of yagi-antennae per single building roof
         // for every building column.
         struct NYAI1x_t {
                
                int32_t ncols;
                int32_t nval;
                int32_t * __restrict nya;
         };
         
         
          // The number of log-periodic dipole antennae per single building roof
         // for every building column.
         struct NLPDAI1x_t {
                
                int32_t ncols;
                int32_t nval;
                int32_t * __restrict nlpda;
         };
         
         
          // The number of parabollic antennae per single building roof
         // for every building column.
         struct NPAI1x_t {
                
                int32_t ncols;
                int32_t nval;
                int32_t * __restrict npa;
         };
         
         
          // The number of cell-phone antennae per single building roof
         // for every building column.
         struct NCPAI1x_t {
                
                int32_t ncols;
                int32_t nval;
                int32_t * __restrict ncpa;
         };
         
         
         // The values of RCS for the flat roof of
         // building column.
         template<typename T>
         struct RCSFRR1x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t   rcsfr;
         };
         
         
          // The values of RCS for the south wall of
         // of every building in the building column.
         template<typename T>
         struct RCSSWR1x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t   rcssw;
         };
         
         
          
          // The values of RCS for the east wall of
         // of every building in the building column.
         template<typename T>
         struct RCSEWR1x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t   rcsew;
         };
         
         
          
          // The values of RCS for the west wall of
         // of every building in the building column.
         template<typename T>
         struct RCSWWR1x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t   rcsww;
         };
         
         
          
          // The values of RCS for the north wall of
         // of every building in the building column.
         template<typename T>
         struct RCSNWR1x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t   rcsnw;
         };
         
         
         // The values of RCS for the south angled roof of
         // of every building in the building column.
         template<typename T>
         struct RCSSARR1x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t   rcssar;
         };
         
         
         // The values of RCS for the south east roof of
         // of every building in the building column.
         template<typename T>
         struct RCSEARR1x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t   rcsear;
         };
         
         
         // The values of RCS for the west angled roof of
         // of every building in the building column.
         template<typename T>
         struct RCSWARR1x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t   rcswar;
         };
         
         
         // The values of RCS for the north angled roof of
         // of every building in the building column.
         template<typename T>
         struct RCSNARR1x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t   rcsnar;
         };
         
         
         // The values of whole building surface area
         // of every building in building column
         template<typename T>
         struct WBSAR1x_t {
                 
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t wbsa;
         };
         
         
         // The values of whole building internal volume
         // of every building in building column
         template<typename T>
         struct WBIVR1x_t {
                 
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D<T>_t wbiv;
         };
         
         
         
         
           
     }// radiolocation



}






























#endif /*__GMS_URBAN_ADT_STAT_H__*/
