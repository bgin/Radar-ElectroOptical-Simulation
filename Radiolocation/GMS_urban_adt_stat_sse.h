

#ifndef __GMS_URBAN_ADT_STAT_SSE_H__
#define __GMS_URBAN_ADT_STAT_SSE_H__ 140120240839


namespace file_info {

     const unsigned int GMS_URBAN_ADT_STAT_SSE_MAJOR = 1;
     const unsigned int GMS_URBAN_ADT_STAT_SSE_MINOR = 0;
     const unsigned int GMS_URBAN_ADT_STAT_SSE_MICRO = 0;
     const unsigned int GMS_URBAN_ADT_STAT_SSE_FULLVER =
       1000U*GMS_URBAN_ADT_STAT_SSE_MAJOR+100U*GMS_URBAN_ADT_STAT_SSE_MINOR+
       10U*GMS_URBAN_ADT_STAT_SSE_MICRO;
     const char * const GMS_URBAN_ADT_STAT_SSE_CREATION_DATE = "14-01-2024 08:39AM +00200 (SUN 14 JAN 2024 GMT+2)";
     const char * const GMS_URBAN_ADT_STAT_SSE_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_URBAN_ADT_STAT_SSE_SYNOPSIS      = "Abstract data types (static allocation of SSE type) representing built-up area for diffraction and scattering calculations."

}

#include <immintrin.h>
#include <cstdint>
#include <memory>
#include "GMS_config"


namespace gms {

       namespace radiolocation {
    
    
            // Building units per column
            // Building units per row
            template<int32_t nbpc,int32_t nbpr>
            struct SSEBRCcount_t {

                 __ATTR_ALIGN__(16) __m128i bpc[nbpc];
                 __ATTR_ALIGN__(16) __m128i bpr[nbpr];
                   constexpr static int32_t NBPC = nbpc;
                   constexpr static int32_t NBPR = nbpr;
                   constexpr __m128i * __restrict__ bpc_beg() { return (std::addressof(bpc[0]));} 
                   constexpr __m128i * __restrict__ bpr_beg() { return (std::addressof(bpr[0]));}
                   int32_t * __restrict__ bpc_iptr() { return ((int32_t* __restrict__)&bpc[0]);}
                   int32_t * __restrict__ bpr_iptr() { return ((int32_t* __restrict__)&bpr[0]);}
                   constexpr int32_t bpc_size() { return (NBPC);}
                   constexpr int32_t bpr_size() { return (NBPR);}
                   
            };
        
            // latitude   values (deg), per building
            // longtitude values (deg), per building
            // Number of latitude   values (deg), per building
            // Number of longtitude values (deg), per building
            template<int32_t nblatd,int32_t nblond>
            struct SSEBLatLondR1x_t {

                 __ATTR_ALIGN__(16)  __m128 blatd[nblatd];
               
                 __ATTR_ALIGN__(16)  __m128 blond[nblond];
                   constexpr static int32_t NBLATD = nblatd;
                   constexpr static int32_t NBLOND = nblond;
                   constexpr __m128 * __restrict__ blatd_beg() { return (std::addressof(blatd[0]));} 
                   constexpr __m128 * __restrict__ blond_beg() { return (std::addressof(blond[0]));}
                   float * __restrict__ blatd_fptr() { return ((float* __restrict__)&blatd[0]);}
                   float * __restrict__ blond_fptr() { return ((float* __restrict__)&blond[0]);}
                   constexpr int32_t blatd_size() { return (NBLATD);}
                   constexpr int32_t blond_size() { return (NBLOND);}
                   
            };
            
            
            // latitude   values (rad), per building
            // longtitude values (rad), per building
            // Number of latitude   values (rad), per building
            // Number of longtitude values (rad), per building
            template<int32_t nblatr,int32_t nblonr>
            struct SSEBLatLonrR1x_t {

                 __ATTR_ALIGN__(16)  __m128 blatr[nblatr];
                
                 __ATTR_ALIGN__(16)  __m128 blonr[nblonr];
                   constexpr static int32_t NBLATR = nblatr;
                   constexpr static int32_t NBLONR = nblonr;
                   constexpr __m128 * __restrict blatr_beg() { return (std::addressof(blatr[0]));} 
                   constexpr __m128 * __restrict blonr_beg() { return (std::addressof(blonr[0]));}
                   float * __restrict blatr_fptr() { return ((float* __restrict)&blatr[0]);}
                   float * __restrict blonr_fptr() { return ((float* __restrict)&blonr[0]);}
                   constexpr int32_t blatd_size() { return (NBLATR);}
                   constexpr int32_t blond_size() { return (NBLONR);}
                   
            };
            
            
            // ellipsoidal (radar waveform irradiating field) cells for building column
            // Number of ellipsoidal (radar waveform irradiating field) cells for building column
            template<int32_t nellpb>
            struct SSEEllpbR1x_t {

                  __ATTR_ALIGN__(16)  __m128 ellpb[nellpb];   
                  constexpr static int32_t NELLPB = nellpb;
                  constexpr __m128 * __restrict ellpb_beg() { return (std::addressof(ellpb[0]));}
                  float * __restrict ellpb_fptr() { ((float* __restrict)&ellpb[0];)}
                  constexpr int32_t ellpb_size() { return (NELLPB);}
            };
            
            
            // ! Parametric equation x (acos(t)) values (building)
            // ! 1st dimension building column, 2nd dimension 'x' parameter values
            // ! Parametric equation y (b(sin(t)) values (building)
            // 1st dimension building column, 2nd dimension 'y' parameter values
            // number of building columns
            template<int32_t nbpc,
                     int32_t npxb,int32_t npyb>
            struct SSEPxybR1x_t {

                   __ATTR_ALIGN__(16)  __m128 pxb[nbpc*npxb];
                   __ATTR_ALIGN__(16)  __m128 pyb[nbpc*npyb];
                   constexpr static int32_t NBPC = nbpc;
                   constexpr static int32_t NPXB = npxb;
                   constexpr static int32_t NPYB = npyb;
                   constexpr __m128 * __restrict pxb_beg() { return (std::addressof(pxb[0]));}
                   constexpr __m128 * __restrict pyb_beg() { return (std::addressof(pyb[0]));}
                   float * __restrict pxb_fptr() { return ((float*)&pxb[0]);}
                   float * __restrict pyb_fptr() { return ((float*)&pyb[0]);}
                   constexpr int32_t pxb_size() { return (NBPC*NPXB);}
                   constexpr int32_t pyb_size() { return (NBPC*NPYB);}
            };
            
            
            // Length, width and an area of every street
             // number of streets
            template<int32_t nstr>
            struct SSESLWAR1x_t {

                    __ATTR_ALIGN__(16)  __m128 lstr[nstr];
                    __ATTR_ALIGN__(16)  __m128 wstr[nstr];
                    __ATTR_ALIGN__(16)  __m128 astr[nstr];
                    constexpr static int32_t NSTR = nstr;
                    constexpr __m128 * __restrict lstr_beg() { return (std::addressof(lstr[0]));}
                    constexpr __m128 * __restrict wstr_beg() { return (std::addressof(wstr[0]));}
                    constexpr __m128 * __restrict astr_beg() { return (std::addressof(astr[0]));}
                    float * __restrict lstr_fptr() { return ((float* __restrict)&lstr[0]);}
                    float * __restrict wstr_fptr() { return ((float* __restrict)&wstr[0]);}
                    float * __restrict astr_fptr() { return ((float* __restrict)&astr[0]);}
                    constexpr int32_t xstr_size() { return (NSTR);}
            };
            
            
            //   ! Moisture of every street (2D array)
            //   ! 2nd dimension humidity values (per street), 1st dimension street numbers
            //    Percent of moist to dry area of evey street at each cell
            // number of streets
            template<int32_t nstr,
                     int32_t nmstr,int32_t npmstr>
            struct SSEMStrR1x_t {

                    __ATTR_ALIGN__(16) __m128 mstr[nstr*nmstr];
                    __ATTR_ALIGN__(16) __m128 pmstr[nstr*npmstr];
                    constexpr static int32_t NSTR  = nstr;
                    constexpr static int32_t NMSTR = nmstr;
                    constexpr static int32_t NPMSTR= npmstr; 
                    constexpr __m128 * __restrict mstr_beg()  { return (std::addressof(mstr[0]));}
                    constexpr __m128 * __restrict pmstr_beg() { return (std::addressof(pmstr[0]));}
                    float * __restrict mstr_fptr() { return ((float* __restrict)&mstr[0]);}
                    float * __restrict pmstr_fptr() { return ((float* __restrict)&pmstr[0]);}
                    constexpr int32_t mstr_size()  { return (NSTR*NMSTR);}
                    constexpr int32_t pmstr_size() { return (NSTR*NPMSTR);}
            };
            
            
            // !Coverage of every street (like: '1' for snow,'2' for mud, '3' for clay, ...etc)
            // number of streets
            template<int32_t nstr>
            struct SSECStrIx_t {

                 __ATTR_ALIGN__(16) __m128i cstr[nstr];
                 constexpr int32_t static NSTR = nstr;
                 constexpr __m128i * __restrict cstr_beg() { return (std::addressof(cstr[0]));}
                 int32_t * __restrict cstr_iptr() { return ((int32_t* __restrict)&cstr[0]);}
                 constexpr int32_t cstr_size()  { return (NSTR);}
            };
            
            
            // Percent of covered to non-covered portion of every street at each irradiated cell
            // Average thickness of each layer (cover) of every street at each irradiated cell
            // Thickness of cover along street (number of values) at each irradiated cell
             // number of streets
            template<int32_t nstr,int32_t npcstr,
                              int32_t natstr,int32_t ntcstr>
            struct SSECDStrR1x_t {

                    __ATTR_ALIGN__(16) __m128 pcstr[nstr*npcstr];
           
                    __ATTR_ALIGN__(16) __m128 atstr[nstr*natstr];
               
                    __ATTR_ALIGN__(16) __m128 tcstr[nstr*ntcstr];
                    constexpr static int32_t NSTR  = mstr;
                    constexpr static int32_t NPCSTR= npcstr;
                    constexpr static int32_t NATSTR= natstr;
                    constexpr static int32_t NTCSTR= ntcstr;
                    constexpr __m128 * __restrict pcstr_beg() { return (std::addressof(pcstr[0]));}
                    constexpr __m128 * __restrict atstr_beg() { return (std::addressof(atstr[0]));}
                    constexpr __m128 * __restrict tcstr_beg() { return (std::addressof(tcstr[0]));}
                    float * __restrict pcstr_fptr() { return ((float* __restrict)&pcstr[0]);}
                    float * __restrict atstr_fptr() { return ((float* __restrict)&atstr[0]);}
                    float * __restrict tcstr_fptr() { return ((float* __restrict)&tcstr[0]);}
                    constexpr int32_t pcstr_size() { return (nstr*npcstr);}
                    constexpr int32_t atstr_size() { return (nstr*natstr);}
                    constexpr int32_t tcstr_size() { return (nstr*ntcstr);}
            };
            
            
            // Mu values for 'clean' street interpolated along the street length at each irradiated cell
            // Eps for 'clean' street street length interpolated at each irradiated cell
             // number of streets
            template<int32_t nstr,
                     int32_t nmustr,int32_t nepstr> 
            struct SSEMEStr1C1x_t {       

                     __ATTR_ALIGN__(16) __m128  murstr[nstr*nmustr];
                   
                     __ATTR_ALIGN__(16) __m128 muistr[nstr*nmustr];
                  
                     __ATTR_ALIGN__(16) __m128 eprstr[nstr*nepstr];
                  
                     __ATTR_ALIGN__(16) __m128 epistr[nstr*nepstr];
                    constexpr static int32_t NSTR   = nstr;
                    constexpr static int32_t NMUSTR = nmustr;
                    constexpr static int32_t NEPSTR = nepstr; 
                    constexpr __m128 * __restrict murstr_beg() { return (std::addressof(murstr[0]));}
                    constexpr __m128 * __restrict muistr_beg() { return (std::addressof(muistr[0]));}
                    constexpr __m128 * __restrict eprstr_beg() { return (std::addressof(eprstr[0]));}
                    constexpr __m128 * __restrict epistr_beg() { return (std::addressof(epistr[0]));}
                    float * __restrict murstr_fptr() { return ((float* __restrict)&murstr[0]);}
                    float * __restrict muistr_fptr() { return ((float* __restrict)&muistr[0]);}
                    float * __restrict eprstr_fptr() { return ((float* __restrict)&eprstr[0]);}
                    float * __restrict epistr_fptr() { return ((float* __restrict)&epistr[0]);}
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
            template<int32_t nstr,
                     int32_t nmustr,int32_t nepstr>
            struct SSEMEStr2C1x_t {

                    __ATTR_ALIGN__(16) __m128  murstr[nstr*nmustr];
                    
                    __ATTR_ALIGN__(16) __m128  muistr[nstr*nmustr];
                   
                    __ATTR_ALIGN__(16) __m128  eprstr[nstr*nepstr];
                    
                    __ATTR_ALIGN__(16) __m128  epistr[nstr*nepstr];  
                    constexpr static int32_t NSTR   = nstr;
                    constexpr static int32_t NMUSTR = nmustr;
                    constexpr static int32_t NEPSTR = nepstr;
                    constexpr __m128 * __restrict murstr_beg() { return (std::addressof(murstr[0]));}
                    constexpr __m128 * __restrict muistr_beg() { return (std::addressof(muistr[0]));}
                    constexpr __m128 * __restrict eprstr_beg() { return (std::addressof(eprstr[0]));}
                    constexpr __m128 * __restrict epistr_beg() { return (std::addressof(epistr[0]));}
                    float * __restrict murstr_fptr() { return ((float* __restrict)&murstr[0]);}
                    float * __restrict muistr_fptr() { return ((float* __restrict)&muistr[0]);}
                    float * __restrict eprstr_fptr() { return ((float* __restrict)&eprstr[0]);}
                    float * __restrict epistr_fptr() { return ((float* __restrict)&epistr[0]);}
                    constexpr int32_t murstr_size() { return (nstr*nmustr);}
                    constexpr int32_t muistr_size() { return (nstr*nmustr);}
                    constexpr int32_t eprstr_size() { return (nstr*nepstr);}
                    constexpr int32_t epistr_size() { return (nstr*nepstr);}
            };
            
            
            // Street curvature parametric equation u-parameter
            // Street curvature parametric equation v-parameter
            // number of streets
            template<int32_t nstr,
                     int32_t nustr,int32_t nvstr>
            struct SSESCrvR1x_t {

                     __ATTR_ALIGN__(16) __m128 ustr[nstr*nustr];
                    
                     __ATTR_ALIGN__(16) __m128 vstr[nstr*nvstr];
                     constexpr static int32_t NSTR = nstr;
                     constexpr static int32_t NUSTR= nustr;
                     constexpr static int32_t NVSTR= nvstr;
                     constexpr __m128 * __restrict ustr_beg() { return (std::addressof(ustr[0]));}
                     constexpr __m128 * __restrict vstr_beg() { return (std::addressof(vstr[0]));}
                     float * __restrict ustr_fptr() { return ((float* __restrict)&ustr[0]);}
                     float * __restrict vstr_fptr() { return ((float* __restrict)&vstr[0]);}
                     constexpr int32_t ustr_size() { return (nstr*nustr);}
                     constexpr int32_t vstr_size() { return (nstr*nvstr);}
            };
            
            
            //  Street surface normal vectors x-components 
            //  along the street length at each irradiated cell
            //  Street surface normal vectors y-components 
            //  along the street length at each irradiated cell
            //  Street surface normal vectors z-components 
            //  along the street length at each irradiated cell
             // number of streets
            template<int32_t nstr,
                     int32_t nx,int32_t ny,int32_t,nz>
            struct SSESNrmR1x_t {

                      __ATTR_ALIGN__(16) __m128  nvx[nstr*nx];
                   
                      __ATTR_ALIGN__(16) __m128  nvy[nstr*ny];
                   
                      __ATTR_ALIGN__(16) __m128  nvz[nstr*nz];
                    constexpr static int32_t NSTR = nstr;
                    constexpr static int32_t NX   = nx;
                    constexpr static int32_t NY   = ny;
                    constexpr static int32_t NZ   = nz;
                    constexpr __m128 * __restrict nvx_beg() { return (std::addressof(nvx[0]));}
                    constexpr __m128 * __restrict nvy_beg() { return (std::addressof(nvy[0]));}
                    constexpr __m128 * __restrict nvz_beg() { return (std::addressof(nvz[0]));}
                    float * __restrict nvx_fptr() { return ((float* __restrict)&nvx[0]);}
                    float * __restrict nvy_fptr() { return ((float* __restrict)&nvy[0]);}
                    float * __restrict nvz_fptr() { return ((float* __restrict)&nvz[0]);}
                    constexpr int32_t nvx_size() { return (nstr*nx);}
                    constexpr int32_t nvy_size() { return (nstr*ny);}
                    constexpr int32_t nvz_size() { return (nstr*nz);}
            };
            
            
            // latitude   values (deg), per street length (at irradiance point)
            // longtitude values (deg), per street length (at irradiance point)
           // number of streets
            template<int32_t nstr,
                     int32_t nlon,int32_t nlat>
            struct SSESIRCDR1x_t {

                    __ATTR_ALIGN__(16) __m128 irlon[nstr*nlon];
                  
                    __ATTR_ALIGN__(16) __m128 irlat[nstr*nlat]; 
                    constexpr static int32_t NSTR   = nstr;
                    constexpr static int32_t NLON   = nlon;
                    constexpr static int32_t NLAT   = nlat;
                    constexpr __m128 * __restrict irlon_beg() { return (std::addressof(irlon[0]));}
                    float * __restrict irlon_fptr() { return ((float* __restrict)&irlon[0]);}
                    constexpr int32_t irlon_size() { return (NSTR*NLON);}
                    constexpr __m128 * __restrict irlat_beg() { return (std::addressof(irlat[0]));}
                    float * __restrict irlat_fptr() { return ((float* __restrict)&irlat[0]);}
                    constexpr int32_t irlat_size() { return (NSTR*NLAT);}
                   
            };
            
            
           // latitude   values (rad), per street length (at irradiance point)
           // longtitude values (rad), per street length (at irradiance point)
           // number of streets
           template<int32_t nstr,
                    int32_t nlon,int32_t nlat>
           struct SSESIRCRR1x_t {

                    __ATTR_ALIGN__(16) __m128 irlon[nstr*nlon];

                    __ATTR_ALIGN__(16) __m128 irlat[nstr*nlat];
                    constexpr static int32_t NSTR = nstr;
                    constexpr static int32_t NLON = nlon;
                    constexpr static int32_t NLAT = nlat; 
                    constexpr __m128 * __restrict irlon_beg() { return (std::addressof(irlon[0]));}
                    float * __restrict irlon_fptr() { return ((float* __restrict)&irlon[0]);}
                    constexpr int32_t irlon_size() { return (NSTR*NLON);}
                    constexpr __m128 * __restrict irlat_beg() { return (std::addressof(irlat[0]));}
                    float * __restrict irlat_fptr() { return ((float* __restrict)&irlat[0]);}
                    constexpr int32_t irlat_size() { return (NSTR*NLAT);}
                
           };
           
           
           // latitude   values (deg), of the building area (at irradiance point)
           // longtitude values (deg), of the building area (at irradiance point)
            // number of buildings
            template<int32_t nbld,
                     int32_t nlon,int32_t nlat>
            struct SSEBIRCDR1x_t {

                     __ATTR_ALIGN__(16) __m128 irlon[nbld*nlon];
                  
                     __ATTR_ALIGN__(16) __m128 irlat[nbld*nlat]; 
                     constexpr static int32_t NBLD = nbld;
                     constexpr static int32_t NLON = nlon;
                     constexpr static int32_t NLAT = nlat;
                     constexpr __m128 * __restrict irlon_beg() { return (std::addressof(irlon[0]));}
                     float * __restrict irlon_fptr() { return ((float* __restrict)&irlon[0]);}
                     constexpr int32_t irlon_size() { return (NBLD*NLON);}
                     constexpr __m128 * __restrict irlat_beg() { return (std::addressof(irlat[0]));}
                     float * __restrict irlat_fptr() { return ((float* __restrict)&irlat[0]);}
                     constexpr int32_t irlat_size() { return (NBLD*NLAT);}
            };
            
            
           // latitude   values (rad), of the building area (at irradiance point)
           // longtitude values (rad), of the building area (at irradiance point)
            // number of buildings
           template<int32_t nbld,
                    int32_t nlon,int32_t nlat>
           struct SSEBIRCRR1x_t {

                     __ATTR_ALIGN__(16) __m128  irlon[nbld*nlon];
                 
                     __ATTR_ALIGN__(16) __m128 irlat[nbld*nlat];
                     constexpr static int32_t NBLD = nbld;
                     constexpr static int32_t NLON = nlon;
                     constexpr static int32_t NLAT = nlat;
                     constexpr __m128 * __restrict irlon_beg() { return (std::addressof(irlon[0]));}
                     float * __restrict irlon_fptr() { return ((float* __restrict)&irlon[0]);}
                     constexpr int32_t irlon_size() { return (NBLD*NLON);}
                     constexpr __m128 * __restrict irlat_beg() { return (std::addressof(irlat[0]));}
                     float * __restrict irlat_fptr() { return ((float* __restrict)&irlat[0]);}
                     constexpr int32_t irlat_size() { return (NBLD*NLAT);}
           }; 
           
           
           // Urban area height map (at single building resolution)
           template<int32_t nx,int32_t ny>
           struct SSEUHMapR1x_t {

                  __ATTR_ALIGN__(16) __m128 hmap[nx*ny];
                  constexpr static int32_t NX = nx;
                  constexpr static int32_t NY = ny;
                  constexpr __m128 * __restrict hmap_beg() { return (std::addressof(hmap[0]));}
                  float * __restrict hmap_fptr() { return ((float* __restrict)&hmap[0]);}
                  constexpr int32_t hmap_size() { return (NX*NY);}
           };
           
           
           // Urban area height map (at single building resolution) -- 1st derivative
           template<int32_t nx,int32_t ny>
           struct SSEUHDxDyR1x_t {
               
                    __ATTR_ALIGN__(16) __m128 hdxdy[nx*ny];

                   constexpr static int32_t NX = nx;
                   constexpr static int32_t NY = ny;  
                   constexpr __m128 * __restrict hdxdy_beg() { return (std::addressof(hdxdy[0]));}
                   float * __restrict hdxdy_fptr() { return ((float* __restrict)&hdxdy[0]);}
                   constexpr int32_t hdxdy_size() { return (NX*NY);}              
           };
           
           
           // Urban area height map (at single building resolution) -- gradient x-component
           // Urban area height map (at single building resolution) -- gradient y-component
           template<int32_t nx,int32_t ny>
           struct SSEUHGradR1x_t {
                 
                  __ATTR_ALIGN__(16) __m128 uhgx[nx];
                  
                  __ATTR_ALIGN__(16) __m128 uhgy[ny];
                  constexpr static int32_t NX = nx;
                  constexpr static int32_t NY = ny;
                  constexpr __m128 * __restrict uhgx_beg() { return (std::addressof(uhgx[0]));}
                  float * __restrict uhgx_fptr() { return ((float* __restrict)&uhgx[0]);}
                  constexpr int32_t uhgx_size() { return (NX);}
                  constexpr __m128 * __restrict uhgy_beg() { return (std::addressof(uhgy[0]));}
                  float * __restrict uhgy_fptr() { return ((float* __restrict)&uhgy[0]);}
                  constexpr int32_t uhgy_size() { return (NY);}  
           };
           
           
           // Smoothing and approximating curve for linearly-piecewise height function (x-coordinate)
           // Smoothing and approximating curve for linearly-piecewise height function (y-coordinate)
           template<int32_t nx,int32_t ny>
           struct SSEXYSMBHR1x_t {
                 
                   __ATTR_ALIGN__(16) __m128 xsmbh[nx];
                  
                   __ATTR_ALIGN__(16) __m128 ysmbh[ny];
                   constexpr static int32_t NX = nx;
                   constexpr static int32_t NY = ny;
                   constexpr __m128 * __restrict xsmbh_beg() { return (std::addressof(xsmbh[0]));}
                   float * __restrict xsmbh_fptr() { return ((float* __restrict)&xsmbh[0]);}
                   constexpr int32_t xsmbh_size() { return (NX);}
                   constexpr __m128 * __restrict ysmbh_beg() { return (std::addressof(ysmbh[0]));}
                   float * __restrict ysmbh_fptr() { return ((float* __restrict)&ysmbh[0]);}
                   constexpr int32_t ysmbh_size() { return (NY);}  
           };
           
           
           // Empty space in-between of buildings (per single column) x number columns
           template<int32_t ncols,int32_t nval>
           struct SSEESBBI1x_t {
                  
                   __ATTR_ALIGN__(16) __m128i esbb[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128i * __restrict esbb_beg() { return (std::addressof(esbb[0]));}
                  int32_t  * __restrict esbb_iptr() { return ((int32_t* __restrict)&esbb[0]);}
                  constexpr int32_t esbb_size() { return (NCOLS*NVAL);}
           };
           
           
           // An area values of in-between buildings empty spaces (per single column) x number columns
           template<int32_t ncols,int32_t nval> 
           struct SSEAESBBR1x_t {
                  
                   __ATTR_ALIGN__(16) __m128 aesbb[ncols*nval];
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;
                   constexpr __m128 * __restrict aesbb_beg() { return (std::addressof(aesbb[0]));}
                   float * __restrict aesbb_fptr() { return ((float* __restrict)&aesbb[0]);}
                   constexpr int32_t aesbb_size() { return (NCOLS*NVAL);}
           };
           
           
           // An area values of each building (per single building column) x number columns
           template<int32_t ncols,int32_t nval>
           struct SSEABCR1x_t {
                
                   __ATTR_ALIGN__(16) __m128  abc[ncols*nvals];
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;
                   constexpr __m128 * __restrict abc_beg() { return (std::addressof(abc[0]));}
                   float * __restrict abc_fptr() { return ((float* __restrict)&abc[0]);}
                   constexpr int32_t abc_size() { return (NCOLS*NVAL);}
           };
           
           
           // Number of south-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct SSESWPCI1x_t {
                  
                   __ATTR_ALIGN__(16) __m128i swpc[ncols];
                   constexpr static int32_t NCOLS = ncols;
                   constexpr __m128i * __restrict swpc_beg() { return (std::addressof(swpc[0]));}
                   int32_t * __restrict swpc_iptr() { return ((int32_t* __restrict)&swpc[0]);}
                   constexpr int32_t swpc_size() { return (NCOLS);}
           };
           
           
           // Number of east-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct SSEEWPCI1x_t {
                  
                 __ATTR_ALIGN__(16) __m128i ewpc[ncols];
                 constexpr static int32_t NCOLS = ncols;
                 constexpr __m128i * __restrict ewpc_beg() { return (std::addressof(ewpc[0]));}
                 int32_t * __restrict ewpc_iptr() { return ((int32_t* __restrict)&ewpc[0]);}
                 constexpr int32_t ewpc_size() { return (NCOLS);} 
           };
           
           
           // Number of west-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct SSEWWPCI1x_t {
                  
                 __ATTR_ALIGN__(16) __m128i  wwpc[ncols];
                 constexpr static int32_t NCOLS = ncols;
                 constexpr __m128i * __restrict wwpc_beg() { return (std::addressof(wwpc[0]));}
                 int32_t * __restrict wwpc_iptr() { return ((int32_t* __restrict)&wwpc[0]);}
                 constexpr int32_t wwpc_size() { return (NCOLS);}
           };
           
           
           // Number of north-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct SSENWPCI1x_t {
                
                __ATTR_ALIGN__(16) __m128i  nwpc[ncols];
                constexpr static int32_t NCOLS = ncols;
                constexpr __m128i * __restrict nwpc_beg() { return (std::addressof(nwpc[0]));}
                int32_t * __restrict nwpc_iptr() { return ((int32_t* __restrict)&nwpc[0]);}
                constexpr int32_t nwpc_size() { return (NCOLS);} 
           };
           
           
           // Number of building roofs per each column x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEBRPCI1x_t {
                 
                 __ATTR_ALIGN__(16) __m128i brpc[ncols*nval]; 
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr __m128i * __restrict brpc_beg() { return (std::addressof(brpc[0]));}
                 int32_t * __restrict brpc_iptr() { return ((int32_t* __restrict)&brpc[0]);}
                 constexpr int32_t brpc_size() { return (NCOLS*NVAL);}
           };
           
           
           //  An area of every building [flat] roof (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEBRAPCR1x_t {
               
                  __ATTR_ALIGN__(16) __m128 brapc[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict brapc_beg() { return (std::addressof(brapc[0]));}
                  float * __restrict brapc_fptr() { return ((float* __restrict)&brapc[0]);}
                  constexpr int32_t brapc_size() { return (NCOLS*NVAL);}
           };
           
           
           // Number of angled roof -- south facing roof wall (per each column) x number of columns
           template<int32_t ncols>
           struct SSESRWCI1x_t {
             
                  __ATTR_ALIGN__(16) __m128i srwc[ncols];
                  constexpr static int32_t NCOLS = ncols;
                   constexpr __m128i * __restrict srwc_beg() { return (std::addressof(srwc[0]));}
                  int32_t * __restrict srwc_iptr() { return ((int32_t* __restrict)&srwc[0]);}
                  constexpr int32_t srwc_size() { return (NCOLS);}
           };
           
           
           // Number of angled roof -- east facing roof wall (per each column) x number of columns
           template<int32_t ncols>
           struct SSEERWCI1x_t {
                
                  __ATTR_ALIGN__(16) __m128i  erwc[ncols];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr __m128i * __restrict erwc_beg() { return (std::addressof(erwc[0]));}
                  int32_t * __restrict erwc_iptr() { return ((int32_t* __restrict)&erwc[0]);}
                  constexpr int32_t erwc_size() { return (NCOLS);}
           };
           
           
           // Number of angled roof -- west facing roof wall (per each column)  x number of columns
           template<int32_t ncols>
           struct SSEWRWCI1x_t {
            
                __ATTR_ALIGN__(16) __m128i  wrwc[ncols];
                constexpr static int32_t NCOLS = ncols;
                constexpr __m128i * __restrict wrwc_beg() { return (std::addressof(wrwc[0]));}
                int32_t * __restrict wrwc_iptr() { return ((int32_t* __restrict)&wrwc[0]);}
                constexpr int32_t wrwc_size() { return (NCOLS);}  
           };
           
           
           // Number angled roof -- north facing roof wall (per each column) x number of columns
           template<int32_t ncols>
           struct SSENRWCI1x_t {
                
                __ATTR_ALIGN__(16) __m128i nrwc[ncols];
                constexpr static int32_t NCOLS = ncols;
                constexpr __m128i * __restrict nrwc_beg() { return (std::addressof(nrwc[0]));}
                int32_t * __restrict nrwc_iptr() { return ((int32_t* __restrict)&nrwc[0]);}
                constexpr int32_t nrwc_size() { return (NCOLS);}
           };
           
           
           // An angled roof inclination (deg) -- south facing roof wall 
           // (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEIDSRWR1x_t {
               
                  __ATTR_ALIGN__(16) __m128 idsrw[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict idsrw_beg() { return (std::addressof(idsrw[0]));}
                  float * __restrict idsrw_fptr() { return ((float* __restrict)&idsrw[0]);}
                  constexpr int32_t idsrw_size() { return (NCOLS*NVAL);}
           };
           
           
           // An angled roof inclination (deg) -- east facing roof wall 
           // (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEIDERWR1x_t {
              
                  __ATTR_ALIGN__(16) __m128 iderw[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict iderw_beg() { return (std::addressof(iderw[0]));}
                  float * __restrict iderw_fptr() { return ((float* __restrict)&iderw[0]);}
                  constexpr int32_t iderw_size() { return (NCOLS*NVAL);}
           };
           
           
           // An angled roof inclination (deg) -- west facing roof wall 
           // (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEIDWRWR1x_t {
           
                  __ATTR_ALIGN__(16) __m128 idwrw[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict idwrw_beg() { return (std::addressof(idwrw[0]));}
                  float * __restrict idwrw_fptr() { return ((float* __restrict)&idwrw[0]);}
                  constexpr int32_t idwrw_size() { return (NCOLS*NVAL);}
           };
           
           
           // An angled roof inclination (rad) -- south facing roof wall 
           // (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEIDNRWR1x_t {
               
                   __ATTR_ALIGN__(16) __m128  idnrw[ncols*nval]; 
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;
                   constexpr __m128 * __restrict idnrw_beg() { return (std::addressof(idnrw[0]));}
                   float * __restrict idnrw_fptr() { return ((float* __restrict)&idnrw[0]);}
                   constexpr int32_t idnrw_size() { return (NCOLS*NVAL);}
           };
           
           
           //  An angled roof inclination (rad) -- south facing roof wall 
           //  (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEIRSRWR1x_t {
                
                __ATTR_ALIGN__(16) __m128 irsrw[ncols*nval]; 
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr __m128 * __restrict irsrw_beg() { return (std::addressof(irsrw[0]));}
                float * __restrict irsrw_fptr() { return ((float* __restrict)&irsrw[0]);}
                constexpr int32_t irsrw_size() { return (NCOLS*NVAL);}
           };
           
           
           //  An angled roof inclination (rad) -- east facing roof wall 
           //  (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEIRERWR1x_t {
                 
                    __ATTR_ALIGN__(16) __m128 irerw[ncols*nval];
                    constexpr static int32_t NCOLS = ncols;
                    constexpr static int32_t NVAL  = nval;
                    constexpr __m128 * __restrict irerw_beg() { return (std::addressof(irerw[0]));}
                    float * __restrict irerw_fptr() { return ((float* __restrict)&irerw[0]);}
                    constexpr int32_t irerw_size() { return (NCOLS*NVAL);}
           };
           
           
           //  An angled roof inclination (rad) -- north facing roof wall 
           //  (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEIRNRWR1x_t {
                   
                    __ATTR_ALIGN__(16) __m128 irnrw[ncols*nval]; 
                    constexpr static int32_t NCOLS = ncols;
                    constexpr static int32_t NVAL  = nval;
                    constexpr __m128 * __restrict irnrw_beg() { return (std::addressof(irnrw[0]));}
                    float * __restrict irnrw_fptr() { return ((float* __restrict)&irnrw[0]);}
                    constexpr int32_t irnrw_size() { return (NCOLS*NVAL);}
           };
           
           
           // An angled roof inclination surface area -- south facing roof wall 
           // (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEISRAR1x_t {
                
                  __ATTR_ALIGN__(16) __m128  isra[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval; 
                  constexpr __m128 * __restrict isra_beg() { return (std::addressof(isra[0]));}
                  float * __restrict isra_fptr() { return ((float* __restrict)&isra[0]);}
                  constexpr int32_t isra_size() { return (NCOLS*NVAL);}
           };
           
           
           // An angled roof inclination surface area -- west facing roof wall 
           // (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEIWRAR1x_t {
            
                   __ATTR_ALIGN__(16) __m128 iwra[ncols*nval];
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;  
                   constexpr __m128 * __restrict iwra_beg() { return (std::addressof(iwra[0]));}
                   float * __restrict iwra_fptr() { return ((float* __restrict)&iwra[0]);}
                   constexpr int32_t iwra_size() { return (NCOLS*NVAL);}
           };
           
           
           // An angled roof inclination surface area -- east facing roof wall 
           // (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEIERAR1x_t {
            
                  __ATTR_ALIGN__(16) __m128 iera[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval; 
                  constexpr __m128 * __restrict iera_beg() { return (std::addressof(iera[0]));}
                  float * __restrict iera_fptr() { return ((float* __restrict)&iera[0]);}
                  constexpr int32_t iera_size() { return (NCOLS*NVAL);}  
           };
           
           
           // An angled roof inclination surface area -- north facing roof wall 
           // (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SSEINRAR1x_t {
               
                    __ATTR_ALIGN__(16) __m128 inra[ncols*nval];
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;  
                   constexpr __m128 * __restrict inra_beg() { return (std::addressof(inra[0]));}
                   float * __restrict inra_fptr() { return ((float* __restrict)&inra[0]);}
                   constexpr int32_t inra_size() { return (NCOLS*NVAL);} 
           };
           
           
           // South wall upper-facing edge inclination (rad) -- 
           // (per each column)  x number of columns
          template<int32_t ncols,int32_t nval>
          struct SSESWUER1x_t {
                 
                   __ATTR_ALIGN__(16) __m128 swue[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval; 
                  constexpr __m128 * __restrict swue_beg() { return (std::addressof(swue[0]));}
                  float * __restrict swue_fptr() { return ((float* __restrict)&swue[0]);}
                  constexpr int32_t swue_size() { return (NCOLS*NVAL);} 
          };
          
          
          // East wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<int32_t ncols,int32_t nval>
          struct SSEEWUER1x_t {
                 
                 __ATTR_ALIGN__(16) __m128 ewue[ncols*nval]; 
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;  
                 constexpr __m128 * __restrict ewue_beg() { return (std::addressof(ewue[0]));}
                 float * __restrict ewue_fptr() { return ((float* __restrict)&ewue[0]);}
                 constexpr int32_t ewue_size() { return (NCOLS*NVAL);} 
          };
          
          
          // West wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<int32_t ncols,int32_t nval>
          struct SSEWWUER1x_t {
              
                  __ATTR_ALIGN__(16) __m128 wwue[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;  
                  constexpr __m128 * __restrict wwue_beg() { return (std::addressof(wwue[0]));}
                  float * __restrict wwue_fptr() { return ((float* __restrict)&wwue[0]);}
                  constexpr int32_t wwue_size() { return (NCOLS*NVAL);} 
          };
          
          
          // North wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<int32_t ncols,int32_t nval>
          struct SSENWUER1x_t {
                
                  __ATTR_ALIGN__(16) __m128 nwue[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;    
                  constexpr __m128 * __restrict nwue_beg() { return (std::addressof(nwue[0]));}
                  float * __restrict nwue_fptr() { return ((float* __restrict)&nwue[0]);}
                  constexpr int32_t nwue_size() { return (NCOLS*NVAL);} 
          };
          
          
          // Shared right edge between the south wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          template<int32_t ncols,int32_t nval>
          struct SSESEWER1x_t {
                  
                  __ATTR_ALIGN__(16) __m128 sewe[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;   
                  constexpr __m128 * __restrict sewe_beg() { return (std::addressof(sewe[0]));}
                  float * __restrict sewe_fptr() { return ((float* __restrict)&sewe[0]);}
                  constexpr int32_t sewe_size() { return (NCOLS*NVAL);}
          };
          
          // Shared left edge between the south wall and west wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          template<int32_t ncols,int32_t nval>
          struct SWWER1x_t {
               
                  __ATTR_ALIGN__(16) __m128 swwe[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;  
                  constexpr __m128 * __restrict swwe_beg() { return (std::addressof(swwe[0]));}
                  float * __restrict swwe_fptr() { return ((float* __restrict)&swwe[0]);}
                  constexpr int32_t swwe_size() { return (NCOLS*NVAL);}
          };
          
          
           // Shared right edge between the north wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct NWEER1x_t {
               
                   __ATTR_ALIGN__(16) __m128 nwee[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval; 
                  constexpr __m128 * __restrict nwee_beg() { return (std::addressof(nwee[0]));}
                  float * __restrict nwee_fptr() { return ((float* __restrict)&nwee[0]);}
                  constexpr int32_t nwee_size() { return (NCOLS*NVAL);}
           };
           
           
           // Shared right edge between the north wall and west wall inclination (rad) 
           // ! -- (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct NWWER1x_t {
              
                   __ATTR_ALIGN__(16) __m128 nwwe[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval; 
                  constexpr __m128 * __restrict nwwe_beg() { return (std::addressof(nwwe[0]));}
                  float * __restrict nwwe_fptr() { return ((float* __restrict)&nwwe[0]);}
                  constexpr int32_t nwwe_size() { return (NCOLS*NVAL);}
           };
           
           
           // Simple cell-based mesh
           template<int32_t nL>
           struct CellMeshR1x_t {
                  
                  // Coordinates (x,y,z) of the center
                  // of Lth cell.
                  __ATTR_ALIGN__(16) __m128 cx[nL];
                  __ATTR_ALIGN__(16) __m128 cy[nL];
                  __ATTR_ALIGN__(16) __m128 cz[nL];
                  __ATTR_ALIGN__(16) __m128 dv[nL];
                  // (X,Y,Z) dimensions of the Ith
                  // rectangular volume cell (this is needed for
                  // the numerical integration)
                  __ATTR_ALIGN__(16) __m128 dx[nL];
                  __ATTR_ALIGN__(16) __m128 dy[nL];
                  __ATTR_ALIGN__(16) __m128 dz[nL];
                  // Number of divisions along the x,y,z
                  int32_t ndiv[3];
                  // Compute numerical integration.
                  bool nint;
                 
                  
           };
           
           
           // South walls surface area (for every building, per column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct SWSAR1x_t {
                
                   __ATTR_ALIGN__(16) __m128 swsa[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict swsa_beg() { return (std::addressof(swsa[0]));}
                   float * __restrict swsa_fptr() { return ((float* __restrict)&swsa[0]);}
                  constexpr int32_t swsa_size() { return (NCOLS*NVAL);}
           };
           
           
           // East walls surface area (for every building, per column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct EWSAR1x_t {
                  
                   __ATTR_ALIGN__(16) __m128 ewsa[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict ewsa_beg() { return (std::addressof(ewsa[0]));}
                  float * __restrict ewsa_fptr() { return ((float* __restrict)&ewsa[0]);}
                  constexpr int32_t ewsa_size() { return (NCOLS*NVAL);}
           };
           
           
           // West walls surface area (for every building, per column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct WWSAR1x_t {
                 
                  __ATTR_ALIGN__(16) __m128  wwsa[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict wwsa_beg() { return (std::addressof(wwsa[0]));}
                  float * __restrict wwwsa_fptr() { return ((float* __restrict)&wwsa[0]);}
                  constexpr int32_t wwsa_size() { return (NCOLS*NVAL);}
           };
           
           
           // North walls surface area (for every building, per column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct NWSAR1x_t {
                     
                  __ATTR_ALIGN__(16) __m128 nwsa[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict begin() { return (std::addressof(nwsa[0]));}
                  float * __restrict nwwsa_fptr() { return ((float* __restrict)&nwsa[0]);} 
                  constexpr int32_t size() { return (ncols*nval)};  
                 
           };
           
           
           // South walls moist/non moist logical (per column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct MNMSWB1x_t {
                    
            
            
                   bool  mnmsw[ncols*nval]; 
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;
                   constexpr bool * mnmsw_beg() { return (std::addressof(mnmsw[0]));}
                   constexpr int32_t mnmsw_size() { return (NCOLS*NVAL);}
           };
           
           
            // East walls moist/non moist logical (per column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct MNMEWB1x_t {
                                 
                    
                   bool  mnmew[ncols*nval]; 
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;
                   constexpr bool * mnmew_beg() { return (std::addressof(mnmew[0]));}
                   constexpr int32_t mnmew_size() { return (NCOLS*NVAL);}
           };
           
           
             // West walls moist/non moist logical (per column) x number of columns
           template<int32_t ncols,int32_t nval>  
           struct MNMWWB1x_t {
           
                     
                   bool  mnmww[ncols*nval];
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval; 
                   constexpr bool * mnmww_beg() { return (std::addressof(mnmww[0]));}
                   constexpr int32_t mnmww_size() { return (NCOLS*NVAL);}
           };
           
           
             // North walls moist/non moist logical (per column) x number of columns
           template<int32_t ncols,int32_t nval>  
           struct MNMNWB1x_t {
                     
                   bool  mnmnw[ncols*nval]; 
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;
                   constexpr bool * mnmnw_beg() { return (std::addressof(mnmnw[0]));}
                   constexpr int32_t mnmnw_size() { return (NCOLS*NVAL);}
           };
           
           
           // ! The values describing the ratio (percentage) of south wall 
                            // ! moisture to dryness (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct MDSWRR1x_t {
                 
                    __ATTR_ALIGN__(16) __m128  mdswr[ncols*nval];
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;
                   constexpr __m128 * __restrict mdswr_beg() { return (std::addressof(mdswr[0]));}
                   float * __restrict mdswr_fptr() { return ((float* __restrict)&mdswr[0]);} 
                   constexpr int32_t mdswr_size() { return (NCOLS*NVAL);}
           };
           
           
           // ! The values describing the ratio (percentage) of east wall 
                            // ! moisture to dryness (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct MDEWRR1x_t {
                    
                    __ATTR_ALIGN__(16) __m128    mdewr[ncols*nval];                 
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;
                   constexpr __m128 * __restrict mdewr_beg() { return (std::addressof(mdewr[0]));}
                   float * __restrict mdewr_fptr() { return ((float* __restrict)&mdewr[0]);} 
                   constexpr int32_t mdewr_size() { return (NCOLS*NVAL);}
                  
           };
           
           
           //  The values describing the ratio (percentage) of west wall 
                             //! moisture to dryness (per each column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct MDWWRR1x_t {
                    
                   __ATTR_ALIGN__(16) __m128    mdwwr[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict mdwwr_beg() { return (std::addressof(mdwwr[0]));}
                  float * __restrict mdwwr_fptr() { return ((float* __restrict)&mdwwr[0]);} 
                  constexpr int32_t mdwwr_size() { return (NCOLS*NVAL);}
           };  
           
           
           //  The values describing the ratio (percentage) of north wall 
                             //! moisture to dryness (per each column) x number of columns              
           template<int32_t ncols,int32_t nval>
           struct MDNWRR1x_t {
                
                   __ATTR_ALIGN__(16) __m128 mdnwr[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict mdnwr_beg() { return (std::addressof(mdnwr[0]));}
                  float * __restrict mdnwr_fptr() { return ((float* __restrict)&mdnwr[0]);} 
                  constexpr int32_t mdnwr_size() { return (NCOLS*NVAL);}
           };   
           
           
           // The logical values of flat roof moistness (being either moist or dry) 
                              // ! (per column) x number of columns    
          template<int32_t ncols,int32_t nval>            
          struct MDRB1x_t {
                 
                 
                 bool  mdr[ncols*nval];
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr bool * mdr_beg() { return (std::addressof(mdr[0]));}
                 constexpr int32_t mdr_size() { return (NCOLS*NVAL);}
          }; 
          
          
          // The values describing the ratio (percentage) of flat roof moisture to dryness 
                              // ! (per each column) x number of columns 
          template<int32_t ncols,int32_t nval>
          struct MDRRR1x_t {
             
                __ATTR_ALIGN__(16) __m128   mdrr[ncols*nval]; 
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr __m128 * __restrict mdrr_beg() { return (std::addressof(mdrr[0]));}
                float * __restrict mdrr_fptr() { return ((float* __restrict)&mdrr[0]);} 
                constexpr int32_t mdrr_size() { return (NCOLS*NVAL);}
          };
          
          
          // The values describing the surface of moist part of the flat roof 
                               // ! (per each column) x number of columns
          template<int32_t ncols,int32_t nval>
          struct MPFRR1x_t {
                 
        
                 __ATTR_ALIGN__(16) __m128    mpfr[ncols*nval]; 
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr __m128 * __restrict mpfr_beg() { return (std::addressof(mpfr[0]));}
                 float * __restrict mpfr_fptr() { return ((float* __restrict)&mpfr[0]);} 
                 constexpr int32_t mpfr_size() { return (NCOLS*NVAL);}
          };     
          
          
          // The values describing the surface of dry part of the flat roof 
                               // ! (per each column) x number of columns  
          template<int32_t ncols,int32_t nval>
          struct DPFRR1x_t {
                 
                   __ATTR_ALIGN__(16) __m128  dpfr[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict dpfr_beg() { return (std::addressof(dpfr[0]));}
                  float * __restrict dpfr_fptr() { return ((float* __restrict)&dpfr[0]);} 
                  constexpr int32_t dpfr_size() { return (NCOLS*NVAL);}
          };   
          
          
          //  The values describing the surface of moist part of the south wall 
                                // ! (per each column) x number of columns   
          template<int32_t ncols,int32_t nval>
          struct MPSWR1x_t {
                
                  __ATTR_ALIGN__(16) __m128   mpsw[ncols*nval]; 
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr __m128 * __restrict mpsw_beg() { return (std::addressof(mpsw[0]));}
                 float * __restrict mpsw_fptr() { return ((float* __restrict)&mpsw[0]);} 
                 constexpr int32_t mpsw_size() { return (NCOLS*NVAL);}
          };  
          
          
          //  The values describing the surface of dry part of the south wall 
                               //  ! (per each column) x number of columns  
          template<int32_t ncols,int32_t nval>
          struct DPSWR1x_t {
               
                   __ATTR_ALIGN__(16) __m128  dpsw[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict dpsw_beg() { return (std::addressof(dpsw[0]));}
                  float * __restrict dpsw_fptr() { return ((float* __restrict)&dpsw[0]);} 
                  constexpr int32_t dpsw_size() { return (NCOLS*NVAL);}
          };   
          
          
          //  The values describing the surface of moist part of the east wall 
                                // ! (per each column) x number of columns
          template<int32_t ncols,int32_t nval>
          struct MPEWR1x_t {
                   
                   __ATTR_ALIGN__(16) __m128 mpew[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict mpew_beg() { return (std::addressof(mpew[0]));}
                  float * __restrict mpew_fptr() { return ((float* __restrict)&mpew[0]);} 
                  constexpr int32_t mpew_size() { return (NCOLS*NVAL);}
          };   
          
         // The values describing the surface of dry part of the east wall 
                                // ! (per each column) x number of columns 
          template<int32_t ncols,int32_t nval>
          struct DPEWR1x_t {
                 
                  __ATTR_ALIGN__(16) __m128  dpew[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict dpew_beg() { return (std::addressof(dpew[0]));}
                  float * __restrict dpew_fptr() { return ((float* __restrict)&dpew[0]);} 
                  constexpr int32_t dpew_size() { return (NCOLS*NVAL);}
          };  
          
          
         // The values describing the surface of moist part of the west wall 
                                 //! (per each column) x number of columns
          template<int32_t ncols,int32_t nval>
          struct MPWWR1x_t {
                 
                   __ATTR_ALIGN__(16) __m128  mpww[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict mpww_beg() { return (std::addressof(mpww[0]));}
                  float * __restrict mpww_fptr() { return ((float* __restrict)&mpww[0]);} 
                  constexpr int32_t mpww_size() { return (NCOLS*NVAL);}
          }; 
          
          
        //  The values describing the surface of dry part of the west wall 
                                 //! (per each column) x number of columns 
          template<int32_t ncols,int32_t nval>
          struct DPWWR1x_t {
                 
                   __ATTR_ALIGN__(16) __m128 dpww[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict dpww_beg() { return (std::addressof(dpww[0]));}
                  float * __restrict dpww_fptr() { return ((float* __restrict)&dpww[0]);} 
                  constexpr int32_t dpww_size() { return (NCOLS*NVAL);}
          };  
          
          
        // The values describing the surface of moist part of the north wall 
                                 //! (per each column) x number of columns
         template<int32_t ncols,int32_t nval>
         struct MPNWR1x_t {
              
                   __ATTR_ALIGN__(16) __m128  mpnw[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict mpnw_beg() { return (std::addressof(mpnw[0]));}
                  float * __restrict mpnw_fptr() { return ((float* __restrict)&mpnw[0]);} 
                  constexpr int32_t mpnw_size() { return (NCOLS*NVAL);}
          }; 
          
          
         //  The values describing the surface of dry part of the north wall 
                                // ! (per each column) x number of columns
         template<int32_t ncols,int32_t nval>
         struct DPNWR1x_t {
                
                 __ATTR_ALIGN__(16) __m128  dpnw[ncols*nval]; 
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr __m128 * __restrict dpnw_beg() { return (std::addressof(dpnw[0]));}
                 float * __restrict dpnw_fptr() { return ((float* __restrict)&dpnw[0]);} 
                 constexpr int32_t dpnw_size() { return (NCOLS*NVAL);}
          }; 
          
          
          // The values describing the surface of moist part of the angled south roof wall
                               // ! (per each column) x number of columns
          template<int32_t ncols,int32_t nval>
          struct MPSARR1x_t {
                  
                  __ATTR_ALIGN__(16) __m128 mpsar[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict mpsar_beg() { return (std::addressof(mpsar[0]));}
                  float * __restrict mpsar_fptr() { return ((float* __restrict)&mpsar[0]);} 
                  constexpr int32_t mpsar_size() { return (NCOLS*NVAL);}
          };
          
          
          // The values describing the surface of dry part of the angled south roof wall
                               // ! (per each column) x number of columns
          template<int32_t ncols,int32_t nval>
          struct DPSARR1x_t {
                  
                 __ATTR_ALIGN__(16) __m128  dpsar[ncols*nval];
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr __m128 * __restrict dpsar_beg() { return (std::addressof(dpsar[0]));}
                 float * __restrict dpsar_fptr() { return ((float* __restrict)&dpsar[0]);} 
                 constexpr int32_t dpsar_size() { return (NCOLS*NVAL);}
          };
          
          
          // The values describing the surface of moist part of the angled east roof wall
                               // ! (per each column) x number of columns 
          template<int32_t ncols,int32_t nval>
          struct MPEARR1x_t {
                   
                 __ATTR_ALIGN__(16) __m128  mpear[ncols*nval];
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr __m128 * __restrict mpear_beg() { return (std::addressof(mpear[0]));}
                 float * __restrict mpear_fptr() { return ((float* __restrict)&mpear[0]);} 
                 constexpr int32_t mpear_size() { return (NCOLS*NVAL);}
          };
          
          
         // The values describing the surface of dry part of the angled east roof wall
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct DPEARR1x_t {
                 
                  __ATTR_ALIGN__(16) __m128 dpear[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict dpear_beg() { return (std::addressof(dpear[0]));}
                  float * __restrict dpear_fptr() { return ((float* __restrict)&dpear[0]);} 
                  constexpr int32_t dpear_size() { return (NCOLS*NVAL);}
          }; 
          
          
          // The values describing the surface of moist part of the angled west roof wall
                               // ! (per each column) x number of columns  
          template<int32_t ncols,int32_t nval>
          struct MPWARR1x_t {
                  
                   __ATTR_ALIGN__(16) __m128 mpwar[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict mpwar_beg() { return (std::addressof(mpwar[0]));}
                  float * __restrict mpwar_fptr() { return ((float* __restrict)&mpwar[0]);} 
                  constexpr int32_t mpwar_size() { return (NCOLS*NVAL);}
          }; 
          
          
           // The values describing the surface of dry part of the angled west roof wall
                               // ! (per each column) x number of columns  
          template<int32_t ncols,int32_t nval>
          struct DPWARR1x_t {
              
                   __ATTR_ALIGN__(16) __m128 dpwar[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict dpwar_beg() { return (std::addressof(dpwar[0]));}
                  float * __restrict dpwar_fptr() { return ((float* __restrict)&dpwar[0]);} 
                  constexpr int32_t dpwar_size() { return (NCOLS*NVAL);}
          }; 
          
          
           // The values describing the surface of moist part of the angled north roof wall
                               // ! (per each column) x number of columns  
          template<int32_t ncols,int32_t nval>
          struct MPNARR1x_t {
                
                   __ATTR_ALIGN__(16) __m128 mpnar[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr __m128 * __restrict mpnar_beg() { return (std::addressof(mpnar[0]));}
                  float * __restrict mpnar_fptr() { return ((float* __restrict)&mpnar[0]);} 
                  constexpr int32_t mpnar_size() { return (NCOLS*NVAL);}
          }; 
          
          
          // The values describing the surface of dry part of the angled north roof wall
                               // ! (per each column) x number of columns  
          template<int32_t ncols,int32_t nval>
          struct DPNARR1x_t {
               
                   __ATTR_ALIGN__(16) __m128 dpnar[ncols*nval];
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;
                   constexpr __m128 * __restrict dpnar_beg() { return (std::addressof(dpnar[0]));}
                   float * __restrict dpnar_fptr() { return ((float* __restrict)&dpnar[0]);} 
                   constexpr int32_t dpnar_size() { return (NCOLS*NVAL);}
          }; 
          
          
         // The values describing the complex permittivity of south walls
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct CESWC1x_t {
                  
                
                  __ATTR_ALIGN__(16) __m128 ceswr[ncols*nval];
 
                  __ATTR_ALIGN__(16) __m128 ceswi[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;  
                  constexpr __m128 * __restrict ceswr_beg() { return (std::addressof(ceswr[0]));}
                  float * __restrict ceswr_fptr() { return ((float* __restrict)&ceswr[0]);} 
                  constexpr int32_t ceswr_size() { return (NCOLS*NVAL);}
                  constexpr __m128 * __restrict ceswi_beg() { return (std::addressof(ceswi[0]));}
                  float * __restrict ceswi_fptr() { return ((float* __restrict)&ceswi[0]);} 
                  constexpr int32_t ceswi_size() { return (NCOLS*NVAL);}                
         };
         
         
         // The values describing the complex permeabillity of south walls
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct CMSWC1x_t {
           
                  __ATTR_ALIGN__(16) __m128 cmswr[ncols*nval];
 
                  __ATTR_ALIGN__(16) __m128 cmswi[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                  
                  constexpr __m128 * __restrict cmswr_beg() { return (std::addressof(cmswr[0]));}
                  float * __restrict cmswr_fptr() { return ((float* __restrict)&ceswr[0]);} 
                  constexpr int32_t cmswr_size() { return (NCOLS*NVAL);}
                  constexpr __m128 * __restrict  cmswi_beg() { return (std::addressof(cmswi[0]));}
                  float * __restrict cmswi_fptr() { return ((float* __restrict)&cmswi[0]);} 
                  constexpr int32_t cmswi_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values describing the complex permittivity of west walls
                               // ! (per each column) x number of columns 
         template<int32_t ncols,int32_t nval>
         struct CEWWC1x_t {
                 
                  __ATTR_ALIGN__(16) __m128 cewwr[ncols*nval];

                  __ATTR_ALIGN__(16) __m128 cewwi[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                  
                  constexpr __128 *  __restrict cewwr_beg() { return (std::addressof(cewwr[0]));}
                  float * __restrict cewwr_fptr() { return ((float* __restrict)&cewwr[0]);} 
                  constexpr int32_t cewwr_size() { return (NCOLS*NVAL);}
                  constexpr __m128 * __restrict cewwi_beg() { return (std::addressof(cewwi[0]));}
                  float * __restrict cewwi_fptr() { return ((float* __restrict)&cewwi[0]);} 
                  constexpr int32_t cewwi_size() { return (NCOLS*NVAL);}
         };
         
         
          // The values describing the complex permeability of west walls
                               // ! (per each column) x number of columns 
         template<int32_t ncols,int32_t nval>
         struct CMWWC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cmwwr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cmwwi[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                      
                  constexpr T * cmwwr_beg() { return (std::addressof(cmwwr[0]));}
                  constexpr int32_t cmwwr_size() { return (NCOLS*NVAL);}
                  constexpr T * cmwwi_beg() { return (std::addressof(cmwwi[0]));}
                  constexpr int32_t cmwwi_size() { return (NCOLS*NVAL);}
         };
         
         
          // The values describing the complex permittivity of east walls
                               // ! (per each column) x number of columns 
         template<int32_t ncols,int32_t nval>
         struct CEEWC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T ceewr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T ceewi[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                    
                  constexpr T * ceewr_beg() { return (std::addressof(ceewr[0]));}
                  constexpr int32_t ceewr_size() { return (NCOLS*NVAL);}
                  constexpr T * ceewi_beg() { return (std::addressof(ceewi[0]));}
                  constexpr int32_t ceewi_size() { return (NCOLS*NVAL);}
         };
         
         
          // The values describing the complex permeability of east walls
                               // ! (per each column) x number of columns 
         template<int32_t ncols,int32_t nval>
         struct CMEWC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cmewr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cmewi[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                  
                  constexpr T * cmewr_beg() { return (std::addressof(cmewr[0]));}
                  constexpr int32_t cmewr_size() { return (NCOLS*NVAL);}
                  constexpr T * cmewi_beg() { return (std::addressof(cmewi[0]));}
                  constexpr int32_t cmewi_size() { return (NCOLS*NVAL);}
         };
         
         
          // The values describing the complex permittivity of north walls
                               // ! (per each column) x number of columns 
         template<int32_t ncols,int32_t nval>
         struct CENWC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cenwr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cenwi[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                      
                  constexpr T * cenwr_beg() { return (std::addressof(cenwr[0]));}
                  constexpr int32_t cenwr_size() { return (NCOLS*NVAL);}
                  constexpr T * cenwi_beg() { return (std::addressof(cenwi[0]));}
                  constexpr int32_t cenwi_size() { return (NCOLS*NVAL);}
         };
         
         
          // The values describing the complex permeability of north walls
                               // ! (per each column) x number of columns 
         template<int32_t ncols,int32_t nval>
         struct CMNWC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cmnwr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cmnwi[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                    
                  constexpr T * cmnwr_beg() { return (std::addressof(cmnwr[0]));}
                  constexpr int32_t cmnarr_size() { return (NCOLS*NVAL);}
                  constexpr T * cmnwi_beg() { return (std::addressof(cmnwi[0]));}
                  constexpr int32_t cmnwi_size() { return (NCOLS*NVAL);}
         };
         
         
          // The values describing the complex permittivity of south angled roof
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct CESARC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cesarr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cesari[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                       
                  constexpr T * cesarr_beg() { return (std::addressof(cesarr[0]));}
                  constexpr int32_t cesarr_size() { return (NCOLS*NVAL);}
                  constexpr T * cesari_beg() { return (std::addressof(cesari[0]));}
                  constexpr int32_t cesari_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values describing the complex permeabillity of south angled roof
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct CMSARC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cmsarr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cmsari[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                    
                  constexpr T * cmsarr_beg() { return (std::addressof(cmsarr[0]));}
                  constexpr int32_t cmsarr_size() { return (NCOLS*NVAL);}
                  constexpr T * cmsari_beg() { return (std::addressof(cmsari[0]));}
                  constexpr int32_t cmsari_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values describing the complex permittivity of east angled roof
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct CEEARC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T ceearr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T ceeari[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                      
                  constexpr T * ceearr_beg() { return (std::addressof(ceearr[0]));}
                  constexpr int32_t ceearr_size() { return (NCOLS*NVAL);}
                  constexpr T * ceeari_beg() { return (std::addressof(ceeari[0]));}
                  constexpr int32_t ceeari_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values describing the complex permeabillity of east angled roof
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct CMEARC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cmearr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cmeari[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                  
                  constexpr T * cmearr_beg() { return (std::addressof(cmearr[0]));}
                  constexpr int32_t cmnarr_size() { return (NCOLS*NVAL);}
                  constexpr T * cmeari_beg() { return (std::addressof(cmeari[0]));}
                  constexpr int32_t cmeari_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values describing the complex permittivity of west angled roof
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct CEWARC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cewarr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cewari[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                     
                  constexpr T * cewarr_beg() { return (std::addressof(cewarr[0]));}
                  constexpr int32_t cewarr_size() { return (NCOLS*NVAL);}
                  constexpr T * cewari_beg() { return (std::addressof(cewari[0]));}
                  constexpr int32_t cewari_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values describing the complex permeabillity of west angled roof
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct CMWARC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cmwarr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cmwari[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                 
                  constexpr T * cmwarr_beg() { return (std::addressof(cmwarr[0]));}
                  constexpr int32_t cwnarr_size() { return (NCOLS*NVAL);}
                  constexpr T * cmwari_beg() { return (std::addressof(cmwari[0]));}
                  constexpr int32_t cmwari_size() { return (NCOLS*NVAL);}
         };
         
         
          // The values describing the complex permittivity of north angled roof
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct CENARC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cenarr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cenari[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                    
                  constexpr T * cenarr_beg() { return (std::addressof(cenarr[0]));}
                  constexpr int32_t cenarr_size() { return (NCOLS*NVAL);}
                  constexpr T * cenari_beg() { return (std::addressof(cenari[0]));}
                  constexpr int32_t cenari_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values describing the complex permeabillity of north angled roof
                               // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct CMNARC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cmnarr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cmnari[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                 
                  constexpr T * cmnarr_beg() { return (std::addressof(cmnarr[0]));}
                  constexpr int32_t cmnarr_size() { return (NCOLS*NVAL);}
                  constexpr T * cmnari_beg() { return (std::addressof(cmnari[0]));}
                  constexpr int32_t cmnari_size() { return (NCOLS*NVAL);}
         };
         
         
         // The components of south walls normal vector
                                    // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct NVSWR1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif            
                   T nvx[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T nvy[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T nvz[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr T * nvx_beg() { return (std::addressof(nvx[0]));}
                  constexpr int32_t nvx_size() { return (NCOLS*NVAL);}
                  constexpr T * nvy_beg() { return (std::addressof(nvy[0]));}
                  constexpr int32_t nvy_size() { return (NCOLS*NVAL);}
                  constexpr T * nvz_beg() { return (std::addressof(nvz[0]));}
                  constexpr int32_t nvz_size() { return (NCOLS*NVAL);}
         };
         
         
        // The components of east walls normal vector
                                    // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct NVEWR1x_t {
                  
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif            
                   T nvx[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T nvy[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T nvz[ncols*nval];   
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;              
                  constexpr T * nvx_beg() { return (std::addressof(nvx[0]));}
                  constexpr int32_t nvx_size() { return (NCOLS*NVAL);}
                  constexpr T * nvy_beg() { return (std::addressof(nvy[0]));}
                  constexpr int32_t nvy_size() { return (NCOLS*NVAL);}
                  constexpr T * nvz_beg() { return (std::addressof(nvz[0]));}
                  constexpr int32_t nvz_size() { return (NCOLS*NVAL);}
         };
         
         
        // The components of west walls normal vector
                                    // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct NVWWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif            
                   T nvx[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T nvy[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T nvz[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                   
                  constexpr T * nvx_beg() { return (std::addressof(nvx[0]));}
                  constexpr int32_t nvx_size() { return (NCOLS*NVAL);}
                  constexpr T * nvy_beg() { return (std::addressof(nvy[0]));}
                  constexpr int32_t nvy_size() { return (NCOLS*NVAL);}
                  constexpr T * nvz_beg() { return (std::addressof(nvz[0]));}
                  constexpr int32_t nvz_size() { return (NCOLS*NVAL);}
         };
         
         
        // The components of north walls normal vector
                                    // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct NVNWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif            
                   T nvx[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T nvy[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T nvz[ncols*nval];     
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                 
                  constexpr T * nvx_beg() { return (std::addressof(nvx[0]));}
                  constexpr int32_t nvx_size() { return (NCOLS*NVAL);}
                  constexpr T * nvy_beg() { return (std::addressof(nvy[0]));}
                  constexpr int32_t nvy_size() { return (NCOLS*NVAL);}
                  constexpr T * nvz_beg() { return (std::addressof(nvz[0]));}
                  constexpr int32_t nvz_size() { return (NCOLS*NVAL);}
         };
         
         
         // The components of each building normal vector
                                    // ! (per each column) x number of columns 
         template<int32_t ncols,int32_t nval>
         struct NVBR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif            
                   T nvx[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T nvy[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T nvz[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                    
                  constexpr T * nvx_beg() { return (std::addressof(nvx[0]));}
                  constexpr int32_t nvx_size() { return (NCOLS*NVAL);}
                  constexpr T * nvy_beg() { return (std::addressof(nvy[0]));}
                  constexpr int32_t nvy_size() { return (NCOLS*NVAL);}
                  constexpr T * nvz_beg() { return (std::addressof(nvz[0]));}
                  constexpr int32_t nvz_size() { return (NCOLS*NVAL);}
         };
         
         
         // The components of each building flat roof normal vector
                                    // ! (per each column) x number of columns 
         template<int32_t ncols,int32_t nval>
         struct NVFRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif            
                   T nvx[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T nvy[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T nvz[ncols*nval];   
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                         
                  constexpr T * nvx_beg() { return (std::addressof(nvx[0]));}
                  constexpr int32_t nvx_size() { return (NCOLS*NVAL);}
                  constexpr T * nvy_beg() { return (std::addressof(nvy[0]));}
                  constexpr int32_t nvy_size() { return (NCOLS*NVAL);}
                  constexpr T * nvz_beg() { return (std::addressof(nvz[0]));}
                  constexpr int32_t nvz_size() { return (NCOLS*NVAL);}
         };
         
         
          // The components of south angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct NVSARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif            
                   T nvx[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T nvy[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T nvz[ncols*nval];   
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                   
                  constexpr T * nvx_beg() { return (std::addressof(nvx[0]));}
                  constexpr int32_t nvx_size() { return (NCOLS*NVAL);}
                  constexpr T * nvy_beg() { return (std::addressof(nvy[0]));}
                  constexpr int32_t nvy_size() { return (NCOLS*NVAL);}
                  constexpr T * nvz_beg() { return (std::addressof(nvz[0]));}
                  constexpr int32_t nvz_size() { return (NCOLS*NVAL);}
         };
         
         
        // The components of east angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct NVEARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif            
                   T nvx[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T nvy[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T nvz[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                       
                  constexpr T * nvx_beg() { return (std::addressof(nvx[0]));}
                  constexpr int32_t nvx_size() { return (NCOLS*NVAL);}
                  constexpr T * nvy_beg() { return (std::addressof(nvy[0]));}
                  constexpr int32_t nvy_size() { return (NCOLS*NVAL);}
                  constexpr T * nvz_beg() { return (std::addressof(nvz[0]));}
                  constexpr int32_t nvz_size() { return (NCOLS*NVAL);}
         };
         
         
        // The components of west angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct NVWARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif            
                   T nvx[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T nvy[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T nvz[ncols*nval];  
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                         
                  constexpr T * nvx_beg() { return (std::addressof(nvx[0]));}
                  constexpr int32_t nvx_size() { return (NCOLS*NVAL);}
                  constexpr T * nvy_beg() { return (std::addressof(nvy[0]));}
                  constexpr int32_t nvy_size() { return (NCOLS*NVAL);}
                  constexpr T * nvz_beg() { return (std::addressof(nvz[0]));}
                  constexpr int32_t nvz_size() { return (NCOLS*NVAL);}
         };
         
         
        // The components of north angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct NVNARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif            
                   T nvx[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T nvy[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T nvz[ncols*nval]; 
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;                         
                  constexpr T * nvx_beg() { return (std::addressof(nvx[0]));}
                  constexpr int32_t nvx_size() { return (NCOLS*NVAL);}
                  constexpr T * nvy_beg() { return (std::addressof(nvy[0]));}
                  constexpr int32_t nvy_size() { return (NCOLS*NVAL);}
                  constexpr T * nvz_beg() { return (std::addressof(nvz[0]));}
                  constexpr int32_t nvz_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values of each south wall height and width
                        // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct HWSWR1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                T  hsw[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                T  wsw[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * hsw_beg() { return (std::addressof(hsw[0]));}
                constexpr int32_t hsw_size() { return (NCOLS*NVAL);}
                constexpr T * wsw_beg() { return (std::addressof(wsw[0]));}
                constexpr int32_t wsw_size() { return (NCOLS*NVAL);}
         };
         
         
          // The values of each east wall height and width
                        // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct HWEWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                T  hew[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                T  wew[ncols*nval];  
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;              
                constexpr T * hew_beg() { return (std::addressof(hew[0]));}
                constexpr int32_t hew_size() { return (NCOLS*NVAL);} 
                constexpr T * wew_beg() { return (std::addressof(wew[0]));}
                constexpr int32_t wew_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values of each west wall height and width
                        // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct HWWWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                T  hww[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                T  www[ncols*nval];      
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;          
                constexpr T * hww_beg() { return (std::addressof(hww[0]));}
                constexpr int32_t hww_size() { return (NCOLS*NVAL);}
                constexpr T * www_beg() { return (std::addressof(www[0]));}
                constexpr int32_t www_size() { return (NCOLS*NVAL);}
         };
         
         
          // The values of each north wall height and width
                        // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct HWNWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                T  hnw[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                T  wnw[ncols*nval];   
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;  
                constexpr T * hnw_beg() { return (std::addressof(hnw[0]));}
                constexpr int32_t hnw_size() { return (NCOLS*NVAL);}
                constexpr T * wnw_beg() { return (std::addressof(wnw[0]));}
                constexpr int32_t wnw_size() { return (NCOLS*NVAL);}                    
               
         };
         
         
         // The values of each flat roof height and width
                        // ! (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct HWFRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                T hfr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                T wfr[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * hfr_beg() { return (std::addressof(hfr[0]));}
                constexpr int32_t hfr_size() { return (NCOLS*NVAL);}
                constexpr T * wfr_beg() { return (std::addressof(wfr[0]));}
                constexpr int32_t wfr_size() { return (NCOLS*NVAL);}
         };
         
         
        // The values of each south non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         template<int32_t ncols,int32_t nval>
         struct HWSNFRR1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                
                T  hsnfr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif 
                T  wsnfr[ncols*nval];
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * hsnfr_beg() { return (std::addressof(hsnfr[0]));}
                constexpr int32_t hsnfr_size() { return (NCOLS*NVAL);}
                constexpr T * wsnfr_beg() { return (std::addressof(wsnfr[0]));}
                constexpr int32_t wsnfr_size() { return (NCOLS*NVAL);}
         };
         
         
        // The values of each east non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  [henfr,wenfr]
         template<int32_t ncols,int32_t nval>
         struct HWENFRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                
                T  henfr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif 
                T  wenfr[ncols*nval];
                
                 int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof                
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr T * henfr_beg() { return (std::addressof(henfr[0]));}
                 constexpr int32_t henfr_size() { return (NCOLS*NVAL);}
                 constexpr T * wenfr_beg() { return (std::addressof(wenfr[0]));}
                 constexpr int32_t wenfr_size() { return (NCOLS*NVAL);}
    };
         
         
        // The values of each west non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  [hwnfr,wwnfr]
         template<int32_t ncols,int32_t nval>
         struct HWWNFRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                
                T  hwnfr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif 
                T  wwnfr[ncols*nval];
                
                 int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof                          
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * hwnfr_beg() { return (std::addressof(hwnfr[0]));}
                constexpr int32_t hwnfr_size() { return (NCOLS*NVAL);}
                constexpr T * wwnfr_beg() { return (std::addressof(wwnfr[0]));}
                constexpr int32_t wwnfr_size() { return (NCOLS*NVAL);}
         };
         
         
        // The values of each north non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  [hnnfr,wnnfr]
         template<typename T>
         struct HWNNFRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                
                T  hnnfr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif 
                T  wnnfr[ncols*nval];
                
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof                        
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * hnnfr_beg() { return (std::addressof(hnnfr[0]));}
                constexpr int32_t hnnfr_size() { return (NCOLS*NVAL);}
                constexpr T * wnnfr_beg() { return (std::addressof(wnnfr[0]));}
                constexpr int32_t wnnfr_size() { return (NCOLS*NVAL);}
         };
         
         
         // Any kind of metallic structure fixed on the roof, e.g., an wire antenna, cylindrical object (ventillation)
         // or parabollic antenna, or similar type (per each column) x number of columns.
         template<int32_t ncols,int32_t nval>
         struct BRMSB1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                bool want[ncols*nval]; // a wire antennae
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                bool pant[ncols*nval]; // a parabollic antennae
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                bool yant[ncols*nval]; // yagi type antennae
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                bool lpda[ncols*nval]; // log-periodic dipole array
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                bool cant[ncols*nval]; // cell phone sector bars antennae
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                bool cylo[ncols*nval]; // any kind of cylindrical (ventillation) object
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * want_beg() { return (std::addressof(want[0]));}
                constexpr int32_t want_size() { return (NCOLS*NVAL);}
                constexpr T * pant_beg() { return (std::addressof(pant[0]));}
                constexpr int32_t pant_size() { return (NCOLS*NVAL);}
                constexpr T * yant_beg() { return (std::addressof(yant[0]));}
                constexpr int32_t yant_size() { return (NCOLS*NVAL);}
                constexpr T * lpda_beg() { return (std::addressof(lpda[0]));}
                constexpr int32_t lpda_size() { return (NCOLS*NVAL);}
                constexpr T * cant_beg() { return (std::addressof(cant[0]));}
                constexpr int32_t cant_size() { return (NCOLS*NVAL);}
                constexpr T * cylo_beg() { return (std::addressof(cylo[0]));}
                constexpr int32_t cylo_size() { return (NCOLS*NVAL);}
         };
         
         
         // The number of ventillation objects per single building roof
         // for every building column.
         template<int32_t ncols,int32_t nval>
         struct NVOI1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                int32_t nvo[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * nvo_beg() { return (std::addressof(nvo[0]));}
                constexpr int32_t nvo_size() { return (NCOLS*NVAL);}
         };
         
         
         // The number of wire antennae per single building roof
         // for every building column.
         template<int32_t ncols,int32_t nval>
         struct NWAI1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                int32_t  nwa[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * nwa_beg() { return (std::addressof(nwa[0]));}
                constexpr int32_t nwa_size() { return (NCOLS*NVAL);}
         };
         
         
          // The number of yagi-antennae per single building roof
         // for every building column.
         template<int32_t ncols,int32_t nval>
         struct NYAI1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                 
                int32_t  nya[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * nya_beg() { return (std::addressof(nya[0]));}
                constexpr int32_t nya_size() { return (NCOLS*NVAL);}
         };
         
         
          // The number of log-periodic dipole antennae per single building roof
         // for every building column.
         template<int32_t ncols,int32_t nval>
         struct NLPDAI1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                 
                int32_t  nlpda[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * nlpda_beg() { return (std::addressof(nlpda[0]));}
                constexpr int32_t nlpda_size() { return (NCOLS*NVAL);}
         };
         
         
          // The number of parabollic antennae per single building roof
         // for every building column.
         template<int32_t ncols,int32_t nval>
         struct NPAI1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                int32_t  npa[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * npa_beg() { return (std::addressof(npa[0]));}
                constexpr int32_t npa_size() { return (NCOLS*NVAL);}
         };
         
         
          // The number of cell-phone antennae per single building roof
         // for every building column.
         template<int32_t ncols,int32_t nval>
         struct NCPAI1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                int32_t  ncpa[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * ncpa_beg() { return (std::addressof(ncpa[0]));}
                constexpr int32_t ncpa_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values of RCS for the flat roof of
         // building column.
         template<int32_t ncols,int32_t nval>
         struct RCSFRR1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                
                 T  rcsfr[ncols*nval];
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr T * rcsfr_beg() { return (std::addressof(rcsfr[0]));}
                 constexpr int32_t rcsfr_size() { return (NCOLS*NVAL);}
         };
         
         
          // The values of RCS for the south wall of
         // of every building in the building column.
         template<int32_t ncols,int32_t nval>
         struct RCSSWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                T rcssw[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * rcssw_beg() { return (std::addressof(rcssw[0]));}
                constexpr int32_t rcssw_size() { return (NCOLS*NVAL);}
         };
         
         
          
          // The values of RCS for the east wall of
         // of every building in the building column.
         template<int32_t ncols,int32_t nval>
         struct RCSEWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                 T rcsew[ncols*nval];
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr T * rcsew_beg() { return (std::addressof(rcsew[0]));}
                 constexpr int32_t rcsew_size() { return (NCOLS*NVAL);}
         };
         
         
          
          // The values of RCS for the west wall of
         // of every building in the building column.
         template<int32_t ncols,int32_t nval>
         struct RCSWWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                 T rcsww[ncols*nval];
                 constexpr static int32_t NCOLS = ncols;
                 constexpr static int32_t NVAL  = nval;
                 constexpr T * rcsww_beg() { return (std::addressof(rcsww[0]));}
                 constexpr int32_t rcsww_size() { return (NCOLS*NVAL);}
         };
         
         
          
          // The values of RCS for the north wall of
         // of every building in the building column.
         template<int32_t ncols,int32_t nval>
         struct RCSNWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                T rcsnw[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * rcsnw_beg() { return (std::addressof(rcsnw[0]));}
                constexpr int32_t rcsnw_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values of RCS for the south angled roof of
         // of every building in the building column.
         template<int32_t ncols,int32_t nval>
         struct RCSSARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                   T rcssar[ncols*nval];
                   constexpr static int32_t NCOLS = ncols;
                   constexpr static int32_t NVAL  = nval;
                   constexpr T * rcssar_beg() { return (std::addressof(rcssar[0]));}
                   constexpr int32_t rcssar_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values of RCS for the south east roof of
         // of every building in the building column.
         template<int32_t ncols,int32_t nval>
         struct RCSEARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                T rcsear[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * rcsear_beg() { return (std::addressof(rcsear[0]));}
                constexpr int32_t rcsear_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values of RCS for the west angled roof of
         // of every building in the building column.
         template<int32_t ncols,int32_t nval>
         struct RCSWARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T rcswar[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr T * rcswar_beg() { return (std::addressof(rcswar[0]));}
                  constexpr int32_t rcswar_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values of RCS for the north angled roof of
         // of every building in the building column.
         template<int32_t ncols,int32_t nval>
         struct RCSNARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                T  rcsnar[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * rcsnar_beg() { return (std::addressof(rcsnar[0]));}
                constexpr int32_t rcsnar_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values of whole building surface area
         // of every building in building column
         template<int32_t ncols,int32_t nval>
         struct WBSAR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                T wbsa[ncols*nval];
                constexpr static int32_t NCOLS = ncols;
                constexpr static int32_t NVAL  = nval;
                constexpr T * wbsa_beg() { return (std::addressof(wbsa[0]));}
                constexpr int32_t wbsa_size() { return (NCOLS*NVAL);}
         };
         
         
         // The values of whole building internal volume
         // of every building in building column
         template<int32_t ncols,int32_t nval>
         struct WBIVR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                 T wbiv[ncols*nval];
                  constexpr static int32_t NCOLS = ncols;
                  constexpr static int32_t NVAL  = nval;
                  constexpr T * wbiv_beg() { return (std::addressof(wbiv[0]));}
                  constexpr int32_t wbiv_size() { return (NCOLS*NVAL);}
         };
         
         
         
         
           
     }// radiolocation



}// gms






























#endif /*__GMS_URBAN_ADT_STAT_H__*/
