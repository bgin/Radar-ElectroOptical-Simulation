

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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                   int32_t bpc[nbpc];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                   int32_t bpr[nbpr];

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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                     T blatd[nblatd];
 #if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                     T blond[nblond];

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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                   T blatr[nblatr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                   T blonr[nblonr];

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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                  T ellpb[nellpb];   
    
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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                    T pxb[nbpc*npxb];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                    T pyb[nbpc*npyb];

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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                     T lstr[nstr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T wstr[nstr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T astr[nstr];

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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                     T mstr[nstr*nmstr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T pmstr[nstr*npmstr];
                    
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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                 int32_t cstr[nstr];

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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                    T pcstr[nstr*npcstr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                    T atstr[nstr*natstr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                    T tcstr[nstr*ntcstr];
                 
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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                     T murstr[nstr*nmustr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T muistr[nstr*nmustr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T eprstr[nstr*nepstr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T epistr[nstr*nepstr];

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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                     T murstr[nstr*nmustr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T muistr[nstr*nmustr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T eprstr[nstr*nepstr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T epistr[nstr*nepstr];  

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
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                     T ustr[nstr*nustr];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T vstr[nstr*nvstr];

                 
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
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                     T nvx[nstr*nx];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T nvy[nstr*ny];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T nvz[nstr*nz];

                    
            };
            
            
            // latitude   values (deg), per street length (at irradiance point)
            // longtitude values (deg), per street length (at irradiance point)
           // number of streets
            template<typename T,int32_t nstr,
                     int32_t nlon,int32_t nlat>
            struct SIRCDR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                    T irlon[nstr*nlon];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                    T irlat[nstr*nlat]; 
 
                   
            };
            
            
           // latitude   values (rad), per street length (at irradiance point)
           // longtitude values (rad), per street length (at irradiance point)
           // number of streets
           template<typename T,int32_t nstr,
                    int32_t nlon,int32_t nlat>
           struct SIRCRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                    T irlon[nstr*nlon];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                    T irlat[nstr*nlat]; 

           };
           
           
           // latitude   values (deg), of the building area (at irradiance point)
           // longtitude values (deg), of the building area (at irradiance point)
            // number of buildings
            template<typename T,int32_t nbld,
                     int32_t nlon,int32_t nlat>
            struct BIRCDR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                     T irlon[nbld*nlon];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                     T irlat[nbld*nlat]; 
                   
            };
            
            
           // latitude   values (rad), of the building area (at irradiance point)
           // longtitude values (rad), of the building area (at irradiance point)
            // number of buildings
           template<typename T,int32_t nbld,
                    int32_t nlon,int32_t nlat>
           struct BIRCRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                    T irlon[nbld*nlon];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                    T irlat[nbld*nlat];

                 
           }; 
           
           
           // Urban area height map (at single building resolution)
           template<typename T,int32_t nx,int32_t ny>
           struct UHMapR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                  T hmap[nx*ny];

                 
           };
           
           
           // Urban area height map (at single building resolution) -- 1st derivative
           template<typename T,int32_t nx,int32_t ny>
           struct UHDxDyR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                 
                   T hdxdy[nx*ny];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif
                   T hdxdy[nx*ny];                 
           };
           
           
           // Urban area height map (at single building resolution) -- gradient x-component
           // Urban area height map (at single building resolution) -- gradient y-component
           template<typename T,int32_t nx,int32_t ny>
           struct UHGradR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T uhgx[nx];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T uhgy[ny];
           };
           
           
           // Smoothing and approximating curve for linearly-piecewise height function (x-coordinate)
           // Smoothing and approximating curve for linearly-piecewise height function (y-coordinate)
           template<typename T,int32_t nx,int32_t ny>
           struct XYSMBHR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                   T xsmbh[nx];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                   T ysmbh[ny];
           };
           
           
           // Empty space in-between of buildings (per single column) x number columns
           template<int32_t ncols,int32_t nval>
           struct ESBBI1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  int32_t esbb[ncols*nval];
           };
           
           
           // An area values of in-between buildings empty spaces (per single column) x number columns
           template<typename T,int32_t ncols,int32_t nval> 
           struct AESBBR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                   T aesbb[ncols*nval];
           };
           
           
           // An area values of each building (per single building column) x number columns
           template<typename T,int32_t ncols,int32_t nval>
           struct ABCR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                   T  abc[ncols*nvals];
           };
           
           
           // Number of south-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct SWPCI1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                   int32_t swpc[ncols];
           };
           
           
           // Number of east-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct EWPCI1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                 int32_t ewpc[ncols];
           };
           
           
           // Number of west-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct WWPCI1x_t {
 #if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                 int32_t  wwpc[ncols];
           };
           
           
           // Number of north-facing walls (per each column) x number of columns
           template<int32_t ncols>
           struct NWPCI1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                int32_t  nwpc[ncols];
           };
           
           
           // Number of building roofs per each column x number of columns
           template<int32_t ncols,int32_t nval>
           struct BRPCI1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                 int32_t brpc; 
           };
           
           
           //  An area of every building [flat] roof (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct BRAPCR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T brapc[ncols*nval];
           };
           
           
           // Number of angled roof -- south facing roof wall (per each column) x number of columns
           template<int32_t ncols>
           struct SRWCI1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  int32_t srwc[ncols];
           };
           
           
           // Number of angled roof -- east facing roof wall (per each column) x number of columns
           template<int32_t ncols>
           struct ERWCI1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  int32_t  erwc[ncols];
           };
           
           
           // Number of angled roof -- west facing roof wall (per each column)  x number of columns
           template<int32_t ncols>
           struct WRWCI1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                int32_t  wrwc[ncols];
           };
           
           
           // Number angled roof -- north facing roof wall (per each column) x number of columns
           template<int32_t ncols>
           struct NRWCI1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                int32_t nrwc[ncols];
           };
           
           
           // An angled roof inclination (deg) -- south facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IDSRWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                  T idsrw[ncols*nval];
           };
           
           
           // An angled roof inclination (deg) -- east facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IDERWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                  T iderw[ncols*nval];
           };
           
           
           // An angled roof inclination (deg) -- west facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IDWRWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                
                  T idwrw[ncols*nval];
           };
           
           
           // An angled roof inclination (rad) -- south facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IDNRWR1x_t {
 #if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                   T idnrw[ncols*nval]; 
           };
           
           
           //  An angled roof inclination (rad) -- south facing roof wall 
           //  (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IRSRWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                T irsrw[ncols*nval]; 
           };
           
           
           //  An angled roof inclination (rad) -- east facing roof wall 
           //  (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IRERWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                    T irerw[ncols*nval];
           };
           
           
           //  An angled roof inclination (rad) -- north facing roof wall 
           //  (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IRNRWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                    T irnrw[ncols*nval]; 
           };
           
           
           // An angled roof inclination surface area -- south facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct ISRAR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                  T isra[ncols*nval];  
           };
           
           
           // An angled roof inclination surface area -- west facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IWRAR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                   T iwra[ncols*nval];  
           };
           
           
           // An angled roof inclination surface area -- east facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct IERAR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                 
                  T iera[ncols*nval];   
           };
           
           
           // An angled roof inclination surface area -- north facing roof wall 
           // (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct INRAR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                   T inra[ncols*nval];   
           };
           
           
           // South wall upper-facing edge inclination (rad) -- 
           // (per each column)  x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct SWUER1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                  T swue[ncols*nval];    
          };
          
          
          // East wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct EWUER1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                 T ewue[ncols*nval];    
          };
          
          
          // West wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct WWUER1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T wwue[ncols*nval];     
          };
          
          
          // North wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct NWUER1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                  T nwue[ncols*nval];       
          };
          
          
          // Shared right edge between the south wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct SEWER1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                 T sewe[ncols*nval];     
          };
          
          // Shared left edge between the south wall and west wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct SWWER1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                  T swwe[ncols*nval];   
          };
          
          
           // Shared right edge between the north wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct NWEER1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                   T nwee[ncols*nval];   
           };
           
           
           // Shared right edge between the north wall and west wall inclination (rad) 
           // ! -- (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct NWWER1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                  T nwee[ncols*nval];   
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
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                  T swsa[ncols*nval];
           };
           
           
           // East walls surface area (for every building, per column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct EWSAR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                 T ewsa[ncols*nval];
           };
           
           
           // West walls surface area (for every building, per column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct WWSAR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                  T wwsa[ncols*nval];
           };
           
           
           // North walls surface area (for every building, per column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct NWSAR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                      
                  T wwsa[ncols*nval]; 
                  T * begin() { return (std::__adressof(wwsa[0]));}
                  constexpr int32_t size() { return (ncols*nval)};  
                 
           };
           
           
           // South walls moist/non moist logical (per column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct MNMSWB1x_t {
                    
                  
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                   bool  mnmsw[ncols*nval]; 
           };
           
           
            // East walls moist/non moist logical (per column) x number of columns
           template<int32_t ncols,int32_t nval>
           struct MNMEWB1x_t {
                                 
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                   bool  mnmew[ncols*nval]; 
           };
           
           
             // West walls moist/non moist logical (per column) x number of columns
           template<int32_t ncols,int32_t nval>  
           struct MNMWWB1x_t {
           
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                   bool  mnmww[ncols*nval]; 
           };
           
           
             // North walls moist/non moist logical (per column) x number of columns
           template<int32_t ncols,int32_t nval>  
           struct MNMNWB1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                   bool  mnmnw[ncols*nval]; 
           };
           
           
           // ! The values describing the ratio (percentage) of south wall 
                            // ! moisture to dryness (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct MDSWRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                   T  mdswr[ncols*nval];
           };
           
           
           // ! The values describing the ratio (percentage) of east wall 
                            // ! moisture to dryness (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct MDEWRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                   T   mdewr[ncols*nval];                 
                 
                  
           };
           
           
           //  The values describing the ratio (percentage) of west wall 
                             //! moisture to dryness (per each column) x number of columns
           template<typename T,int32_t ncols,int32_t nval>
           struct MDWWRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                  T   mdwwr[ncols*nval];
           };  
           
           
           //  The values describing the ratio (percentage) of north wall 
                             //! moisture to dryness (per each column) x number of columns              
           template<typename T,int32_t ncols,int32_t nval>
           struct MDNWRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                  T mdnwr[ncols*nval];
           };   
           
           
           // The logical values of flat roof moistness (being either moist or dry) 
                              // ! (per column) x number of columns    
          template<int32_t ncols,int32_t nval>            
          struct MDRB1x_t {
                 
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                 bool  mdr[ncols*nval];
          }; 
          
          
          // The values describing the ratio (percentage) of flat roof moisture to dryness 
                              // ! (per each column) x number of columns 
          template<typename T,int32_t ncols,int32_t nval>
          struct MDRRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                T   mdrr[ncols*nval]; 
          };
          
          
          // The values describing the surface of moist part of the flat roof 
                               // ! (per each column) x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct MPFRR1x_t {
                  
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                      
                 T   mpfr[ncols*nval]; 
          };     
          
          
          // The values describing the surface of dry part of the flat roof 
                               // ! (per each column) x number of columns  
          template<typename T,int32_t ncols,int32_t nval>
          struct DPFRR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                      
                  T   dpfr[ncols*nval]; 
          };   
          
          
          //  The values describing the surface of moist part of the south wall 
                                // ! (per each column) x number of columns   
          template<typename T,int32_t ncols,int32_t nval>
          struct MPSWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                 T  mpsw[ncols*nval]; 
          };  
          
          
          //  The values describing the surface of dry part of the south wall 
                               //  ! (per each column) x number of columns  
          template<typename T,int32_t ncols,int32_t nval>
          struct DPSWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                  T  dpsw[ncols*nval]; 
          };   
          
          
          //  The values describing the surface of moist part of the east wall 
                                // ! (per each column) x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct MPEWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                  T mpew[ncols*nval]; 
          };   
          
         // The values describing the surface of dry part of the east wall 
                                // ! (per each column) x number of columns 
          template<typename T,int32_t ncols,int32_t nval>
          struct DPEWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                  T  dpew[ncols*nval]; 
          };  
          
          
         // The values describing the surface of moist part of the west wall 
                                 //! (per each column) x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct MPWWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                  T mpww[ncols*nval]; 
          }; 
          
          
        //  The values describing the surface of dry part of the west wall 
                                 //! (per each column) x number of columns 
          template<typename T,int32_t ncols,int32_t nval>
          struct DPWWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                  T dpww[ncols*nval]; 
          };  
          
          
        // The values describing the surface of moist part of the north wall 
                                 //! (per each column) x number of columns
         template<typename T,int32_t ncols,int32_t nval>
         struct MPNWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T  mpnw[ncols*nval]; 
          }; 
          
          
         //  The values describing the surface of dry part of the north wall 
                                // ! (per each column) x number of columns
         template<typename T,int32_t ncols,int32_t nval>
         struct DPNWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                 T dpnw[ncols*nval]; 
          }; 
          
          
          // The values describing the surface of moist part of the angled south roof wall
                               // ! (per each column) x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct MPSARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                  T mpsar[ncols*nval];
          };
          
          
          // The values describing the surface of dry part of the angled south roof wall
                               // ! (per each column) x number of columns
          template<typename T,int32_t ncols,int32_t nval>
          struct DPSARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                 T dpsar[ncols*nval];
          };
          
          
          // The values describing the surface of moist part of the angled east roof wall
                               // ! (per each column) x number of columns 
          template<typename T,int32_t ncols,int32_t nval>
          struct MPEARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                 T  mpear[ncols*nval];
          };
          
          
         // The values describing the surface of dry part of the angled east roof wall
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
         struct DPEARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                     
                  T dpear[ncols*nval];
          }; 
          
          
          // The values describing the surface of moist part of the angled west roof wall
                               // ! (per each column) x number of columns  
          template<typename T,int32_t ncols,int32_t nval>
          struct MPWARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                  T mpwar[ncols*nval];
          }; 
          
          
           // The values describing the surface of dry part of the angled west roof wall
                               // ! (per each column) x number of columns  
          template<typename T,int32_t ncols,int32_t nval>
          struct DPWARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                  T dpwar[ncols*nval];
          }; 
          
          
           // The values describing the surface of moist part of the angled north roof wall
                               // ! (per each column) x number of columns  
          template<typename T,int32_t ncols,int32_t nval>
          struct MPNARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                  T mpnar[ncols*nval];
          }; 
          
          
          // The values describing the surface of dry part of the angled north roof wall
                               // ! (per each column) x number of columns  
          template<typename T,int32_t ncols,int32_t nval>
          struct DPNARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                    
                   T dpnar[ncols*nval];
          }; 
          
          
         // The values describing the complex permittivity of south walls
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
         struct CESWC1x_t {
                  
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T ceswr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T ceswi[ncols*nval];                  
         };
         
         
         // The values describing the complex permeabillity of south walls
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
         struct CMSWC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cmswr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cmswi[ncols*nval];                    
                  
         };
         
         
         // The values describing the complex permittivity of west walls
                               // ! (per each column) x number of columns 
         template<typename T,int32_t ncols,int32_t nval>
         struct CEWWC1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                       
                  T cewwr[ncols*nval];
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif     
                  T cewwi[ncols*nval];                   
                  
         };
         
         
          // The values describing the complex permeability of west walls
                               // ! (per each column) x number of columns 
         template<typename T,int32_t ncols,int32_t nval>
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
                  
         };
         
         
          // The values describing the complex permittivity of east walls
                               // ! (per each column) x number of columns 
         template<typename T,int32_t ncols,int32_t nval>
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
                  
         };
         
         
          // The values describing the complex permeability of east walls
                               // ! (per each column) x number of columns 
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
          // The values describing the complex permittivity of north walls
                               // ! (per each column) x number of columns 
         template<typename T,int32_t ncols,int32_t nval>
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
                  
         };
         
         
          // The values describing the complex permeability of north walls
                               // ! (per each column) x number of columns 
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
          // The values describing the complex permittivity of south angled roof
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                  
         };
         
         
         // The values describing the complex permeabillity of south angled roof
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                
         };
         
         
         // The values describing the complex permittivity of east angled roof
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
         // The values describing the complex permeabillity of east angled roof
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                  
         };
         
         
         // The values describing the complex permittivity of west angled roof
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
         // The values describing the complex permeabillity of west angled roof
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
          // The values describing the complex permittivity of north angled roof
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
         // The values describing the complex permeabillity of north angled roof
                               // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
         // The components of south walls normal vector
                                    // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
         };
         
         
        // The components of east walls normal vector
                                    // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
        // The components of west walls normal vector
                                    // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
        // The components of north walls normal vector
                                    // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                  
         };
         
         
         // The components of each building normal vector
                                    // ! (per each column) x number of columns 
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
         // The components of each building flat roof normal vector
                                    // ! (per each column) x number of columns 
         template<typename T,int32_t ncols,int32_t nval>
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
                
         };
         
         
          // The components of south angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
        // The components of east angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
        // The components of west angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                 
         };
         
         
        // The components of north angled roof normal vector
                                    // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
                  
         };
         
         
         // The values of each south wall height and width
                        // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
         };
         
         
          // The values of each east wall height and width
                        // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
               
         };
         
         
         // The values of each west wall height and width
                        // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
              
         };
         
         
          // The values of each north wall height and width
                        // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
               
         };
         
         
         // The values of each flat roof height and width
                        // ! (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
         };
         
         
        // The values of each south non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         template<typename T,int32_t ncols,int32_t nval>
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
         };
         
         
        // The values of each east non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  [henfr,wenfr]
         template<typename T,int32_t ncols,int32_t nval>
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
         
    };
         
         
        // The values of each west non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  [hwnfr,wwnfr]
         template<typename T,int32_t ncols,int32_t nval>
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
         };
         
         
         // The values of RCS for the flat roof of
         // building column.
         template<typename T,int32_t ncols,int32_t nval>
         struct RCSFRR1x_t {
                
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                
                 T  rcsfr[ncols*nval];
         };
         
         
          // The values of RCS for the south wall of
         // of every building in the building column.
         template<typename T,int32_t ncols,int32_t nval>
         struct RCSSWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                T rcssw[ncols*nval];
         };
         
         
          
          // The values of RCS for the east wall of
         // of every building in the building column.
         template<typename T,int32_t ncols,int32_t nval>
         struct RCSEWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                 T rcsew[ncols*nval];
         };
         
         
          
          // The values of RCS for the west wall of
         // of every building in the building column.
         template<typename T,int32_t ncols,int32_t nval>
         struct RCSWWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                 T rcsww[ncols*nval];
         };
         
         
          
          // The values of RCS for the north wall of
         // of every building in the building column.
         template<typename T,int32_t ncols,int32_t nval>
         struct RCSNWR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                T rcsnw[ncols*nval];
         };
         
         
         // The values of RCS for the south angled roof of
         // of every building in the building column.
         template<typename T,int32_t ncols,int32_t nval>
         struct RCSSARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                   T rcssar[ncols*nval];
         };
         
         
         // The values of RCS for the south east roof of
         // of every building in the building column.
         template<typename T,int32_t ncols,int32_t nval>
         struct RCSEARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                T rcsear[ncols*nval];
         };
         
         
         // The values of RCS for the west angled roof of
         // of every building in the building column.
         template<typename T,int32_t ncols,int32_t nval>
         struct RCSWARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                  T rcswar[ncols*nval];
         };
         
         
         // The values of RCS for the north angled roof of
         // of every building in the building column.
         template<typename T,int32_t ncols,int32_t nval>
         struct RCSNARR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                  
                T  rcsnar[ncols*nval];
         };
         
         
         // The values of whole building surface area
         // of every building in building column
         template<typename T,int32_t ncols,int32_t nval>
         struct WBSAR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                T wbsa[ncols*nval];
         };
         
         
         // The values of whole building internal volume
         // of every building in building column
         template<typename T,int32_t ncols,int32_t nval>
         struct WBIVR1x_t {
#if defined (__AVX512F__)                   
                   __ATTR_ALIGN__(64)                    
#elif defined (__AVX__) || defined (__AVX2__)                
                   __ATTR_ALIGN__(32) 
#elif defined (__SSE__)
                   __ATTR_ALIGN__(16)                  
#endif                   
                 T wbiv[ncols*nval];
         };
         
         
         
         
           
     }// radiolocation



}// gms






























#endif /*__GMS_URBAN_ADT_STAT_H__*/
