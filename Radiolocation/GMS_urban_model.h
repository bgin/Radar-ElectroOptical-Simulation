

#ifndef __GMS_URBAN_MODEL_H__
#define __GMS_URBAN_MODEL_H__


namespace file_info {

     const unsigned int GMS_URBAN_MODEL_MAJOR = 1;
     const unsigned int GMS_URBAN_MODEL_MINOR = 0;
     const unsigned int GMS_URBAN_MODEL_MICRO = 0;
     const unsigned int GMS_URBAN_MODEL_FULLVER =
       1000U*GMS_URBAN_MODEL_MAJOR+100U*GMS_URBAN_MODEL_MINOR+
       10U*GMS_URBAN_MODEL_MICRO;
     const char * const GMS_URBAN_MODEL_CREATION_DATE = "23-12-2023 12:09 +00200 (SAT 23 DEC 2023 GMT+2)";
     const char * const GMS_URBAN_MODEL_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_URBAN_MODEL_SYNOPSIS      = "Abstract data types representing built-up area for diffraction and scattering calculations."

}

#include <cstdint>
#include "GMS_config"
#include "GMS_dyn_containers.hpp"

namespace gms {

       namespace radiolocation {
    
    
            // Building units per column
            // Building units per row
            struct BRCcount_t {
                   
                   int32_t nbpc;
                   int32_t nbpr;
                   int32_t * __restrict bpc;
                   int32_t * __restrict bpr;
            };
        
            // latitude   values (deg), per building
            // longtitude values (deg), per building
            template<typename T>
            struct BLatLondR1x_t {
                   
                   // Number of latitude   values (deg), per building
                   std::size_t nblatd;
                   // Number of longtitude values (deg), per building
                   std::size_t nblond;
                   DC1D<T>_t   blatd;
                   DC1D<T>_t   blond;
            };
            
            
            // latitude   values (rad), per building
            // longtitude values (rad), per building
            template<typename T>
            struct BLatLonrR1x_t {
                   
                   // Number of latitude   values (rad), per building
                   std::size_t nblatr;
                   // Number of longtitude values (rad), per building
                   std::size_t nblonr;
                   DC1D<T>_t   blatr;
                   DC1D<T>_t   blonr;
            };
            
            
            // ellipsoidal (radar waveform irradiating field) cells for building column
            template<typename T>
            struct EllpbR1x_t {
                    
                    // Number of ellipsoidal (radar waveform irradiating field) cells for building column
                    std::size_t nellpb;
                    DC1D<T>_t   ellpb;
            };
            
            
            // ! Parametric equation x (acos(t)) values (building)
            // ! 1st dimension building column, 2nd dimension 'x' parameter values
            // ! Parametric equation y (b(sin(t)) values (building)
            // 1st dimension building column, 2nd dimension 'y' parameter values
            template<typename T>
            struct PxybR1x_t {
                   
                   // number of building collumns
                   std::size_t nbpc;
                   std::size_t npxb;
                   std::size_t npyb
                   DC1D<T>_t   pxb;
                   DC1D<T>_t   pyb;
            };
            
            
            // Length, width and an area of every street
            template<typename T>
            struct SLWAR1x_t {
                    
                    // number of streets
                    std::size_t nstr;
                    DC1D<T>_t   lstr;
                    DC1D<T>_t   wstr;
                    DC1D<T>_t   astr;
            };
            
            
            //   ! Moisture of every street (2D array)
            //   ! 2nd dimension humidity values (per street), 1st dimension street numbers
            //    Percent of moist to dry area of evey street at each cell
            template<typename T>
            struct MStrR1x_t {
                    
                    // number of streets
                    std::size_t nstr;
                    std::size_t nmstr;
                    std::size_t npmstr;
                    DC1D<T>_t   mstr;
                    DC1D<T>_t   pmstr;
            };
            
            
            // !Coverage of every street (like: '1' for snow,'2' for mud, '3' for clay, ...etc)
            struct CStrR1x_t {
                    
                   // number of streets
                    std::size_t nstr; 
                    int32_t * __restrict cstr;
            };
            
            
            // Percent of covered to non-covered portion of every street at each irradiated cell
            // Average thickness of each layer (cover) of every street at each irradiated cell
            // Thickness of cover along street (number of values) at each irradiated cell
            template<typename T>
            struct CDStrR1x_t {
                    
                    // number of streets
                    std::size_t nstr;
                    std::size_t npcstr;
                    std::size_t natstr;
                    std::size_t ntcstr; 
                    DC1D<T>_t   pcstr;
                    DC1D<T>_t   atstr;
                    DC1D<T>_t   tcstr;
            };
            
            
            // Mu values for 'clean' street interpolated along the street length at each irradiated cell
            // Eps for 'clean' street street length interpolated at each irradiated cell
            template<typename T> {
            struct MEStr1C1x_t        
                     // number of streets
                    std::size_t nstr; 
                    std::size_t nmustr;
                    std::size_t nepstr;
                    DC2D<T>_t   mustr;
                    DC2D<T>_t   epstr;
            };
            
            
            // Mu for covered (i.e. by mud,snow,clay, ..etc) 
            // street interpolated along the street length at each irradiated cell
            // Eps for covered (i.e. by mud,snow,clay, ..etc) 
            // street  length interpolated at each irradiated cell
            template<typename T>
            struct MEStr2C1x_t {
                     
                      // number of streets
                    std::size_t nstr; 
                    std::size_t nmustr;
                    std::size_t nepstr;
                    DC2D<T>_t   mustr;
                    DC2D<T>_t   epstr; 
            };
            
            
            // Street curvature parametric equation u-parameter
            // Street curvature parametric equation v-parameter
            template<typename T>
            struct SCrvR1x_t {
                     
                     // number of streets
                    std::size_t nstr; 
                    std::size_t nustr;
                    std::size_t nvstr;
                    DC1D<T>_t   ustr;
                    DC1D<T>_t   vstr; 
            };
            
            
            //  Street surface normal vectors x-components 
            //  along the street length at each irradiated cell
            //  Street surface normal vectors y-components 
            //  along the street length at each irradiated cell
            //  Street surface normal vectors z-components 
            //  along the street length at each irradiated cell
            template<typename T>
            struct SNrmR1x_t {
                    
                     // number of streets
                    std::size_t nstr; 
                    std::size_t nx;
                    std::size_t ny;
                    std::size_t nz;
                    DC1D<T>_t   nvx;
                    DC1D<T>_t   nvy;
                    DC1D<T>_t   nvz;  
            };
            
            
            // latitude   values (deg), per street length (at irradiance point)
            // longtitude values (deg), per street length (at irradiance point)
            template<typename T>
            struct SIRCDR1x_t {
                     
                    // number of streets
                    std::size_t nstr; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D<T>_t   irlon;
                    DC1D<T>_t   irlat;  
            } 
            
            
           // latitude   values (rad), per street length (at irradiance point)
           // longtitude values (rad), per street length (at irradiance point)
           template<typename T>
           struct SIRCRR1x_t {
                    
                   // number of streets
                    std::size_t nstr; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D<T>_t   irlon;
                    DC1D<T>_t   irlat;  
           };
           
           
           // latitude   values (deg), of the building area (at irradiance point)
           // longtitude values (deg), of the building area (at irradiance point)
            template<typename T>
            struct BIRCDR1x_t {
                     
                    // number of buildings
                    std::size_t nbld; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D<T>_t   irlon;
                    DC1D<T>_t   irlat;  
            };
            
            
           // latitude   values (rad), of the building area (at irradiance point)
           // longtitude values (rad), of the building area (at irradiance point)
           template<typename T>
           struct BIRCRR1x_t {
                    
                   // number of buildings
                    std::size_t nbld; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D<T>_t   irlon;
                    DC1D<T>_t   irlat;  
           }; 
           
           
           // Urban area height map (at single building resolution)
           template<typename T>
           struct UHMapR1x_t {
                  
                  std::size_t nx;
                  std::size_t ny;
                  DC1D<T>_t   hmap; 
           };
           
           
           // Urban area height map (at single building resolution) -- 1st derivative
           template<typename T>
           struct UHDxDyR1x_t {
                  
                  std::size_t nx;
                  std::size_t ny;
                  DC1D<T>_t   hdxdy;
           };
           
           
           // Urban area height map (at single building resolution) -- gradient x-component
           // Urban area height map (at single building resolution) -- gradient y-component
           template<typename T>
           struct UHGradR1x_t {
                  
                  std::size_t nx;
                  std::size_t ny;
                  DC1D<T>_t   uhgx;
                  DC1D<T>_t   uhgy;
           };
     }



}






























#endif /*__GMS_URBAN_MODEL_H__*/
