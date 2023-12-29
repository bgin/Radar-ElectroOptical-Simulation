

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
           
           
           // Smoothing and approximating curve for linearly-piecewise height function (x-coordinate)
           // Smoothing and approximating curve for linearly-piecewise height function (y-coordinate)
           template<typename T>
           struct XYSMBHR1x_t {
                   
                  std::size_t nx;
                  std::size_t ny;
                  DC1D<T>_t   xsmbh;
                  DC1D<T>_t   ysmbh
           };
           
           
           // Empty space in-between of buildings (per single column) x number columns
           struct ESBBI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict esbb;
           };
           
           
           // An area values of in-between buildings empty spaces (per single column) x number columns
           template<typename T> 
           struct AESBBR1x_t {
                  
                  std::size_t ncols
                  std::size_t nval;
                  DC1D<T>_t   aesbb;
           };
           
           
           // An area values of each building (per single building column) x number columns
           template<typename T>
           struct ABCR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   abc;
           };
           
           
           // Number of south-facing walls (per each column) x number of columns
           struct SWPCI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict swpc;
           };
           
           
           // Number of east-facing walls (per each column) x number of columns
           struct EWPCI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict ewpc;
           };
           
           
           // Number of west-facing walls (per each column) x number of columns
           struct WWPCI1x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict wwpc;
           };
           
           
           // Number of north-facing walls (per each column) x number of columns
           struct NWPCI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict nwpc;
           };
           
           
           // Number of building roofs per each column x number of columns
           struct BRPCI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict brpc; 
           };
           
           
           //  An area of every building [flat] roof (per each column) x number of columns
           template<typename T>
           struct BRAPCR1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D<T>_t   brapc;
           };
           
           
           // Number of angled roof -- south facing roof wall (per each column) x number of columns
           struct SRWCI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict srwc;
           };
           
           
           // Number of angled roof -- east facing roof wall (per each column) x number of columns
           struct ERWCI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict erwc;
           };
           
           
           // Number of angled roof -- west facing roof wall (per each column)  x number of columns
           struct WRWCI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict wrwc;
           };
           
           
           // Number angled roof -- north facing roof wall (per each column) x number of columns
           struct NRWCI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict nrwc;
           };
           
           
           // An angled roof inclination (deg) -- south facing roof wall 
           // (per each column) x number of columns
           template<typename T>
           struct IDSRWR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   idsrw;
           };
           
           
           // An angled roof inclination (deg) -- east facing roof wall 
           // (per each column) x number of columns
           template<typename T>
           struct IDERWR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   iderw;
           };
           
           
           // An angled roof inclination (deg) -- west facing roof wall 
           // (per each column) x number of columns
           template<typename T>
           struct IDWRWR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   idwrw;
           };
           
           
           // An angled roof inclination (rad) -- south facing roof wall 
           // (per each column) x number of columns
           template<typename T>
           struct IDNRWR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   idnrw; 
           };
           
           
           //  An angled roof inclination (rad) -- south facing roof wall 
           //  (per each column) x number of columns
           template<typename T>
           struct IRSRWR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   irsrw; 
           };
           
           
           //  An angled roof inclination (rad) -- east facing roof wall 
           //  (per each column) x number of columns
           template<typename T>
           struct IRERWR1x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   irerw;
           };
           
           
           //  An angled roof inclination (rad) -- north facing roof wall 
           //  (per each column) x number of columns
           template<typename T>
           struct IRNRWR1x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   irnrw; 
           };
           
           
           // An angled roof inclination surface area -- south facing roof wall 
           // (per each column) x number of columns
           template<typename T>
           struct ISRAR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   isra;  
           };
           
           
           // An angled roof inclination surface area -- west facing roof wall 
           // (per each column) x number of columns
           template<typename T>
           struct IWRAR1x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   iwra;  
           };
           
           
           // An angled roof inclination surface area -- east facing roof wall 
           // (per each column) x number of columns
           template<typename T>
           struct IERAR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   iera;   
           };
           
           
           // An angled roof inclination surface area -- north facing roof wall 
           // (per each column) x number of columns
           template<typename T>
           struct INRAR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   inra;  
           };
           
           
           // South wall upper-facing edge inclination (rad) -- 
           // (per each column)  x number of columns
          template<typename T>
          struct SWUER1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   swue;   
          };
          
          
          // East wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<typename T>
          struct EWUER1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   ewue;   
          };
          
          
          // West wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<typename T>
          struct WWUER1x_t {
                  
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   wwue;   
          };
          
          
          // North wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          template<typename T>
          struct NWUER1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   nwue;    
          };
          
          
          // Shared right edge between the south wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          template<typename T>
          struct SEWER1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   sewe; 
          };
          
          // Shared left edge between the south wall and west wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          template<typename T>
          struct SWWER1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   swwe;
          };
          
          
           // Shared right edge between the north wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
           template<typename T>
           struct NWEER1x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   nwee;
           };
           
           
           // Shared right edge between the north wall and west wall inclination (rad) 
           // ! -- (per each column) x number of columns
           template<typename T>
           struct NWWER1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   nwwe; 
           };
           
           
           // Simple cell-based mesh
           template<typename T>
           struct CellMeshR1x_t {
                  
                  // Number of divisions along the x,y,z
                  int32_t ndiv[3];
                  // Number of cells
                  std::size_t nL;
                  // Compute numerical integration.
                  bool nint;
                  // Coordinates (x,y,z) of the center
                  // of Lth cell.
                  DC1D<T>_t  cx;
                  DC1D<T>_t  cy;
                  DC1D<T>_t  cz;
                  DC1D<T>_t  dv;
                  // (X,Y,Z) dimensions of the Ith
                  // rectangular volume cell (this is needed for
                  // the numerical integration)
                  DC1D<T>_t  dx;
                  DC1D<T>_t  dy;
                  DC1D<T>_t  dz;
           };
           
           
           // South walls surface area (for every building, per column) x number of columns
           template<typename T>
           struct SWSAR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   swsa;
           };
           
           
           // East walls surface area (for every building, per column) x number of columns
           template<typename T>
           struct EWSAR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   ewsa; 
           };
           
           
           // West walls surface area (for every building, per column) x number of columns
           template<typename T>
           struct WWSAR1x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   wwsa;  
           };
           
           
           // North walls surface area (for every building, per column) x number of columns
           template<typename T>
           struct NWSAR1x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D<T>_t   nwsa;   
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
         
         
           
     }// radiolocation



}






























#endif /*__GMS_URBAN_MODEL_H__*/
