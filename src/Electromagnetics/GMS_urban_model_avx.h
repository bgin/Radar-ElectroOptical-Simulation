

#ifndef __GMS_URBAN_MODEL_AVX_H__
#define __GMS_URBAN_MODEL_AVX_H__ 311220230816


namespace file_info {

     const unsigned int GMS_URBAN_MODEL_AVX_MAJOR = 1;
     const unsigned int GMS_URBAN_MODEL_AVX_MINOR = 0;
     const unsigned int GMS_URBAN_MODEL_AVX_MICRO = 0;
     const unsigned int GMS_URBAN_MODEL_AVX_FULLVER =
       1000U*GMS_URBAN_MODEL_AVX_MAJOR+100U*GMS_URBAN_MODEL_AVX_MINOR+
       10U*GMS_URBAN_MODEL_AVX_MICRO;
     const char * const GMS_URBAN_MODEL_AVX_CREATION_DATE = "31-12-2023 08:21AM +00200 (SUN 31 DEC 2023 GMT+2)";
     const char * const GMS_URBAN_MODEL_AVX_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_URBAN_MODEL_AVX_SYNOPSIS      = "Abstract data (AVX-based) types representing built-up area for diffraction and scattering calculations."

}

#include <cstdint>
#include "GMS_config"
#include "GMS_dyn_containers_avx.hpp"

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
            
            struct BLatLondR8x_t {
                   
                   // Number of latitude   values (deg), per building
                   std::size_t nblatd;
                   // Number of longtitude values (deg), per building
                   std::size_t nblond;
                   DC1D_m256_t  blatd;
                   DC1D_m256_t  blond;
            };
            
            
            // latitude   values (rad), per building
            // longtitude values (rad), per building
            
            struct BLatLonrR8x_t {
                   
                   // Number of latitude   values (rad), per building
                   std::size_t nblatr;
                   // Number of longtitude values (rad), per building
                   std::size_t nblonr;
                   DC1D_m256_t  blatr;
                   DC1D_m256_t  blonr;
            };
            
            
            // ellipsoidal (radar waveform irradiating field) cells for building column
            
            struct EllpbR8x_t {
                    
                    // Number of ellipsoidal (radar waveform irradiating field) cells for building column
                    std::size_t nellpb;
                    DC1D_m256_t  ellpb;
            };
            
            
            // ! Parametric equation x (acos(t)) values (building)
            // ! 1st dimension building column, 2nd dimension 'x' parameter values
            // ! Parametric equation y (b(sin(t)) values (building)
            // 1st dimension building column, 2nd dimension 'y' parameter values
            
            struct PxybR8x_t {
                   
                   // number of building collumns
                   std::size_t nbpc;
                   std::size_t npxb;
                   std::size_t npyb
                   DC1D_m256_t  pxb;
                   DC1D_m256_t  pyb;
            };
            
            
            // Length, width and an area of every street
            
            struct SLWAR8x_t {
                    
                    // number of streets
                    std::size_t nstr;
                    DC1D_m256_t  lstr;
                    DC1D_m256_t  wstr;
                    DC1D_m256_t  astr;
            };
            
            
            //   ! Moisture of every street (2D array)
            //   ! 2nd dimension humidity values (per street), 1st dimension street numbers
            //    Percent of moist to dry area of evey street at each cell
            
            struct MStrR8x_t {
                    
                    // number of streets
                    std::size_t nstr;
                    std::size_t nmstr;
                    std::size_t npmstr;
                    DC1D_m256_t  mstr;
                    DC1D_m256_t  pmstr;
            };
            
            
            // !Coverage of every street (like: '1' for snow,'2' for mud, '3' for clay, ...etc)
            struct CStrR8x_t {
                    
                   // number of streets
                    std::size_t nstr; 
                    int32_t * __restrict cstr;
            };
            
            
            // Percent of covered to non-covered portion of every street at each irradiated cell
            // Average thickness of each layer (cover) of every street at each irradiated cell
            // Thickness of cover along street (number of values) at each irradiated cell
            
            struct CDStrR8x_t {
                    
                    // number of streets
                    std::size_t nstr;
                    std::size_t npcstr;
                    std::size_t natstr;
                    std::size_t ntcstr; 
                    DC1D_m256_t  pcstr;
                    DC1D_m256_t  atstr;
                    DC1D_m256_t  tcstr;
            };
            
            
            // Mu values for 'clean' street interpolated along the street length at each irradiated cell
            // Eps for 'clean' street street length interpolated at each irradiated cell
            
            struct MEStr1C8x_t   {     
                     // number of streets
                    std::size_t nstr; 
                    std::size_t nmustr;
                    std::size_t nepstr;
                    DC1D_ymm8c4_t   mustr;
                    DC1D_ymm8c4_t   epstr;
            };
            
            
            // Mu for covered (i.e. by mud,snow,clay, ..etc) 
            // street interpolated along the street length at each irradiated cell
            // Eps for covered (i.e. by mud,snow,clay, ..etc) 
            // street  length interpolated at each irradiated cell
            
            struct MEStr2C8x_t {
                     
                      // number of streets
                    std::size_t nstr; 
                    std::size_t nmustr;
                    std::size_t nepstr;
                    DC1D_ymm8c4_t   mustr;
                    DC1D_ymm8c4_t   epstr; 
            };
            
            
            // Street curvature parametric equation u-parameter
            // Street curvature parametric equation v-parameter
            
            struct SCrvR8x_t {
                     
                     // number of streets
                    std::size_t nstr; 
                    std::size_t nustr;
                    std::size_t nvstr;
                    DC1D_m256_t  ustr;
                    DC1D_m256_t  vstr; 
            };
            
            
            //  Street surface normal vectors x-components 
            //  along the street length at each irradiated cell
            //  Street surface normal vectors y-components 
            //  along the street length at each irradiated cell
            //  Street surface normal vectors z-components 
            //  along the street length at each irradiated cell
            
            struct SNrmR8x_t {
                    
                     // number of streets
                    std::size_t nstr; 
                    std::size_t nx;
                    std::size_t ny;
                    std::size_t nz;
                    DC1D_m256_t  nvx;
                    DC1D_m256_t  nvy;
                    DC1D_m256_t  nvz;  
            };
            
            
            // latitude   values (deg), per street length (at irradiance point)
            // longtitude values (deg), per street length (at irradiance point)
            
            struct SIRCDR8x_t {
                     
                    // number of streets
                    std::size_t nstr; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D_m256_t  irlon;
                    DC1D_m256_t  irlat;  
            } 
            
            
           // latitude   values (rad), per street length (at irradiance point)
           // longtitude values (rad), per street length (at irradiance point)
           
           struct SIRCRR8x_t {
                    
                   // number of streets
                    std::size_t nstr; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D_m256_t  irlon;
                    DC1D_m256_t  irlat;  
           };
           
           
           // latitude   values (deg), of the building area (at irradiance point)
           // longtitude values (deg), of the building area (at irradiance point)
            
            struct BIRCDR8x_t {
                     
                    // number of buildings
                    std::size_t nbld; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D_m256_t  irlon;
                    DC1D_m256_t  irlat;  
            };
            
            
           // latitude   values (rad), of the building area (at irradiance point)
           // longtitude values (rad), of the building area (at irradiance point)
           
           struct BIRCRR8x_t {
                    
                   // number of buildings
                    std::size_t nbld; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D_m256_t  irlon;
                    DC1D_m256_t  irlat;  
           }; 
           
           
           // Urban area height map (at single building resolution)
           
           struct UHMapR8x_t {
                  
                  std::size_t nx;
                  std::size_t ny;
                  DC1D_m256_t  hmap; 
           };
           
           
           // Urban area height map (at single building resolution) -- 1st derivative
           
           struct UHDxDyR8x_t {
                  
                  std::size_t nx;
                  std::size_t ny;
                  DC1D_m256_t  hdxdy;
           };
           
           
           // Urban area height map (at single building resolution) -- gradient x-component
           // Urban area height map (at single building resolution) -- gradient y-component
           
           struct UHGradR8x_t {
                  
                  std::size_t nx;
                  std::size_t ny;
                  DC1D_m256_t  uhgx;
                  DC1D_m256_t  uhgy;
           };
           
           
           // Smoothing and approximating curve for linearly-piecewise height function (x-coordinate)
           // Smoothing and approximating curve for linearly-piecewise height function (y-coordinate)
           
           struct XYSMBHR8x_t {
                   
                  std::size_t nx;
                  std::size_t ny;
                  DC1D_m256_t  xsmbh;
                  DC1D_m256_t  ysmbh
           };
           
           
           // Empty space in-between of buildings (per single column) x number columns
           struct ESBBI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict esbb;
           };
           
           
           // An area values of in-between buildings empty spaces (per single column) x number columns
            
           struct AESBBR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval;
                  DC1D_m256_t  aesbb;
           };
           
           
           // An area values of each building (per single building column) x number columns
           
           struct ABCR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  abc;
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
           
           struct BRAPCR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  brapc;
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
           
           struct IDSRWR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  idsrw;
           };
           
           
           // An angled roof inclination (deg) -- east facing roof wall 
           // (per each column) x number of columns
           
           struct IDERWR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  iderw;
           };
           
           
           // An angled roof inclination (deg) -- west facing roof wall 
           // (per each column) x number of columns
           
           struct IDWRWR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  idwrw;
           };
           
           
           // An angled roof inclination (rad) -- south facing roof wall 
           // (per each column) x number of columns
           
           struct IDNRWR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  idnrw; 
           };
           
           
           //  An angled roof inclination (rad) -- south facing roof wall 
           //  (per each column) x number of columns
           
           struct IRSRWR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  irsrw; 
           };
           
           
           //  An angled roof inclination (rad) -- east facing roof wall 
           //  (per each column) x number of columns
           
           struct IRERWR8x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  irerw;
           };
           
           
           //  An angled roof inclination (rad) -- north facing roof wall 
           //  (per each column) x number of columns
           
           struct IRNRWR8x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  irnrw; 
           };
           
           
           // An angled roof inclination surface area -- south facing roof wall 
           // (per each column) x number of columns
           
           struct ISRAR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  isra;  
           };
           
           
           // An angled roof inclination surface area -- west facing roof wall 
           // (per each column) x number of columns
           
           struct IWRAR8x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  iwra;  
           };
           
           
           // An angled roof inclination surface area -- east facing roof wall 
           // (per each column) x number of columns
           
           struct IERAR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  iera;   
           };
           
           
           // An angled roof inclination surface area -- north facing roof wall 
           // (per each column) x number of columns
           
           struct INRAR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  inra;  
           };
           
           
           // South wall upper-facing edge inclination (rad) -- 
           // (per each column)  x number of columns
          
          struct SWUER8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  swue;   
          };
          
          
          // East wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          
          struct EWUER8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  ewue;   
          };
          
          
          // West wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          
          struct WWUER8x_t {
                  
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  wwue;   
          };
          
          
          // North wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          
          struct NWUER8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  nwue;    
          };
          
          
          // Shared right edge between the south wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          
          struct SEWER8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  sewe; 
          };
          
          // Shared left edge between the south wall and west wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          
          struct SWWER8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  swwe;
          };
          
          
           // Shared right edge between the north wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
           
           struct NWEER8x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  nwee;
           };
           
           
           // Shared right edge between the north wall and west wall inclination (rad) 
           // ! -- (per each column) x number of columns
           
           struct NWWER8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  nwwe; 
           };
           
           
           // Simple cell-based mesh
           
           struct CellMeshR8x_t {
                  
                  // Number of divisions along the x,y,z
                  int32_t ndiv[3];
                  // Number of cells
                  std::size_t nL;
                  // Compute numerical integration.
                  bool nint;
                  // Coordinates (x,y,z) of the center
                  // of Lth cell.
                  DC1D_m256_t cx;
                  DC1D_m256_t cy;
                  DC1D_m256_t cz;
                  DC1D_m256_t dv;
                  // (X,Y,Z) dimensions of the Ith
                  // rectangular volume cell (this is needed for
                  // the numerical integration)
                  DC1D_m256_t dx;
                  DC1D_m256_t dy;
                  DC1D_m256_t dz;
           };
           
           
           // South walls surface area (for every building, per column) x number of columns
           
           struct SWSAR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  swsa;
           };
           
           
           // East walls surface area (for every building, per column) x number of columns
           
           struct EWSAR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  ewsa; 
           };
           
           
           // West walls surface area (for every building, per column) x number of columns
           
           struct WWSAR8x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  wwsa;  
           };
           
           
           // North walls surface area (for every building, per column) x number of columns
           
           struct NWSAR8x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m256_t  nwsa;   
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
           
           struct MDSWRR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mdswr;
           };
           
           
           // ! The values describing the ratio (percentage) of east wall 
                            // ! moisture to dryness (per each column) x number of columns
           
           struct MDEWRR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mdewr;
           };
           
           
           //  The values describing the ratio (percentage) of west wall 
                             //! moisture to dryness (per each column) x number of columns
           
           struct MDWWRR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mdwwr;
           };  
           
           
           //  The values describing the ratio (percentage) of north wall 
                             //! moisture to dryness (per each column) x number of columns              
           
           struct MDNWRR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mdnwr;
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
          
          struct MDRRR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mdrr; 
          };
          
          
          // The values describing the surface of moist part of the flat roof 
                               // ! (per each column) x number of columns
          
          struct MPFRR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mpfr; 
          };     
          
          
          // The values describing the surface of dry part of the flat roof 
                               // ! (per each column) x number of columns  
          
          struct DPFRR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  dpfr; 
          };   
          
          
          //  The values describing the surface of moist part of the south wall 
                                // ! (per each column) x number of columns   
          
          struct MPSWR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mpsw; 
          };  
          
          
          //  The values describing the surface of dry part of the south wall 
                               //  ! (per each column) x number of columns  
          
          struct DPSWR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  dpsw; 
          };   
          
          
          //  The values describing the surface of moist part of the east wall 
                                // ! (per each column) x number of columns
          
          struct MPEWR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mpew; 
          };   
          
         // The values describing the surface of dry part of the east wall 
                                // ! (per each column) x number of columns 
          
          struct DPEWR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  dpew; 
          };  
          
          
         // The values describing the surface of moist part of the west wall 
                                 //! (per each column) x number of columns
          
          struct MPWWR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mpww; 
          }; 
          
          
        //  The values describing the surface of dry part of the west wall 
                                 //! (per each column) x number of columns 
          
          struct DPWWR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  dpww; 
          };  
          
          
        // The values describing the surface of moist part of the north wall 
                                 //! (per each column) x number of columns
         
         struct MPNWR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mpnw; 
          }; 
          
          
         //  The values describing the surface of dry part of the north wall 
                                // ! (per each column) x number of columns
         
         struct DPNWR8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  dpnw; 
          }; 
          
          
          // The values describing the surface of moist part of the angled south roof wall
                               // ! (per each column) x number of columns
          
          struct MPSARR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mpsar;
          };
          
          
          // The values describing the surface of dry part of the angled south roof wall
                               // ! (per each column) x number of columns
          
          struct DPSARR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  dpsar;
          };
          
          
          // The values describing the surface of moist part of the angled east roof wall
                               // ! (per each column) x number of columns 
          
          struct MPEARR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mpear;
          };
          
          
         // The values describing the surface of dry part of the angled east roof wall
                               // ! (per each column) x number of columns  
         
         struct DPEARR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  dpear;
          }; 
          
          
          // The values describing the surface of moist part of the angled west roof wall
                               // ! (per each column) x number of columns  
          
          struct MPWARR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mpwar;
          }; 
          
          
           // The values describing the surface of dry part of the angled west roof wall
                               // ! (per each column) x number of columns  
          
          struct DPWARR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  dpwar;
          }; 
          
          
           // The values describing the surface of moist part of the angled north roof wall
                               // ! (per each column) x number of columns  
          
          struct MPNARR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  mpnar;
          }; 
          
          
          // The values describing the surface of dry part of the angled north roof wall
                               // ! (per each column) x number of columns  
          
          struct DPNARR8x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m256_t  dpnar;
          }; 
          
          
         // The values describing the complex permittivity of south walls
                               // ! (per each column) x number of columns  
         
         struct CESWC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cesw;
         };
         
         
         // The values describing the complex permeabillity of south walls
                               // ! (per each column) x number of columns  
         
         struct CMSWC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cmsw;
         };
         
         
         // The values describing the complex permittivity of west walls
                               // ! (per each column) x number of columns 
         
         struct CEWWC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   ceww;
         };
         
         
          // The values describing the complex permeability of west walls
                               // ! (per each column) x number of columns 
         
         struct CMWWC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cmww;
         };
         
         
          // The values describing the complex permittivity of east walls
                               // ! (per each column) x number of columns 
         
         struct CEEWC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   ceew;
         };
         
         
          // The values describing the complex permeability of east walls
                               // ! (per each column) x number of columns 
         
         struct CMEWC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cmew;
         };
         
         
          // The values describing the complex permittivity of north walls
                               // ! (per each column) x number of columns 
         
         struct CENWC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cenw;
         };
         
         
          // The values describing the complex permeability of north walls
                               // ! (per each column) x number of columns 
         
         struct CMNWC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cmnw;
         };
         
         
          // The values describing the complex permittivity of south angled roof
                               // ! (per each column) x number of columns  
         
         struct CESARC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cesar;
         };
         
         
         // The values describing the complex permeabillity of south angled roof
                               // ! (per each column) x number of columns  
         
         struct CMSARC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cmsar;
         };
         
         
         // The values describing the complex permittivity of east angled roof
                               // ! (per each column) x number of columns  
         
         struct CEEARC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   ceear;
         };
         
         
         // The values describing the complex permeabillity of east angled roof
                               // ! (per each column) x number of columns  
         
         struct CMEARC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cmear;
         };
         
         
         // The values describing the complex permittivity of west angled roof
                               // ! (per each column) x number of columns  
         
         struct CEWARC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cewar;
         };
         
         
         // The values describing the complex permeabillity of west angled roof
                               // ! (per each column) x number of columns  
         
         struct CMWARC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cmwar;
         };
         
         
          // The values describing the complex permittivity of north angled roof
                               // ! (per each column) x number of columns  
         
         struct CENARC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cenar;
         };
         
         
         // The values describing the complex permeabillity of north angled roof
                               // ! (per each column) x number of columns  
         
         struct CMNARC8x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_ymm8c4_t   cmnar;
         };
         
         
         // The components of south walls normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVSWR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m256_t  nvx;
                  DC1D_m256_t  nvy;
                  DC1D_m256_t  nvz;
         };
         
         
        // The components of east walls normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVEWR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m256_t  nvx;
                  DC1D_m256_t  nvy;
                  DC1D_m256_t  nvz;
         };
         
         
        // The components of west walls normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVWWR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m256_t  nvx;
                  DC1D_m256_t  nvy;
                  DC1D_m256_t  nvz;
         };
         
         
        // The components of north walls normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVNWR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m256_t  nvx;
                  DC1D_m256_t  nvy;
                  DC1D_m256_t  nvz;
         };
         
         
         // The components of each building normal vector
                                    // ! (per each column) x number of columns 
         
         struct NVBR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m256_t  nvx;
                  DC1D_m256_t  nvy;
                  DC1D_m256_t  nvz;  
         };
         
         
         // The components of each building flat roof normal vector
                                    // ! (per each column) x number of columns 
         
         struct NVFRR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m256_t  nvx;
                  DC1D_m256_t  nvy;
                  DC1D_m256_t  nvz;   
         };
         
         
          // The components of south angled roof normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVSARR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m256_t  nvx;
                  DC1D_m256_t  nvy;
                  DC1D_m256_t  nvz;
         };
         
         
        // The components of east angled roof normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVEARR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m256_t  nvx;
                  DC1D_m256_t  nvy;
                  DC1D_m256_t  nvz;
         };
         
         
        // The components of west angled roof normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVWARR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m256_t  nvx;
                  DC1D_m256_t  nvy;
                  DC1D_m256_t  nvz;
         };
         
         
        // The components of north angled roof normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVNARR8x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m256_t  nvx;
                  DC1D_m256_t  nvy;
                  DC1D_m256_t  nvz;
         };
         
         
         // The values of each south wall height and width
                        // ! (per each column) x number of columns  
         
         struct HWSWR8x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D_m256_t  hsw;
                DC1D_m256_t  wsw;
         };
         
         
          // The values of each east wall height and width
                        // ! (per each column) x number of columns  
         
         struct HWEWR8x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D_m256_t  hew;
                DC1D_m256_t  wew;
         };
         
         
         // The values of each west wall height and width
                        // ! (per each column) x number of columns  
         
         struct HWWWR8x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D_m256_t  hww;
                DC1D_m256_t  www;
         };
         
         
          // The values of each north wall height and width
                        // ! (per each column) x number of columns  
         
         struct HWNWR8x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D_m256_t  hnw;
                DC1D_m256_t  wnw;
         };
         
         
         // The values of each flat roof height and width
                        // ! (per each column) x number of columns  
         
         struct HWFRR8x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D_m256_t  hfr;
                DC1D_m256_t  wfr;
         };
         
         
        // The values of each south non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         
         struct HWSNFRR8x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D_m256_t  hsnfr;
                DC1D_m256_t  wsnfr;
         };
         
         
        // The values of each east non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         
         struct HWENFRR8x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D_m256_t  henfr;
                DC1D_m256_t  wenfr;
         };
         
         
        // The values of each west non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         
         struct HWWNFRR8x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D_m256_t  hwnfr;
                DC1D_m256_t  wwnfr;
         };
         
         
        // The values of each north non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         
         struct HWNNFRR8x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D_m256_t  hnnfr;
                DC1D_m256_t  wnnfr;
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
         
         struct RCSFRR8x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t  rcsfr;
         };
         
         
          // The values of RCS for the south wall of
         // of every building in the building column.
         
         struct RCSSWR8x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t  rcssw;
         };
         
         
          
          // The values of RCS for the east wall of
         // of every building in the building column.
         
         struct RCSEWR8x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t  rcsew;
         };
         
         
          
          // The values of RCS for the west wall of
         // of every building in the building column.
         
         struct RCSWWR8x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t  rcsww;
         };
         
         
          
          // The values of RCS for the north wall of
         // of every building in the building column.
         
         struct RCSNWR8x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t  rcsnw;
         };
         
         
         // The values of RCS for the south angled roof of
         // of every building in the building column.
         
         struct RCSSARR8x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t  rcssar;
         };
         
         
         // The values of RCS for the south east roof of
         // of every building in the building column.
         
         struct RCSEARR8x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t  rcsear;
         };
         
         
         // The values of RCS for the west angled roof of
         // of every building in the building column.
         
         struct RCSWARR8x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t  rcswar;
         };
         
         
         // The values of RCS for the north angled roof of
         // of every building in the building column.
         
         struct RCSNARR8x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t  rcsnar;
         };
         
         
         // The values of whole building surface area
         // of every building in building column
         
         struct WBSAR8x_t {
                 
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t wbsa;
         };
         
         
         // The values of whole building internal volume
         // of every building in building column
         
         struct WBIVR8x_t {
                 
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m256_t wbiv;
         };
         
         
         
         
           
     }// radiolocation



}






























#endif /*__GMS_URBAN_MODEL_AVX_H__*/
