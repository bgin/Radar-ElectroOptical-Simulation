

#ifndef __GMS_URBAN_MODEL_AVX512_H__
#define __GMS_URBAN_MODEL_AVX512_H__ 311220230825


namespace file_info {

     const unsigned int GMS_URBAN_MODEL_AVX512_MAJOR = 1;
     const unsigned int GMS_URBAN_MODEL_AVX512_MINOR = 0;
     const unsigned int GMS_URBAN_MODEL_AVX512_MICRO = 0;
     const unsigned int GMS_URBAN_MODEL_AVX512_FULLVER =
       1000U*GMS_URBAN_MODEL_AVX512_MAJOR+100U*GMS_URBAN_MODEL_AVX512_MINOR+
       10U*GMS_URBAN_MODEL_AVX512_MICRO;
     const char * const GMS_URBAN_MODEL_AVX512_CREATION_DATE = "31-12-2023 08:25AM +00200 (SUN 31 DEC 2023 GMT+2)";
     const char * const GMS_URBAN_MODEL_AVX512_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_URBAN_MODEL_AVX512_SYNOPSIS      = "Abstract data (AVX512-based) types representing built-up area for diffraction and scattering calculations."

}

#include <cstdint>
#include "GMS_config"
#include "GMS_dyn_containers_avx512.hpp"

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
            
            struct BLatLondR16x_t {
                   
                   // Number of latitude   values (deg), per building
                   std::size_t nblatd;
                   // Number of longtitude values (deg), per building
                   std::size_t nblond;
                   DC1D_m512_t  blatd;
                   DC1D_m512_t  blond;
            };
            
            
            // latitude   values (rad), per building
            // longtitude values (rad), per building
            
            struct BLatLonrR16x_t {
                   
                   // Number of latitude   values (rad), per building
                   std::size_t nblatr;
                   // Number of longtitude values (rad), per building
                   std::size_t nblonr;
                   DC1D_m512_t  blatr;
                   DC1D_m512_t  blonr;
            };
            
            
            // ellipsoidal (radar waveform irradiating field) cells for building column
            
            struct EllpbR16x_t {
                    
                    // Number of ellipsoidal (radar waveform irradiating field) cells for building column
                    std::size_t nellpb;
                    DC1D_m512_t  ellpb;
            };
            
            
            // ! Parametric equation x (acos(t)) values (building)
            // ! 1st dimension building column, 2nd dimension 'x' parameter values
            // ! Parametric equation y (b(sin(t)) values (building)
            // 1st dimension building column, 2nd dimension 'y' parameter values
            
            struct PxybR16x_t {
                   
                   // number of building collumns
                   std::size_t nbpc;
                   std::size_t npxb;
                   std::size_t npyb
                   DC1D_m512_t  pxb;
                   DC1D_m512_t  pyb;
            };
            
            
            // Length, width and an area of every street
            
            struct SLWAR16x_t {
                    
                    // number of streets
                    std::size_t nstr;
                    DC1D_m512_t  lstr;
                    DC1D_m512_t  wstr;
                    DC1D_m512_t  astr;
            };
            
            
            //   ! Moisture of every street (2D array)
            //   ! 2nd dimension humidity values (per street), 1st dimension street numbers
            //    Percent of moist to dry area of evey street at each cell
            
            struct MStrR16x_t {
                    
                    // number of streets
                    std::size_t nstr;
                    std::size_t nmstr;
                    std::size_t npmstr;
                    DC1D_m512_t  mstr;
                    DC1D_m512_t  pmstr;
            };
            
            
            // !Coverage of every street (like: '1' for snow,'2' for mud, '3' for clay, ...etc)
            struct CStrR16x_t {
                    
                   // number of streets
                    std::size_t nstr; 
                    int32_t * __restrict cstr;
            };
            
            
            // Percent of covered to non-covered portion of every street at each irradiated cell
            // Average thickness of each layer (cover) of every street at each irradiated cell
            // Thickness of cover along street (number of values) at each irradiated cell
            
            struct CDStrR16x_t {
                    
                    // number of streets
                    std::size_t nstr;
                    std::size_t npcstr;
                    std::size_t natstr;
                    std::size_t ntcstr; 
                    DC1D_m512_t  pcstr;
                    DC1D_m512_t  atstr;
                    DC1D_m512_t  tcstr;
            };
            
            
            // Mu values for 'clean' street interpolated along the street length at each irradiated cell
            // Eps for 'clean' street street length interpolated at each irradiated cell
            
            struct MEStr1C16x_t   {     
                     // number of streets
                    std::size_t nstr; 
                    std::size_t nmustr;
                    std::size_t nepstr;
                    DC1D_zmm16c4_t   mustr;
                    DC1D_zmm16c4_t   epstr;
            };
            
            
            // Mu for covered (i.e. by mud,snow,clay, ..etc) 
            // street interpolated along the street length at each irradiated cell
            // Eps for covered (i.e. by mud,snow,clay, ..etc) 
            // street  length interpolated at each irradiated cell
            
            struct MEStr2C16x_t {
                     
                      // number of streets
                    std::size_t nstr; 
                    std::size_t nmustr;
                    std::size_t nepstr;
                    DC1D_zmm16c4_t   mustr;
                    DC1D_zmm16c4_t   epstr; 
            };
            
            
            // Street curvature parametric equation u-parameter
            // Street curvature parametric equation v-parameter
            
            struct SCrvR16x_t {
                     
                     // number of streets
                    std::size_t nstr; 
                    std::size_t nustr;
                    std::size_t nvstr;
                    DC1D_m512_t  ustr;
                    DC1D_m512_t  vstr; 
            };
            
            
            //  Street surface normal vectors x-components 
            //  along the street length at each irradiated cell
            //  Street surface normal vectors y-components 
            //  along the street length at each irradiated cell
            //  Street surface normal vectors z-components 
            //  along the street length at each irradiated cell
            
            struct SNrmR16x_t {
                    
                     // number of streets
                    std::size_t nstr; 
                    std::size_t nx;
                    std::size_t ny;
                    std::size_t nz;
                    DC1D_m512_t  nvx;
                    DC1D_m512_t  nvy;
                    DC1D_m512_t  nvz;  
            };
            
            
            // latitude   values (deg), per street length (at irradiance point)
            // longtitude values (deg), per street length (at irradiance point)
            
            struct SIRCDR16x_t {
                     
                    // number of streets
                    std::size_t nstr; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D_m512_t  irlon;
                    DC1D_m512_t  irlat;  
            } 
            
            
           // latitude   values (rad), per street length (at irradiance point)
           // longtitude values (rad), per street length (at irradiance point)
           
           struct SIRCRR16x_t {
                    
                   // number of streets
                    std::size_t nstr; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D_m512_t  irlon;
                    DC1D_m512_t  irlat;  
           };
           
           
           // latitude   values (deg), of the building area (at irradiance point)
           // longtitude values (deg), of the building area (at irradiance point)
            
            struct BIRCDR16x_t {
                     
                    // number of buildings
                    std::size_t nbld; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D_m512_t  irlon;
                    DC1D_m512_t  irlat;  
            };
            
            
           // latitude   values (rad), of the building area (at irradiance point)
           // longtitude values (rad), of the building area (at irradiance point)
           
           struct BIRCRR16x_t {
                    
                   // number of buildings
                    std::size_t nbld; 
                    std::size_t nlon;
                    std::size_t nlat;
                    DC1D_m512_t  irlon;
                    DC1D_m512_t  irlat;  
           }; 
           
           
           // Urban area height map (at single building resolution)
           
           struct UHMapR16x_t {
                  
                  std::size_t nx;
                  std::size_t ny;
                  DC1D_m512_t  hmap; 
           };
           
           
           // Urban area height map (at single building resolution) -- 1st derivative
           
           struct UHDxDyR16x_t {
                  
                  std::size_t nx;
                  std::size_t ny;
                  DC1D_m512_t  hdxdy;
           };
           
           
           // Urban area height map (at single building resolution) -- gradient x-component
           // Urban area height map (at single building resolution) -- gradient y-component
           
           struct UHGradR16x_t {
                  
                  std::size_t nx;
                  std::size_t ny;
                  DC1D_m512_t  uhgx;
                  DC1D_m512_t  uhgy;
           };
           
           
           // Smoothing and approximating curve for linearly-piecewise height function (x-coordinate)
           // Smoothing and approximating curve for linearly-piecewise height function (y-coordinate)
           
           struct XYSMBHR16x_t {
                   
                  std::size_t nx;
                  std::size_t ny;
                  DC1D_m512_t  xsmbh;
                  DC1D_m512_t  ysmbh
           };
           
           
           // Empty space in-between of buildings (per single column) x number columns
           struct ESBBI1x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  int32_t * __restrict esbb;
           };
           
           
           // An area values of in-between buildings empty spaces (per single column) x number columns
            
           struct AESBBR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval;
                  DC1D_m512_t  aesbb;
           };
           
           
           // An area values of each building (per single building column) x number columns
           
           struct ABCR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  abc;
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
           
           struct BRAPCR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  brapc;
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
           
           struct IDSRWR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  idsrw;
           };
           
           
           // An angled roof inclination (deg) -- east facing roof wall 
           // (per each column) x number of columns
           
           struct IDERWR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  iderw;
           };
           
           
           // An angled roof inclination (deg) -- west facing roof wall 
           // (per each column) x number of columns
           
           struct IDWRWR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  idwrw;
           };
           
           
           // An angled roof inclination (rad) -- south facing roof wall 
           // (per each column) x number of columns
           
           struct IDNRWR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  idnrw; 
           };
           
           
           //  An angled roof inclination (rad) -- south facing roof wall 
           //  (per each column) x number of columns
           
           struct IRSRWR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  irsrw; 
           };
           
           
           //  An angled roof inclination (rad) -- east facing roof wall 
           //  (per each column) x number of columns
           
           struct IRERWR16x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  irerw;
           };
           
           
           //  An angled roof inclination (rad) -- north facing roof wall 
           //  (per each column) x number of columns
           
           struct IRNRWR16x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  irnrw; 
           };
           
           
           // An angled roof inclination surface area -- south facing roof wall 
           // (per each column) x number of columns
           
           struct ISRAR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  isra;  
           };
           
           
           // An angled roof inclination surface area -- west facing roof wall 
           // (per each column) x number of columns
           
           struct IWRAR16x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  iwra;  
           };
           
           
           // An angled roof inclination surface area -- east facing roof wall 
           // (per each column) x number of columns
           
           struct IERAR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  iera;   
           };
           
           
           // An angled roof inclination surface area -- north facing roof wall 
           // (per each column) x number of columns
           
           struct INRAR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  inra;  
           };
           
           
           // South wall upper-facing edge inclination (rad) -- 
           // (per each column)  x number of columns
          
          struct SWUER16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  swue;   
          };
          
          
          // East wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          
          struct EWUER16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  ewue;   
          };
          
          
          // West wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          
          struct WWUER16x_t {
                  
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  wwue;   
          };
          
          
          // North wall upper-facing edge inclination (rad) -- 
          // (per each column)  x number of columns
          
          struct NWUER16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  nwue;    
          };
          
          
          // Shared right edge between the south wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          
          struct SEWER16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  sewe; 
          };
          
          // Shared left edge between the south wall and west wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
          
          struct SWWER16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  swwe;
          };
          
          
           // Shared right edge between the north wall and east wall inclination (rad) 
                                    // ! -- (per each column) x number of columns
           
           struct NWEER16x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  nwee;
           };
           
           
           // Shared right edge between the north wall and west wall inclination (rad) 
           // ! -- (per each column) x number of columns
           
           struct NWWER16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  nwwe; 
           };
           
           
           // Simple cell-based mesh
           
           struct CellMeshR16x_t {
                  
                  // Number of divisions along the x,y,z
                  int32_t ndiv[3];
                  // Number of cells
                  std::size_t nL;
                  // Compute numerical integration.
                  bool nint;
                  // Coordinates (x,y,z) of the center
                  // of Lth cell.
                  DC1D_m512_t cx;
                  DC1D_m512_t cy;
                  DC1D_m512_t cz;
                  DC1D_m512_t dv;
                  // (X,Y,Z) dimensions of the Ith
                  // rectangular volume cell (this is needed for
                  // the numerical integration)
                  DC1D_m512_t dx;
                  DC1D_m512_t dy;
                  DC1D_m512_t dz;
           };
           
           
           // South walls surface area (for every building, per column) x number of columns
           
           struct SWSAR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  swsa;
           };
           
           
           // East walls surface area (for every building, per column) x number of columns
           
           struct EWSAR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  ewsa; 
           };
           
           
           // West walls surface area (for every building, per column) x number of columns
           
           struct WWSAR16x_t {
                   
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  wwsa;  
           };
           
           
           // North walls surface area (for every building, per column) x number of columns
           
           struct NWSAR16x_t {
                    
                   std::size_t ncols;
                   std::size_t nval;
                   DC1D_m512_t  nwsa;   
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
           
           struct MDSWRR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mdswr;
           };
           
           
           // ! The values describing the ratio (percentage) of east wall 
                            // ! moisture to dryness (per each column) x number of columns
           
           struct MDEWRR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mdewr;
           };
           
           
           //  The values describing the ratio (percentage) of west wall 
                             //! moisture to dryness (per each column) x number of columns
           
           struct MDWWRR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mdwwr;
           };  
           
           
           //  The values describing the ratio (percentage) of north wall 
                             //! moisture to dryness (per each column) x number of columns              
           
           struct MDNWRR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mdnwr;
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
          
          struct MDRRR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mdrr; 
          };
          
          
          // The values describing the surface of moist part of the flat roof 
                               // ! (per each column) x number of columns
          
          struct MPFRR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mpfr; 
          };     
          
          
          // The values describing the surface of dry part of the flat roof 
                               // ! (per each column) x number of columns  
          
          struct DPFRR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  dpfr; 
          };   
          
          
          //  The values describing the surface of moist part of the south wall 
                                // ! (per each column) x number of columns   
          
          struct MPSWR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mpsw; 
          };  
          
          
          //  The values describing the surface of dry part of the south wall 
                               //  ! (per each column) x number of columns  
          
          struct DPSWR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  dpsw; 
          };   
          
          
          //  The values describing the surface of moist part of the east wall 
                                // ! (per each column) x number of columns
          
          struct MPEWR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mpew; 
          };   
          
         // The values describing the surface of dry part of the east wall 
                                // ! (per each column) x number of columns 
          
          struct DPEWR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  dpew; 
          };  
          
          
         // The values describing the surface of moist part of the west wall 
                                 //! (per each column) x number of columns
          
          struct MPWWR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mpww; 
          }; 
          
          
        //  The values describing the surface of dry part of the west wall 
                                 //! (per each column) x number of columns 
          
          struct DPWWR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  dpww; 
          };  
          
          
        // The values describing the surface of moist part of the north wall 
                                 //! (per each column) x number of columns
         
         struct MPNWR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mpnw; 
          }; 
          
          
         //  The values describing the surface of dry part of the north wall 
                                // ! (per each column) x number of columns
         
         struct DPNWR16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  dpnw; 
          }; 
          
          
          // The values describing the surface of moist part of the angled south roof wall
                               // ! (per each column) x number of columns
          
          struct MPSARR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mpsar;
          };
          
          
          // The values describing the surface of dry part of the angled south roof wall
                               // ! (per each column) x number of columns
          
          struct DPSARR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  dpsar;
          };
          
          
          // The values describing the surface of moist part of the angled east roof wall
                               // ! (per each column) x number of columns 
          
          struct MPEARR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mpear;
          };
          
          
         // The values describing the surface of dry part of the angled east roof wall
                               // ! (per each column) x number of columns  
         
         struct DPEARR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  dpear;
          }; 
          
          
          // The values describing the surface of moist part of the angled west roof wall
                               // ! (per each column) x number of columns  
          
          struct MPWARR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mpwar;
          }; 
          
          
           // The values describing the surface of dry part of the angled west roof wall
                               // ! (per each column) x number of columns  
          
          struct DPWARR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  dpwar;
          }; 
          
          
           // The values describing the surface of moist part of the angled north roof wall
                               // ! (per each column) x number of columns  
          
          struct MPNARR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  mpnar;
          }; 
          
          
          // The values describing the surface of dry part of the angled north roof wall
                               // ! (per each column) x number of columns  
          
          struct DPNARR16x_t {
                   
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_m512_t  dpnar;
          }; 
          
          
         // The values describing the complex permittivity of south walls
                               // ! (per each column) x number of columns  
         
         struct CESWC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cesw;
         };
         
         
         // The values describing the complex permeabillity of south walls
                               // ! (per each column) x number of columns  
         
         struct CMSWC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cmsw;
         };
         
         
         // The values describing the complex permittivity of west walls
                               // ! (per each column) x number of columns 
         
         struct CEWWC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   ceww;
         };
         
         
          // The values describing the complex permeability of west walls
                               // ! (per each column) x number of columns 
         
         struct CMWWC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cmww;
         };
         
         
          // The values describing the complex permittivity of east walls
                               // ! (per each column) x number of columns 
         
         struct CEEWC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   ceew;
         };
         
         
          // The values describing the complex permeability of east walls
                               // ! (per each column) x number of columns 
         
         struct CMEWC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cmew;
         };
         
         
          // The values describing the complex permittivity of north walls
                               // ! (per each column) x number of columns 
         
         struct CENWC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cenw;
         };
         
         
          // The values describing the complex permeability of north walls
                               // ! (per each column) x number of columns 
         
         struct CMNWC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cmnw;
         };
         
         
          // The values describing the complex permittivity of south angled roof
                               // ! (per each column) x number of columns  
         
         struct CESARC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cesar;
         };
         
         
         // The values describing the complex permeabillity of south angled roof
                               // ! (per each column) x number of columns  
         
         struct CMSARC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cmsar;
         };
         
         
         // The values describing the complex permittivity of east angled roof
                               // ! (per each column) x number of columns  
         
         struct CEEARC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   ceear;
         };
         
         
         // The values describing the complex permeabillity of east angled roof
                               // ! (per each column) x number of columns  
         
         struct CMEARC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cmear;
         };
         
         
         // The values describing the complex permittivity of west angled roof
                               // ! (per each column) x number of columns  
         
         struct CEWARC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cewar;
         };
         
         
         // The values describing the complex permeabillity of west angled roof
                               // ! (per each column) x number of columns  
         
         struct CMWARC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cmwar;
         };
         
         
          // The values describing the complex permittivity of north angled roof
                               // ! (per each column) x number of columns  
         
         struct CENARC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cenar;
         };
         
         
         // The values describing the complex permeabillity of north angled roof
                               // ! (per each column) x number of columns  
         
         struct CMNARC16x_t {
                  
                  std::size_t ncols;
                  std::size_t nval;
                  DC1D_zmm16c4_t   cmnar;
         };
         
         
         // The components of south walls normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVSWR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m512_t  nvx;
                  DC1D_m512_t  nvy;
                  DC1D_m512_t  nvz;
         };
         
         
        // The components of east walls normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVEWR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m512_t  nvx;
                  DC1D_m512_t  nvy;
                  DC1D_m512_t  nvz;
         };
         
         
        // The components of west walls normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVWWR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m512_t  nvx;
                  DC1D_m512_t  nvy;
                  DC1D_m512_t  nvz;
         };
         
         
        // The components of north walls normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVNWR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m512_t  nvx;
                  DC1D_m512_t  nvy;
                  DC1D_m512_t  nvz;
         };
         
         
         // The components of each building normal vector
                                    // ! (per each column) x number of columns 
         
         struct NVBR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m512_t  nvx;
                  DC1D_m512_t  nvy;
                  DC1D_m512_t  nvz;  
         };
         
         
         // The components of each building flat roof normal vector
                                    // ! (per each column) x number of columns 
         
         struct NVFRR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m512_t  nvx;
                  DC1D_m512_t  nvy;
                  DC1D_m512_t  nvz;   
         };
         
         
          // The components of south angled roof normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVSARR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m512_t  nvx;
                  DC1D_m512_t  nvy;
                  DC1D_m512_t  nvz;
         };
         
         
        // The components of east angled roof normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVEARR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m512_t  nvx;
                  DC1D_m512_t  nvy;
                  DC1D_m512_t  nvz;
         };
         
         
        // The components of west angled roof normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVWARR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m512_t  nvx;
                  DC1D_m512_t  nvy;
                  DC1D_m512_t  nvz;
         };
         
         
        // The components of north angled roof normal vector
                                    // ! (per each column) x number of columns  
         
         struct NVNARR16x_t {
                  
                  std::size_t ncols
                  std::size_t nval
                  DC1D_m512_t  nvx;
                  DC1D_m512_t  nvy;
                  DC1D_m512_t  nvz;
         };
         
         
         // The values of each south wall height and width
                        // ! (per each column) x number of columns  
         
         struct HWSWR16x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D_m512_t  hsw;
                DC1D_m512_t  wsw;
         };
         
         
          // The values of each east wall height and width
                        // ! (per each column) x number of columns  
         
         struct HWEWR16x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D_m512_t  hew;
                DC1D_m512_t  wew;
         };
         
         
         // The values of each west wall height and width
                        // ! (per each column) x number of columns  
         
         struct HWWWR16x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D_m512_t  hww;
                DC1D_m512_t  www;
         };
         
         
          // The values of each north wall height and width
                        // ! (per each column) x number of columns  
         
         struct HWNWR16x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D_m512_t  hnw;
                DC1D_m512_t  wnw;
         };
         
         
         // The values of each flat roof height and width
                        // ! (per each column) x number of columns  
         
         struct HWFRR16x_t {
                
                std::size_t ncols;
                std::size_t nval;
                DC1D_m512_t  hfr;
                DC1D_m512_t  wfr;
         };
         
         
        // The values of each south non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         
         struct HWSNFRR16x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D_m512_t  hsnfr;
                DC1D_m512_t  wsnfr;
         };
         
         
        // The values of each east non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         
         struct HWENFRR16x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D_m512_t  henfr;
                DC1D_m512_t  wenfr;
         };
         
         
        // The values of each west non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         
         struct HWWNFRR16x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D_m512_t  hwnfr;
                DC1D_m512_t  wwnfr;
         };
         
         
        // The values of each north non-flat roof (either triangular, squared or rectangulart) 
        // height and width (base) (per each column) x number of columns  
         
         struct HWNNFRR16x_t {
                
                std::size_t ncols;
                std::size_t nval;
                int32_t     type; // 0 triangular roof , 1 squared roof , 2 reactangular roof
                DC1D_m512_t  hnnfr;
                DC1D_m512_t  wnnfr;
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
         
         struct RCSFRR16x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t  rcsfr;
         };
         
         
          // The values of RCS for the south wall of
         // of every building in the building column.
         
         struct RCSSWR16x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t  rcssw;
         };
         
         
          
          // The values of RCS for the east wall of
         // of every building in the building column.
         
         struct RCSEWR16x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t  rcsew;
         };
         
         
          
          // The values of RCS for the west wall of
         // of every building in the building column.
         
         struct RCSWWR16x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t  rcsww;
         };
         
         
          
          // The values of RCS for the north wall of
         // of every building in the building column.
         
         struct RCSNWR16x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t  rcsnw;
         };
         
         
         // The values of RCS for the south angled roof of
         // of every building in the building column.
         
         struct RCSSARR16x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t  rcssar;
         };
         
         
         // The values of RCS for the south east roof of
         // of every building in the building column.
         
         struct RCSEARR16x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t  rcsear;
         };
         
         
         // The values of RCS for the west angled roof of
         // of every building in the building column.
         
         struct RCSWARR16x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t  rcswar;
         };
         
         
         // The values of RCS for the north angled roof of
         // of every building in the building column.
         
         struct RCSNARR16x_t {
                
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t  rcsnar;
         };
         
         
         // The values of whole building surface area
         // of every building in building column
         
         struct WBSAR16x_t {
                 
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t wbsa;
         };
         
         
         // The values of whole building internal volume
         // of every building in building column
         
         struct WBIVR16x_t {
                 
                 std::size_t ncols;
                 std::size_t nval;
                 DC1D_m512_t wbiv;
         };
         
         
         
         
           
     }// radiolocation



}






























#endif /*__GMS_URBAN_MODEL_AVX512_H__*/
