
#include "GMS_config.fpp"


module antenna_model_adt



!===================================================================================85
 !---------------------------- DESCRIPTION ------------------------------------------85
 !
 !
 !
 !          Module  name:
 !                         antenna_model_adt
 !          
 !          Purpose:
 !                        Module types for 'antenna_sensor'  implementation.
 !                        Various characteristics of different antenna types  
 !                        Based mainly on book titled (rus):          
 !                        Проектирование антенно фидерных устройств. Жук М.С. Молочков Ю.Б
 !          History:
 !                        Date: 18-11-2023
 !                        Time: 09:53 GMT+2
 !                        
 !          Version:
 !
 !                      Major: 1
 !                      Minor: 0
 !                      Micro: 0
 !
 !          Author:  
 !                      Bernard Gingold
 !          
 !                 
 !          References:
 !         
 !                      Проектирование антенно фидерных устройств. Жук М.С. Молочков Ю.Б      
 !         
 !          E-mail:
 !                  
 !                      beniekg@gmail.com
!==================================================================================85
    ! Tab:5 col - Type and etc.. definitions
    ! Tab:10,11 col - Type , function and subroutine code blocks.
   
    use mod_kinds,    only : i4,sp,dp

    public
    implicit none

    ! Major version
    integer(kind=i4),  parameter :: ANTENNA_MODEL_ADT_MAJOR = 1
    ! Minor version
    integer(kind=i4),  parameter :: ANTENNA_MODEL_ADT_MINOR = 0
    ! Micro version
    integer(kind=i4),  parameter :: ANTENNA_MODEL_ADT_MICRO = 0
    ! Full version
    integer(kind=i4),  parameter :: ANTENNA_MODEL_ADT_FULLVER =   &
            1000*ANTENNA_MODEL_ADT_MAJOR+100*ANTENNA_MODEL_ADT_MINOR+10*ANTENNA_MODEL_ADT_MICRO
    ! Module creation date
    character(*),        parameter :: ANTENNA_MODEL_ADT_CREATE_DATE = "18-11-2023 09:53 +00200 (SAT 18 NOV 20223 GMT+2)"
    ! Module build date
    character(*),        parameter :: ANTENNA_MODEL_ADT_BUILD_DATE  = __DATE__ " " __TIME__
    ! Module author info
    character(*),        parameter :: ANTENNA_MODEL_ADT_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com"
    ! Short description
    character(*),        parameter :: ANTENNA_MODEL_ADT_SYNOPSIS    = "ADT describing various antenna types characteristics -- module based."

    
    ! Array holding a number of points for various antenna characteristics computation.
    integer(kind=i4) :: ne      ! number of  Electric fields (for phased array radiating element)
    integer(kind=i4) :: nepts1  !x-dim
    integer(kind=i4) :: nepts2  !y-dim
    integer(kind=i4) :: nepts3  !z-dim
    
    integer(kind=i4) :: nm      ! number of  Magnetic fields (for phased array radiating element)
    integer(kind=i4) :: nmpts1  !x-dim
    integer(kind=i4) :: nmpts2  !y-dim
    integer(kind=i4) :: nmpts3  !z-dim
    
    integer(kind=i4) :: nesp    ! number of Electric fields spherical-coordinates (for phased array radiating element)
    integer(kind=i4) :: neptsph ! phi-coordinate
    integer(kind=i4) :: neptsth ! theta-coordinate
    
    integer(kind=i4) :: nmsp    ! number of Magnetic fields spherical-coordinates (for phased array radiating element)
    integer(kind=i4) :: nmptsph ! phi-coordinate
    integer(kind=i4) :: nmptsth ! theta-coordinate
    
    integer(kind=i4) :: nje     ! number of surface currents (electric) (for phased array radiating element)
    integer(kind=i4) :: njepts1 !x-dim
    integer(kind=i4) :: njepts2 !y-dim
    integer(kind=i4) :: njepts3 !z-dim 
    
    integer(kind=i4) :: njm     ! number of surface currents (magnetic) (for phased array radiating element)
    integer(kind=i4) :: njmpts1 !x-dim
    integer(kind=i4) :: njmpts2 !y-dim
    integer(kind=i4) :: njmpts3 !z-dim
    
    integer(kind=i4) :: nhe     ! number of Hertz (electric) vectors
    integer(kind=i4) :: nhepts1 !x-dim
    integer(kind=i4) :: nhepts2 !y-dim
    integer(kind=i4) :: nhepts3 !z-dim
    
    integer(kind=i4) :: nhm     ! number of Hertz (magnetic) vectors
    integer(kind=i4) :: nhmpts1 !x-dim
    integer(kind=i4) :: nhmpts2 !y-dim
    integer(kind=i4) :: nhmpts3 !z-dim
    
    integer(kind=i4) :: nith    ! number of theta integration points
    integer(kind=i4) :: niph    ! number of phi integration points
    integer(kind=i4) :: noth    ! number of theta observation points
    integer(kind=i4) :: noph    ! number of phi obsertvation points
    
    integer(kind=i4) :: nsur    ! number of R components of spherical coordinate unit vector
    integer(kind=i4) :: nsup    ! number of Phi components of spherical coordinate unit vector
    integer(kind=i4) :: nsut    ! number of Theta components of spherical coordinate unit vector
    
    integer(kind=i4) :: nrx    ! number of x components of normal vectors
    integer(kind=i4) :: nry    ! number of y components of normal vectors
    integer(kind=i4) :: nrz    ! number of z components of normal vectors
    
    integer(kind=i4) :: nirx    ! number of Rx elements of integration points.
    integer(kind=i4) :: niry    ! number of Ry elements of integration points.
    integer(kind=i4) :: nirz    ! number of Rz elements of integration points.
    
    integer(kind=i4) :: norx    ! number of Rx elements of observation points.
    integer(kind=i4) :: nory    ! number of Ry elements of observation points.
    integer(kind=i4) :: norz    ! number of Rz elements of observation points.
    
    integer(kind=i4) :: npsip   ! number of Psi function's phi values
    integer(kind=i4) :: npsit   ! number of Psi function's theta values
    integer(kind=i4) :: npsi    ! number of Psi functions for phased arrays.
    
    integer(kind=i4) :: nswapp  ! numbers of terms for spherical wave approximation
    integer(kind=i4) :: ncost   ! number of values of cos(theta) argument (formulae: 2.52,2.53,...etc)
    integer(kind=i4) :: nrho    ! distances from the antenna apperture to the coordinate system origin (2.53)
    integer(kind=i4) :: nsint   ! number of sin(theta) values (2.62)
    integer(kind=i4) :: nsinp   ! number of sin(phi) values (2.62)
    integer(kind=i4) :: cosp    ! number of cos(phi) values (2.62)
    integer(kind=i4) :: nx      ! number of 'x' values antenna apperture (2.62)
    integer(kind=i4) :: ny      ! number of 'y' values antenna apperture (2.62)
    integer(kind=i4) :: nL      ! number of 'L' sizes of apperture (2.65)
    integer(kind=i4) :: ngam    ! number of wavelengths (2.65)
    integer(kind=i4) :: nftf265 ! number of values i.e. F(theta) for radiation pattern (2.65)
    integer(kind=i4) :: nftf268 ! number of values i.e. F(theta) for radiation pattern (2.68) 
    integer(kind=i4) :: nu      ! number of 'u' values (2.70)
    integer(kind=i4) :: mf271   ! number of 'm' terms (2.71)
    integer(kind=i4) :: nf271   ! value of 'n' (2.71)
    integer(kind=i4) :: nftf271 ! number of values i.e. F(theta) for radiation pattern (2.71)
    integer(kind=i4) :: nftf269 ! number of values i.e. f(2x/L) for radiation pattern (2.69)
    integer(kind=i4) :: nftf274 ! number of values i.e. f(2x/L) for radiation pattern (2.74)
    integer(kind=i4) :: nftf275 ! number of values i.e. F(theta) for radiation pattern (2.75)
    integer(kind=i4) :: nftf278 ! number of values i.e. f(x) for radiation pattern (2.78)
    integer(kind=i4) :: nNxf283 ! number of Nx function component (2.83)
    integer(kind=i4) :: nrf283  ! number of 'r' coordinate of apperture (2.83)
    integer(kind=i4) :: npf283  ! number of 'phi' coordinate of apperture (2.83)
    integer(kind=i4) :: nuf286  ! number of 'u' values (2.86)
    integer(kind=i4) :: nfuf285 ! number of values i.e. F(u)vfor radiation pattern (2.85)
    integer(kind=i4) :: nauf286 ! number of values of function (2.86)
    integer(kind=i4) :: nrhf286 ! number of values 'rho' argument (2.86)
    integer(kind=i4) :: nrf287  ! number of values rho coordinate f(rho,phi) (2.87)
    integer(kind=i4) :: npf287  ! number of values phi coordinate f(rho,phi) (2.87)
    integer(kind=i4) :: nftf291 ! number of values i.e. F(theta) for radiation pattern (2.91)
    integer(kind=i4) :: nmf291  ! number of m-1 values (2.91)
    integer(kind=i4) :: nfrf293 ! number of values i.e. F(rho) for radiation pattern (2.93)
    integer(kind=i4) :: nmf293  ! number of terms (2.93)
    integer(kind=i4) :: nrf293  ! number of 'r' coordinate values (2.93)
    integer(kind=i4) :: nfxf294 ! number of values i.e. f(x,y) for amplitude field of aperture (2.94) 
    
    
    ! Complex Electric Field (complex-single)
    ! First dimension nth field, second dimension number of sample points
    complex(kind=sp), dimension(:,:), allocatable :: exc4
    complex(kind=sp), dimension(:,:), allocatable :: eyc4
    complex(kind=sp), dimension(:,:), allocatable :: ezc4
    !dir$ attributes align : 64 :: exc4
    !dir$ attributes align : 64 :: eyc4
    !dir$ attributes align : 64 :: ezc4
    
    ! Complex Electric Field (complex-double)
     ! First dimension nth field, second dimension number of sample points
    complex(kind=dp), dimension(:,:), allocatable :: exc8
    complex(kind=dp), dimension(:,:), allocatable :: eyc8
    complex(kind=dp), dimension(:,:), allocatable :: ezc8
    !dir$ attributes align : 64 :: exc8
    !dir$ attributes align : 64 :: eyc8
    !dir$ attributes align : 64 :: ezc8
    
    ! Complex Electric Field (real-single) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=sp), dimension(:,:), allocatable :: ex4r
    real(kind=sp), dimension(:,:), allocatable :: ex4i
    real(kind=sp), dimension(:,:), allocatable :: ey4r
    real(kind=sp), dimension(:,:), allocatable :: ey4i
    real(kind=sp), dimension(:,:), allocatable :: ez4r
    real(kind=sp), dimension(:,:), allocatable :: ez4i
    !dir$ attributes align : 64 :: ex4r
    !dir$ attributes align : 64 :: ex4i
    !dir$ attributes align : 64 :: ey4r
    !dir$ attributes align : 64 :: ey4i
    !dir$ attributes align : 64 :: ez4r
    !dir$ attributes align : 64 :: ez4i

    ! Complex Electric Field (real-double) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=dp), dimension(:,:), allocatable :: ex8r
    real(kind=dp), dimension(:,:), allocatable :: ex8i
    real(kind=dp), dimension(:,:), allocatable :: ey8r
    real(kind=dp), dimension(:,:), allocatable :: ey8i
    real(kind=dp), dimension(:,:), allocatable :: ez8r
    real(kind=dp), dimension(:,:), allocatable :: ez8i
    !dir$ attributes align : 64 :: ex8r
    !dir$ attributes align : 64 :: ex8i
    !dir$ attributes align : 64 :: ey8r
    !dir$ attributes align : 64 :: ey8i
    !dir$ attributes align : 64 :: ez8r
    !dir$ attributes align : 64 :: ez8i
    
    ! Complex Electric Field (complex-single) (theta,phi coordinates)
    ! First dimension nth field, second dimension number of sample points 
    complex(kind=sp), dimension(:,:), allocatable :: etc4
    complex(kind=sp), dimension(:,:), allocatable :: epc4
    !dir$ attributes align : 64 :: etc4
    !dir$ attributes align : 64 :: epc4
    
    ! Complex Electric Field (complex-double) (theta,phi coordinates)
    ! First dimension nth field, second dimension number of sample points 
    complex(kind=dp), dimension(:,:), allocatable :: etc8
    complex(kind=dp), dimension(:,:), allocatable :: epc8
    !dir$ attributes align : 64 :: etc8
    !dir$ attributes align : 64 :: epc8
    
    ! Complex Electric Field (real-single) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=sp), dimension(:,:), allocatable :: etr4
    real(kind=sp), dimension(:,:), allocatable :: eti4
    real(kind=sp), dimension(:,:), allocatable :: epr4
    real(kind=sp), dimension(:,:), allocatable :: epi4
    !dir$ attributes align : 64 :: etr4
    !dir$ attributes align : 64 :: eti4
    !dir$ attributes align : 64 :: epr4
    !dir$ attributes align : 64 :: epi4
    
    
    ! Complex Electric Field (real-double) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=dp), dimension(:,:), allocatable :: etr8
    real(kind=dp), dimension(:,:), allocatable :: eti8
    real(kind=dp), dimension(:,:), allocatable :: epr8
    real(kind=dp), dimension(:,:), allocatable :: epi8
    !dir$ attributes align : 64 :: etr8
    !dir$ attributes align : 64 :: eti8
    !dir$ attributes align : 64 :: epr8
    !dir$ attributes align : 64 :: epi8
    
       
    ! Complex Magnetic Field (complex-single)
    ! First dimension nth field, second dimension number of sample points
    complex(kind=sp), dimension(:,:), allocatable :: hxc4
    complex(kind=sp), dimension(:,:), allocatable :: hyc4
    complex(kind=sp), dimension(:,:), allocatable :: hzc4
    !dir$ attributes align : 64 :: hxc4
    !dir$ attributes align : 64 :: hyc4
    !dir$ attributes align : 64 :: hzc4
    
    ! Complex Magnetic Field (complex-double)
     ! First dimension nth field, second dimension number of sample points
    complex(kind=dp), dimension(:,:), allocatable :: hxc8
    complex(kind=dp), dimension(:,:), allocatable :: hyc8
    complex(kind=dp), dimension(:,:), allocatable :: hzc8
    !dir$ attributes align : 64 :: hxc8
    !dir$ attributes align : 64 :: hyc8
    !dir$ attributes align : 64 :: hzc8
    
    ! Complex Magnetic Field (real-single) (decomposed)
     ! First dimension nth field, second dimension number of sample points
    real(kind=sp), dimension(:,:), allocatable :: hx4r
    real(kind=sp), dimension(:,:), allocatable :: hx4i
    real(kind=sp), dimension(:,:), allocatable :: hy4r
    real(kind=sp), dimension(:,:), allocatable :: hy4i
    real(kind=sp), dimension(:,:), allocatable :: hz4r
    real(kind=sp), dimension(:,:), allocatable :: hz4i
    !dir$ attributes align : 64 :: hx4r
    !dir$ attributes align : 64 :: hx4i
    !dir$ attributes align : 64 :: hy4r
    !dir$ attributes align : 64 :: hy4i
    !dir$ attributes align : 64 :: hz4r
    !dir$ attributes align : 64 :: hz4i

    ! Complex Magnetic Field (real-double) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=dp), dimension(:,:), allocatable :: hx8r
    real(kind=dp), dimension(:,:), allocatable :: hx8i
    real(kind=dp), dimension(:,:), allocatable :: hy8r
    real(kind=dp), dimension(:,:), allocatable :: hy8i
    real(kind=dp), dimension(:,:), allocatable :: hz8r
    real(kind=dp), dimension(:,:), allocatable :: hz8i
    !dir$ attributes align : 64 :: hx8r
    !dir$ attributes align : 64 :: hx8i
    !dir$ attributes align : 64 :: hy8r
    !dir$ attributes align : 64 :: hy8i
    !dir$ attributes align : 64 :: hz8r
    !dir$ attributes align : 64 :: hz8i
    
    ! Complex Magnetic Field (complex-single) (theta,phi coordinates)
    ! First dimension nth field, second dimension number of sample points 
    complex(kind=sp), dimension(:,:), allocatable :: htc4
    complex(kind=sp), dimension(:,:), allocatable :: hpc4
    !dir$ attributes align : 64 :: htc4
    !dir$ attributes align : 64 :: hpc4
    
    ! Complex Magnetic Field (complex-double) (theta,phi coordinates)
    ! First dimension nth field, second dimension number of sample points 
    complex(kind=dp), dimension(:,:), allocatable :: htc8
    complex(kind=dp), dimension(:,:), allocatable :: hpc8
    !dir$ attributes align : 64 :: htc8
    !dir$ attributes align : 64 :: hpc8
    
    ! Complex Magnetic Field (real-single) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=sp), dimension(:,:), allocatable :: htr4
    real(kind=sp), dimension(:,:), allocatable :: hti4
    real(kind=sp), dimension(:,:), allocatable :: hpr4
    real(kind=sp), dimension(:,:), allocatable :: hpi4
    !dir$ attributes align : 64 :: htr4
    !dir$ attributes align : 64 :: hti4
    !dir$ attributes align : 64 :: hpr4
    !dir$ attributes align : 64 :: hpi4
    
    
    ! Complex Magnetic Field (real-double) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=dp), dimension(:,:), allocatable :: htr8
    real(kind=dp), dimension(:,:), allocatable :: hti8
    real(kind=dp), dimension(:,:), allocatable :: hpr8
    real(kind=dp), dimension(:,:), allocatable :: hpi8
    !dir$ attributes align : 64 :: htr8
    !dir$ attributes align : 64 :: hti8
    !dir$ attributes align : 64 :: hpr8
    !dir$ attributes align : 64 :: hpi8
    
    ! Antenna apperture surface vector function 'N'
    ! complex (single)
    complex(kind=sp), dimension(:,:), allocatable :: nxc4
    complex(kind=sp), dimension(:,:), allocatable :: nyc4
    complex(kind=sp), dimension(:,:), allocatable :: nzc4
    !dir$ attributes align : 64 :: nxc4
    !dir$ attributes align : 64 :: nyc4
    !dir$ attributes align : 64 :: nzc4
    
    ! Antenna apperture surface vector function 'N'
    ! complex (double)
    complex(kind=dp), dimension(:,:), allocatable :: nxc8
    complex(kind=dp), dimension(:,:), allocatable :: nyc8
    complex(kind=dp), dimension(:,:), allocatable :: nzc8
    !dir$ attributes align : 64 :: nxc8
    !dir$ attributes align : 64 :: nyc8
    !dir$ attributes align : 64 :: nzc8
    
    ! Antenna apperture surface vector function 'N'
    ! complex (real-single) 
    real(kind=sp), dimension(:,:), allocatable :: nx4r
    real(kind=sp), dimension(:,:), allocatable :: nx4i
    real(kind=sp), dimension(:,:), allocatable :: ny4r
    real(kind=sp), dimension(:,:), allocatable :: ny4i
    real(kind=sp), dimension(:,:), allocatable :: nz4r
    real(kind=sp), dimension(:,:), allocatable :: nz4i
    !dir$ attributes align : 64 :: nx4r
    !dir$ attributes align : 64 :: nx4i
    !dir$ attributes align : 64 :: ny4r
    !dir$ attributes align : 64 :: ny4i
    !dir$ attributes align : 64 :: nz4r
    !dir$ attributes align : 64 :: nz4i

    ! Antenna apperture surface vector function 'N'
    ! complex (real-double) 
    real(kind=dp), dimension(:,:), allocatable :: nx8r
    real(kind=dp), dimension(:,:), allocatable :: nx8i
    real(kind=dp), dimension(:,:), allocatable :: ny8r
    real(kind=dp), dimension(:,:), allocatable :: ny8i
    real(kind=dp), dimension(:,:), allocatable :: nz8r
    real(kind=dp), dimension(:,:), allocatable :: nz8i
    !dir$ attributes align : 64 :: nx8r
    !dir$ attributes align : 64 :: nx8i
    !dir$ attributes align : 64 :: ny8r
    !dir$ attributes align : 64 :: ny8i
    !dir$ attributes align : 64 :: nz8r
    !dir$ attributes align : 64 :: nz8i
    
     
    ! Antenna apperture surface vector function 'N' (theta,phi) coordinate
    ! complex (single)
    complex(kind=sp), dimension(:,:), allocatable :: ntc4
    complex(kind=sp), dimension(:,:), allocatable :: npc4
    !dir$ attributes align : 64 :: ntc4
    !dir$ attributes align : 64 :: npc4
    
    ! Antenna apperture surface vector function 'N' (theta,phi) coordinate
    ! complex (double)
    complex(kind=dp), dimension(:,:), allocatable :: ntc8
    complex(kind=dp), dimension(:,:), allocatable :: npc8
    !dir$ attributes align : 64 :: ntc8
    !dir$ attributes align : 64 :: npc8
    
      
    ! Antenna apperture surface vector function 'N' (theta,phi) coordinate
    ! complex (real-single)
    real(kind=sp), dimension(:,:), allocatable :: ntr4
    real(kind=sp), dimension(:,:), allocatable :: nti4
    real(kind=sp), dimension(:,:), allocatable :: npr4
    real(kind=sp), dimension(:,:), allocatable :: npi4
    !dir$ attributes align : 64 :: ntr4
    !dir$ attributes align : 64 :: nti4
    !dir$ attributes align : 64 :: npr4
    !dir$ attributes align : 64 :: npi4
    
    
    ! Antenna apperture surface vector function 'N' (theta,phi) coordinate
    ! complex (real-double)
    real(kind=dp), dimension(:,:), allocatable :: ntr8
    real(kind=dp), dimension(:,:), allocatable :: nti8
    real(kind=dp), dimension(:,:), allocatable :: npr8
    real(kind=dp), dimension(:,:), allocatable :: npi8
    !dir$ attributes align : 64 :: ntr8
    !dir$ attributes align : 64 :: nti8
    !dir$ attributes align : 64 :: npr8
    !dir$ attributes align : 64 :: npi8
    
    
    ! Complex Electric Current (complex-single)
    ! First dimension nth field, second dimension number of sample points
    complex(kind=sp), dimension(:,:), allocatable :: jexc4
    complex(kind=sp), dimension(:,:), allocatable :: jeyc4
    complex(kind=sp), dimension(:,:), allocatable :: jezc4
    !dir$ attributes align : 64 :: jexc4
    !dir$ attributes align : 64 :: jeyc4
    !dir$ attributes align : 64 :: jezc4
    
    ! Complex Electric Current (complex-double)
     ! First dimension nth field, second dimension number of sample points
    complex(kind=dp), dimension(:,:), allocatable :: jexc8
    complex(kind=dp), dimension(:,:), allocatable :: jeyc8
    complex(kind=dp), dimension(:,:), allocatable :: jezc8
    !dir$ attributes align : 64 :: jexc8
    !dir$ attributes align : 64 :: jeyc8
    !dir$ attributes align : 64 :: jezc8
    
    ! Complex Electric Current (real-single) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=sp), dimension(:,:), allocatable :: jex4r
    real(kind=sp), dimension(:,:), allocatable :: jex4i
    real(kind=sp), dimension(:,:), allocatable :: jey4r
    real(kind=sp), dimension(:,:), allocatable :: jey4i
    real(kind=sp), dimension(:,:), allocatable :: jez4r
    real(kind=sp), dimension(:,:), allocatable :: jez4i
    !dir$ attributes align : 64 :: jex4r
    !dir$ attributes align : 64 :: jex4i
    !dir$ attributes align : 64 :: jey4r
    !dir$ attributes align : 64 :: jey4i
    !dir$ attributes align : 64 :: jez4r
    !dir$ attributes align : 64 :: jez4i

    ! Complex Electric Current (real-double) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=dp), dimension(:,:), allocatable :: jex8r
    real(kind=dp), dimension(:,:), allocatable :: jex8i
    real(kind=dp), dimension(:,:), allocatable :: jey8r
    real(kind=dp), dimension(:,:), allocatable :: jey8i
    real(kind=dp), dimension(:,:), allocatable :: jez8r
    real(kind=dp), dimension(:,:), allocatable :: jez8i
    !dir$ attributes align : 64 :: jex8r
    !dir$ attributes align : 64 :: jex8i
    !dir$ attributes align : 64 :: jey8r
    !dir$ attributes align : 64 :: jey8i
    !dir$ attributes align : 64 :: jez8r
    !dir$ attributes align : 64 :: jez8i
    
    ! Complex Magnetic Current (complex-single)
    ! First dimension nth field, second dimension number of sample points
    complex(kind=sp), dimension(:,:), allocatable :: jhxc4
    complex(kind=sp), dimension(:,:), allocatable :: jhyc4
    complex(kind=sp), dimension(:,:), allocatable :: jhzc4
    !dir$ attributes align : 64 :: jhxc4
    !dir$ attributes align : 64 :: jhyc4
    !dir$ attributes align : 64 :: jhzc4
    
    ! Complex Magnetic Current (complex-double)
     ! First dimension nth field, second dimension number of sample points
    complex(kind=dp), dimension(:,:), allocatable :: jhxc8
    complex(kind=dp), dimension(:,:), allocatable :: jhyc8
    complex(kind=dp), dimension(:,:), allocatable :: jhzc8
    !dir$ attributes align : 64 :: jhxc8
    !dir$ attributes align : 64 :: jhyc8
    !dir$ attributes align : 64 :: jhzc8
    
    ! Complex Magnetic Current (real-single) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=sp), dimension(:,:), allocatable :: jhx4r
    real(kind=sp), dimension(:,:), allocatable :: jhx4i
    real(kind=sp), dimension(:,:), allocatable :: jhy4r
    real(kind=sp), dimension(:,:), allocatable :: jhy4i
    real(kind=sp), dimension(:,:), allocatable :: jhz4r
    real(kind=sp), dimension(:,:), allocatable :: jhz4i
    !dir$ attributes align : 64 :: jhx4r
    !dir$ attributes align : 64 :: jhx4i
    !dir$ attributes align : 64 :: jhy4r
    !dir$ attributes align : 64 :: jhy4i
    !dir$ attributes align : 64 :: jhz4r
    !dir$ attributes align : 64 :: jhz4i

    ! Complex Magnetic Current (real-double) (decomposed)
    ! First dimension nth field, second dimension number of sample points
    real(kind=dp), dimension(:,:), allocatable :: jhx8r
    real(kind=dp), dimension(:,:), allocatable :: jhx8i
    real(kind=dp), dimension(:,:), allocatable :: jhy8r
    real(kind=dp), dimension(:,:), allocatable :: jhy8i
    real(kind=dp), dimension(:,:), allocatable :: jhz8r
    real(kind=dp), dimension(:,:), allocatable :: jhz8i
    !dir$ attributes align : 64 :: jhx8r
    !dir$ attributes align : 64 :: jhx8i
    !dir$ attributes align : 64 :: jhy8r
    !dir$ attributes align : 64 :: jhy8i
    !dir$ attributes align : 64 :: jhz8r
    !dir$ attributes align : 64 :: jhz8i
    
    ! Normalized function Psi (power of radiated field)
    real(kind=sp), dimension(:,:), allocatable :: psir4
    !dir$ attributes align : 64 :: psir4
    ! Normalized function Psi (power of radiated field)
    real(kind=dp), dimension(:,:), allocatable :: psir8
    !dir$ attributes align : 64 :: psir8 
    
    ! An ansemble of Psi functions e.g. (phased arrays)
    real(kind=sp), dimension(:,:,:), allocatable :: psir4a
    !dir$ attributes align : 64 :: psir4a
    real(kind=dp), dimension(:,:,:), allocatable :: psir8a
    !dir$ attributes align : 64 :: psir8a
    
    
    
























end module antenna_model_adt
