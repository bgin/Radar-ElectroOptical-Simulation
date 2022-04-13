
module radar_atmos_loss


 !===================================================================================85
 !---------------------------- DESCRIPTION ------------------------------------------85
 !
 !
 !
 !          Module  name:
 !                         'radar_atmos_loss'
 !          
 !          Purpose:
 !
 !                         Atmospheric loss calculation of radar signal propagation.
 !          History:
 !                        Date: 11-04-2022
 !                        Time: 14:08 GMT+2
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
 !          Re_sperences:
 !         
 !                   (Artech House Radar Library) Sergey A. Leonov, Alexander I. Leonov - Handbook of Computer Simulation in Radio Engineering, Communications and Radar-Artech House (2001)
 !         
 !          E-mail:
 !                  
 !                      beniekg@gmail.com
 !==================================================================================85
    ! Tab:5 col - Type and etc.. definitions
    ! Tab:10,11 col - Type , function and subroutine code blocks.

     use mod_kinds, only : i4,sp,dp
     
     implicit none
     public
     !=====================================================59
     !  File and module in_spormation:
     !  version,creation and build date, author,description
     !=====================================================59

     ! Major version
     integer(kind=i4),  parameter :: RADAR_ATMOS_LOSS_MAJOR = 1
     ! Minor version
     integer(kind=i4),  parameter :: RADAR_ATMOS_LOSS_MINOR = 0
     ! Micro version
     integer(kind=i4),  parameter :: RADAR_ATMOS_LOSS_MICRO = 0
     ! Full version
     integer(kind=i4),  parameter :: RADAR_ATMOS_LOSS_FULLVER =   &
            1000*RADAR_ATMOS_LOSS_MAJOR+100*RADAR_ATMOS_LOSS_MINOR+10*RADAR_ATMOS_LOSS_MICRO
     ! Module creation date
     character(*),        parameter :: RADAR_ATMOS_LOSS_CREATE_DATE = "20-12-2021 15:54 +00200 (MON 20 DEC 2021 GMT+2)"
     ! Module build date
     character(*),        parameter :: RADAR_ATMOS_LOSS_BUILD_DATE  = __DATE__ " " __TIME__
     ! Module author in_spo
     character(*),        parameter :: RADAR_ATMOS_LOSS_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com"
     ! Short description
     character(*),        parameter :: RADAR_ATMOS_LOSS_SYNOPSIS    = "Radar atmospheric loss eqautions implementation"

     !==============================!
     !        Constant data         !
     !==============================!

     ! Oxygen resonance _sprequencies (Blake Model usage) 
     real(sp), dimension(48), parameter, private :: f_N_plus = [6.2648_sp,56.2648_sp,58.4466_sp,58.4466_sp,   &
	                                                       59.5910_sp,59.5910_sp,60.4348_sp,60.4348_sp,  &
							       61.1506_sp,61.1506_sp,61.8002_sp,61.8002_sp,  &
							       62.4212_sp,62.4212_sp,62.9980_sp,62.9980_sp,  &
							       63.5658_sp,63.5658_sp,64.1272_sp,64.1272_sp,  &
							       64.6779_sp,64.6779_sp,65.2240_sp,65.2240_sp,  &
							       65.7626_sp,65.7626_sp,66.2978_sp,66.2978_sp,  &
							       66.8313_sp,66.8313_sp,67.3627_sp,67.3627_sp,  &
							       67.8923_sp,67.8923_sp,68.4205_sp,68.4205_sp,  &
							       68.9478_sp,68.9478_sp,69.4741_sp,69.4741_sp,  &
							       70.0_sp,70.0_sp,70.5249_sp,70.55249_sp,       &
							       71.0497_sp,71.0497_sp,0.0_sp,0.0_sp,0.0_sp]

     real(sp), dimension(48), parameter, private :: f_N_minus = [118.7505_sp,118.7505_sp,62.4862_sp,62.4862_sp,  &
	                                                        60.3061_sp,60.3061_sp,59.1642_sp,59.1642_sp,    &
							        58.3239_sp,58.3239_sp,57.6125_sp,57.6125_sp,    &
							        56.9682_sp,56.9682_sp,56.3634_sp,56.3634_sp,    &
							        55.7839_sp,55.7839_sp,55.2214_sp,55.2214_sp,    &
							        54.6728_sp,54.6728_sp,54.1294_sp,54.1294_sp,    &
							        53.5960_sp,53.5960_sp,53.0695_sp,53.0695_sp,    &
							        52.5458_sp,52.5458_sp,52.0259_sp,52.0259_sp,    &
							        51.5091_sp,51.5091_sp,50.9949_sp,50.9949_sp,    &
							        50.4830_sp,50.4830_sp,49.9730_sp,49.9730_sp,    &
							        49.4648_sp,49.4648_sp,48.9582_sp,48.9582_sp,    &
							        48.4530_sp,48.4530_sp,0.0_sp,0.0_sp,0.0_sp]

     real(sp), dimension(48), parameter, private :: Z        =  [1.0_sp,0.0_sp,1.0_sp,0.0_sp,1.0_sp,0.0_sp,1.0_sp, &
                                                                0.0_sp,1.0_sp,0.0_sp,1.0_sp,0.0_sp,1.0_sp,0.0_sp, &
                                                                1.0_sp,0.0_sp,1.0_sp,0.0_sp,1.0_sp,0.0_sp,1.0_sp, &
                                                                0.0_sp,1.0_sp,0.0_sp,1.0_sp,0.0_sp,1.0_sp,0.0_sp, &
                                                                1.0_sp,0.0_sp,1.0_sp,0.0_sp,1.0_sp,0.0_sp,1.0_sp, &
                                                                0.0_sp,1.0_sp,0.0_sp,1.0_sp,0.0_sp,1.0_sp,0.0_sp, &
                                                                1.0_sp,0.0_sp,1.0_sp,0.0_sp,0.0_sp,0.0_sp]


      real(dp), dimension(48), parameter, private :: f_N_plusr8 = [6.2648_dp,56.2648_dp,58.4466_dp,58.4466_dp,   &
	                                                          59.5910_dp,59.5910_dp,60.4348_dp,60.4348_dp,  &
							          61.1506_dp,61.1506_dp,61.8002_dp,61.8002_dp,  &
							          62.4212_dp,62.4212_dp,62.9980_dp,62.9980_dp,  &
							          63.5658_dp,63.5658_dp,64.1272_dp,64.1272_dp,  &
							          64.6779_dp,64.6779_dp,65.2240_dp,65.2240_dp,  &
							          65.7626_dp,65.7626_dp,66.2978_dp,66.2978_dp,  &
							          66.8313_dp,66.8313_dp,67.3627_dp,67.3627_dp,  &
							          67.8923_dp,67.8923_dp,68.4205_dp,68.4205_dp,  &
							          68.9478_dp,68.9478_dp,69.4741_dp,69.4741_dp,  &
							          70.0_dp,70.0_dp,70.5249_dp,70.55249_dp,       &
							          71.0497_dp,71.0497_dp,0.0_dp,0.0_dp,0.0_dp]

     real(dp), dimension(48), parameter, private :: f_N_minusr8 = [118.7505_dp,118.7505_dp,62.4862_dp,62.4862_dp,  &
	                                                          60.3061_dp,60.3061_dp,59.1642_dp,59.1642_dp,    &
							          58.3239_dp,58.3239_dp,57.6125_dp,57.6125_dp,    &
							          56.9682_dp,56.9682_dp,56.3634_dp,56.3634_dp,    &
							          55.7839_dp,55.7839_dp,55.2214_dp,55.2214_dp,    &
							          54.6728_dp,54.6728_dp,54.1294_dp,54.1294_dp,    &
							          53.5960_dp,53.5960_dp,53.0695_dp,53.0695_dp,    &
							          52.5458_dp,52.5458_dp,52.0259_dp,52.0259_dp,    &
							          51.5091_dp,51.5091_dp,50.9949_dp,50.9949_dp,    &
							          50.4830_dp,50.4830_dp,49.9730_dp,49.9730_dp,    &
							          49.4648_dp,49.4648_dp,48.9582_dp,48.9582_dp,    &
							          48.4530_dp,48.4530_dp,0.0_dp,0.0_dp,0.0_dp]

     real(dp), dimension(48), parameter, private :: Zr8        =  [1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp, &
                                                                  0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp, &
                                                                  1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp, &
                                                                  0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp, &
                                                                  1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp, &
                                                                  0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp, &
                                                                  1.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]


     

   contains

     

    subroutine blake_atmos_loss_r4_1(h_a,_sp,theta,R_m,K,L1a)
        !dir$ attributes align : 32 :: blake_atmos_loss_r4_1
        !dir$ optimize : 3 
        !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: blake_atmos_loss_r4_1
        real(sp),                intent(in)  :: h_a   ! m, target heigth
        real(sp),                intent(in)  :: _sp     ! Mhz, radar _sprequency
        real(sp),                intent(in)  :: theta ! deg, angle
        real(sp),                intent(in)  :: R_m   ! m, target range
        integer(i4),             intent(in)  :: K     ! Number o_sp integration points (steps)
        real(sp),                intent(out) :: L1a   ! dB, signal strength loss
        ! Locals
        real(sp), dimension(48), parameter :: u_plus = [ 2.5_sp,       &
			                                 4.6666665_sp, &
	                                                 6.7500000_sp, &
						         8.8000002_sp, &
			                                 10.8333330_sp,&
							 12.8571424_sp,&
							 14.8750000_sp,&
							 16.8888893_sp,&
							 18.8999996_sp,&
							 20.9090900_sp,&
							 22.9166660_sp,&
							 24.9230766_sp,&
							 26.9285717_sp,&
							 28.9333324_sp,&
							 30.9375000_sp,&
							 32.9411774_sp,&
							 34.9444427_sp,&
							 36.9473686_sp,&
							 38.9500008_sp,&
							 40.9523811_sp,&
                                                         42.9545441_sp,&
                                                         44.9565201_sp,&
                                                         46.9583321_sp,&
                                                         48.9599991_sp,&
                                                         50.9615402_sp,&
                                                         52.9629631_sp,&
                                                         54.9642868_sp,&
                                                         56.9655190_sp,&
                                                         58.9666672_sp,&
                                                         60.9677429_sp,&
                                                         62.9687500_sp,&
                                                         64.9696960_sp,&
                                                         66.9705887_sp,&
                                                         68.9714279_sp,&
                                                         70.9722214_sp,&
                                                         72.9729767_sp,&
                                                         74.9736862_sp,&
                                                         76.9743576_sp,&
                                                         78.9749985_sp,&
                                                         80.9756088_sp,&
                                                         82.9761887_sp,&
                                                         84.9767456_sp,&
                                                         86.9772720_sp,&
                                                         88.9777756_sp,&
                                                         90.9782639_sp,&
                                                         0.0_sp,0.0_sp,0.0_sp]
         real(sp), dimension(48), parameter :: u_minus = [2.0000000_sp, &
                                                          4.5000000_sp, &
                                                          6.6666665_sp, &
                                                          8.7500000_sp, &
                                                          10.8000002_sp,&
                                                          12.8333330_sp,&
                                                          14.8571424_sp,&
                                                          16.8750000_sp,&
                                                          18.8888893_sp,&
                                                          20.8999996_sp,&
                                                          22.9090900_sp,&
                                                          24.9166660_sp,&
                                                          26.9230766_sp,&
                                                          28.9285717_sp,&
                                                          30.9333324_sp,&
                                                          32.9375000_sp,&
                                                          34.9411774_sp,&
                                                          36.9444427_sp,&
                                                          38.9473686_sp,&
                                                          40.9500008_sp,&
                                                          42.9523811_sp,&
                                                          44.9545441_sp,&
                                                          46.9565201_sp,&
                                                          48.9583321_sp,&
                                                          50.9599991_sp,&
                                                          52.9615402_sp,&
                                                          54.9629631_sp,&
                                                          56.9642868_sp,&
                                                          58.9655190_sp,&
                                                          60.9666672_sp,&
                                                          62.9677429_sp,&
                                                          64.9687500_sp,&
                                                          66.9696960_sp,&
                                                          68.9705887_sp,&
                                                          70.9714279_sp,&
                                                          72.9722214_sp,&
                                                          74.9729767_sp,&
                                                          76.9736862_sp,&
                                                          78.9743576_sp,&
                                                          80.9749985_sp,&
                                                          82.9756088_sp,&
                                                          84.9761887_sp,&
                                                          86.9767456_sp,&
                                                          88.9772720_sp,&
                                                          90.9777756_sp,&
							  0.0_sp,0.0_sp,0.0_sp]
       real(sp), dimension(48), parameter :: u_0 = [9.0000000_sp, &
                                                    11.6666670_sp,&
                                                    15.1666670_sp,&
                                                    18.8999996_sp,&
                                                    22.7333336_sp,&
                                                    26.6190472_sp,&
                                                    30.5357151_sp,&
                                                    34.4722214_sp,&
                                                    38.4222221_sp,&
                                                    42.3818169_sp,&
                                                    46.3484840_sp,&
                                                    50.3205147_sp,&
                                                    54.2967033_sp,&
                                                    58.2761917_sp,&
                                                    62.2583351_sp,&
                                                    66.2426453_sp,&
                                                    70.2287598_sp,&
                                                    74.2163773_sp,&
                                                    78.2052612_sp,&
                                                    82.1952362_sp,&
                                                    86.1861496_sp,&
                                                    90.1778641_sp,&
                                                    94.1702881_sp,&
                                                    98.1633301_sp,&
                                                    102.1569214_sp,&
                                                    106.1509933_sp,&
                                                    110.1455002_sp,&
                                                    114.1403961_sp,&
                                                    118.1356354_sp,&
                                                    122.1311798_sp,&
                                                    126.1270142_sp,&
                                                    130.1231079_sp,&
                                                    134.1194305_sp,&
                                                    138.1159668_sp,&
                                                    142.1127014_sp,&
                                                    146.1096039_sp,&
                                                    150.1066895_sp,&
                                                    154.1039124_sp,&
                                                    158.1012878_sp,&
                                                    162.0987854_sp,&
                                                    166.0964050_sp,&
                                                    170.0941315_sp,&
                                                    174.0919647_sp,&
                                                    178.0899048_sp,&
                                                    182.0879211_sp,&
                                                    0.0_sp,0.0_sp,0.0_sp]
        real(sp), dimension(48), parameter :: ser_n = [2.0000000_sp, &
                                                       6.0000000_sp, &
                                                       12.0000000_sp,&
                                                       20.0000000_sp,&
                                                       30.0000000_sp,&
                                                       42.0000000_sp,&
                                                       56.0000000_sp,&
                                                       72.0000000_sp,&
                                                       90.0000000_sp,&
                                                       110.0000000_sp,&
                                                       132.0000000_sp,&
                                                       156.0000000_sp,&
                                                       182.0000000_sp,&
                                                       210.0000000_sp,&
                                                       240.0000000_sp,&
                                                       272.0000000_sp,&
                                                       306.0000000_sp,&
                                                       342.0000000_sp,&
                                                       380.0000000_sp,&
                                                       420.0000000_sp,&
                                                       462.0000000_sp,&
                                                       506.0000000_sp,&
                                                       552.0000000_sp,&
                                                       600.0000000_sp,&
                                                       650.0000000_sp,&
                                                       702.0000000_sp,&
                                                       756.0000000_sp,&
                                                       812.0000000_sp,&
                                                       870.0000000_sp,&
                                                       930.0000000_sp,&
                                                       992.0000000_sp,&
                                                       1056.0000000_sp,&
                                                       1122.0000000_sp,&
                                                       1190.0000000_sp,&
                                                       1260.0000000_sp,&
                                                       1332.0000000_sp,&
                                                       1406.0000000_sp,&
                                                       1482.0000000_sp,&
                                                       1560.0000000_sp,&
                                                       1640.0000000_sp,&
                                                       1722.0000000_sp,&
                                                       1806.0000000_sp,&
                                                       1892.0000000_sp,&
                                                       1980.0000000_sp,&
                                                       2070.0000000_sp,&
                                                       0.0_sp,0.0_sp,0.0_sp]
        
        real(sp), parameter :: N    = 0.000313_sp
        real(sp), parameter :: c_e  = -0.149_sp     !1/km a decay constant
        real(sp), parameter :: n0   = 1.000313_sp   !refractive index
        real(sp), parameter :: r_km = 6370.0_sp     !Earth radius
        real(sp), parameter :: r_m  = 6370000.0_sp  !Earth radius
        real(sp), parameter :: a_e  = 8493333.0_sp  !Earth effective radius
	real(sp), parameter :: p0   = 1013.25_sp    !atmospheric pressure constant
        real(sp), parameter :: alf1 = 5.2561222_sp   !tropospheric model constants
	real(sp), parameter :: alf2 = 0.034164794_sp  !as above
	real(sp), parameter :: alf3 = 11.388265_sp    !as above
	real(sp), parameter :: T0   = 300.0_sp      !standard temperature, K
	real(sp), parameter :: C    = 2.0058_sp     !absorption coeff const
	real(sp), parameter :: z    = 0.017453292519943295769236907685_sp !deg-to-rad (PI/180)
	real(sp), parameter :: zth  = z*theta
	real(sp), parameter :: f_ghz= f*1000.0_sp
	real(sp), parameter :: fghz2= f_ghz*f_ghz
	real(sp), parameter :: czth = R_m*cos(zth)
	real(sp), parameter :: h_m  = R_m*sin(zth)+((czth*czth)/(2.0_sp*a_e))
	real(sp), parameter :: h_km = h_m/1000.0_sp
        real(sp), parameter :: delh = h_km/real(K,kind=sp)
        integer(i4), parameter :: n45 = 45 ! inner loop iterations
	real(sp), automatic :: T_k  = 0.0_sp  !//K, atmos temperature
	real(sp), automatic :: P_k  = 0.0_sp  !//mb, atmos pressure
	real(sp), automatic :: g_k  = 0.0_sp  !// line breadth constant parameter
        real(sp), automatic :: gamma_k = 0.0_sp !// dB/km, absorption coefficient
	real(sp), automatic :: S1   = 0.0_sp  !// absorption loss integral coeff
	real(sp), automatic :: S2   = 0.0_sp  !// absorption loss integral coeff
	real(sp), automatic :: S3   = 0.0_sp   !// absorption loss integral coeff
	real(sp), automatic :: m_k  = 0.0_sp   !// refractive index model
        real(sp), automatic :: gamma0  = 0.0_sp
	real(sp), automatic :: h0   = 0.0_sp
        real(sp), automatic :: m0   = 0.0_sp
        !
	real(sp), volatile  :: u_plus_preload  = u_plus(1)
	real(sp), volatile  :: u_minus_preload = u_minus(1)
	real(sp), volatile  :: u_0_preload     = u_0(1)
        real(sp), volatile  :: ser_n_preload   = ser_n(1)
        real(sp), automatic :: tK,h_k,h_gm,h_gkm,delfk,F0_k,delfk2, &
                               Sigma_K,t0,Sig1_plus,Sig2_plus,      &
                               t1,Sig1_minus,Sig2_minus,F_n_plus,   &
                               F_n_minus,t2,A1,A2,A3,En,T_k3,n0czth,&
                               cterm,term1,term12
        integer(i4) :: i,j
        
        if(K<=1) then
           L1a = 0.0_sp
           return
        end if

        do i = 1, K
             tK    = real(i,kind=sp)
             h_k   = h_a*0.001_sp+k*delh !//km, current height
             h_gm  = r_m*h_k*1000.0_sp/(r_m+h_k*1000.0_sp) !//m, geopotential altitude
             ! Atmosphere temp
             if(h_gkm<=11.0_sp) then
                 T_k = 288.16_sp-0.0065_sp*h_k*1000.0_sp
             else if(h_gkm>11.0_sp .and. h_gkm<25.0_sp) then
                 T_k = 216.66_sp 
             else
                 T_k = 216.66_sp+0.003_sp*(h_k-25.0_sp)*1000.0_sp
             end if

             if(h_gkm<=11.0_sp) then
                 P_k = p0*(T_k*0.003470294280955024986118822876_sp)**alf1
             else if(h_gkm>11.0_sp .and. h_gkm<25.0_sp) then
                 t0 = 226.32_sp/T_k
                 P_k = t0*exp(-alf2*(h_k-11.0_sp)*1000.0_sp)
             else
                 P_k = 24.886_sp*(216.66_sp/T_k)**alf3
             end if

             if(h_k<=8.0_sp) then
                 g_k = 0.640_sp
             else if(h_k>8.0_sp .and. h_k<=25.0_sp) then
                 g_k = 0.640_sp+0.04218_sp*(h_k-8.0_sp)
	     else		           
                 g_k = 1.357_sp
             end if

             delfk = g_k*(P_k/p0)*(T0/T_k)  !// line-breadth constant
	     F0_k  = delfk/((delfk*delfk)+fghz2) !// nonresonant contribution
             delfk2 = delfk*delfk
             Sigma_K      = 0.0_sp
             do j = 1, n45
                  t0  = f_N_plus(j)
                  Sig1_plus  = delfk/((t0-f_ghz)*(t0-f_ghz))+delfk2
                  Sig2_plus  = delfk/((t0+f_ghz)*(t0+f_ghz))+delfk2
		  t1         = f_N_minus(j)
		  Sig1_minus = delfk/((t1-f_ghz)*(t1-f_ghz))+delfk2
		  Sig2_minus = delfk/((t1+f_ghz)*(t1+f_ghz))+delfk2
		  F_N_plus   = Sig1_plus+Sig2_plus
		  F_N_minus  = Sig1_minus+Sig2_minus
		  t2         = u_0(j)
		  A1         = u_plus(j)*F_N_plus
		  A2         = u_minus(j)*F_N_minus
		  En         = 2.06844_sp*ser_n(j)
		  Sigma_k    = Sigma_k + Z(j)*((A1+A2+t2*F0_k)*exp(-(En/T_k)))
             end do
             T_k3   = 1.0_sp/(T_k*T_k*T_k)
	     if(i==0) then
                gamma0          = C*P_k*T_k3*fghz2*Sigma_k
		m0              = 1.0_sp+N*exp(c_e*h_k)
		h0              = h_a*0.001_sp+k*delh
	     end if
	     gamma_k            = C*P_k*T_k3*fghz2*Sigma_k
	     m_k                = 1.0_sp+N*exp(c_e*h_k)
	      
	     n0czth = n0*czth
	     term1  = n0czth/(m_k*(1.0_sp+h_k/r_km))
	     term12 = term1*term1
             S2                 = gamma_k/sqrt(1.0f-term12)
	     S3                 S3 + S2
         end do
         
	term1 = n0czth/(m0*(1.0_sp+h0/r_km))
	term12= term1*term1
	S1 = gamma0/sqrt(1.0_sp-term12)
        L1a = 2.0_sp*(S1+S2+2.0_sp*S3)*delh*0.5f
    end subroutine blake_atmos_loss_r4_1
     
     

     


















end module radar_atmos_loss
