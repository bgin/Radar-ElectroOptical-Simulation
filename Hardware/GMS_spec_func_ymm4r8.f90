

#include "GMS_config.fpp"

!/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!*/

module spec_funcs_ymm4r8


!===================================================================================85
 !---------------------------- DESCRIPTION ------------------------------------------85
 !
 !
 !
 !          Module  name:
 !                         spec_funcs_ymm4r8
 !          
 !          Purpose:
 !                       Various vectorized special functions.
 !                        
 !          History:
 !                        Date: 08-28-2023
 !                        Time: 16:09 GMT+2
 !                        
 !          Version:
 !
 !                      Major: 1
 !                      Minor: 0
 !                      Micro: 0
 !
 !          Author:  
 !                      Vectorized by Bernard Gingold (based on different authors work)
 !                      The details stated by the specific function description.
 !                 
 !          References:
 !         
 !                      Provided at the specific function description
 !         
 !          E-mail:
 !                  
 !                      beniekg@gmail.com
!==================================================================================85
    ! Tab:5 col - Type and etc.. definitions
    ! Tab:10,11 col - Type , function and subroutine code blocks.

    use mod_kinds,    only : i4,dp
    use mod_vectypes, only : YMM4r8_t,Mask4_t
    
    public
    implicit none
    
    
      ! Major version
     integer(kind=i4),  parameter :: SPEC_FUNCS_YMM4R8_MAJOR = 1
     ! Minor version
     integer(kind=i4),  parameter :: SPEC_FUNCS_YMM4R8_MINOR = 0
     ! Micro version
     integer(kind=i4),  parameter :: SPEC_FUNCS_YMMY4R8_MICRO = 0
     ! Full version
     integer(kind=i4),  parameter :: SPEC_FUNCS_YMMY4R8_FULLVER =   &
            1000*SPEC_FUNCS_YMM4R8_MAJOR+100*SPEC_FUNCS_YMM4R8_MINOR+10*SPEC_FUNCS_YMM4R8_MICRO
     ! Module creation date
     character(*),        parameter :: SPEC_FUNCS_YMM4R8_CREATE_DATE = "28-08-2022 06:11 +00200 (MON 28 AUG 2023 GMT+2)"
     ! Module build date
     character(*),        parameter :: SPEC_FUNCS_YMM4R8_BUILD_DATE  = __DATE__ " " __TIME__
     ! Module author info
     character(*),        parameter :: SPEC_FUNCS_YMMY4R8_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com"
     ! Short description
     character(*),        parameter :: SPEC_FUNCS_YMM4R8_SYNOPSIS    = "Vectorized various special functions" 
     
     
     !! calcei_zm8r8 constants arrays (saved).
     
     !dir$ attributes align : 64 :: calcei_a
     !dir$ attributes align : 64 :: calcei_b
     !dir$ attributes align : 64 :: calcei_c
     !dir$ attributes align : 64 :: calcei_d
     !dir$ attributes align : 64 :: calcei_e
     !dir$ attributes align : 64 :: calcei_f
     !dir$ attributes align : 64 :: calcei_plg
     !dir$ attributes align : 64 :: calcei_qlg
     !dir$ attributes align : 64 :: calcei_p
     !dir$ attributes align : 64 :: calcei_q
     !dir$ attributes align : 64 :: calcei_r
     !dir$ attributes align : 64 :: calcei_s
     !dir$ attributes align : 64 :: calcei_p1
     !dir$ attributes align : 64 :: calcei_q1
     !dir$ attributes align : 64 :: calcei_p2
     !dir$ attributes align : 64 :: calcei_q2
     type(YMM4r8_t), dimension(0:6), save :: calcei_a = [YMM4r8_t(1.1669552669734461083368e+2_dp),   &
	                                                  YMM4r8_t(2.1500672908092918123209e+3_dp),   &
                                                          YMM4r8_t(1.5924175980637303639884e+4_dp),   &
                                                          YMM4r8_t(8.9904972007457256553251e+4_dp),   &
                                                          YMM4r8_t(1.5026059476436982420737e+5_dp),   &
                                                          YMM4r8_t(-1.4815102102575750838086e+5_dp),  &
                                                          YMM4r8_t(5.0196785185439843791020e+0_dp)]
      type(YMM4r8_t), dimension(0:6), save :: calcei_b = [YMM4r8_t(4.0205465640027706061433e+1_dp),   & 
                                                          YMM4r8_t(7.5043163907103936624165e+2_dp),   &
                                                          YMM4r8_t(8.1258035174768735759855e+3_dp),   & 
                                                          YMM4r8_t(5.2440529172056355429883e+4_dp),   & 
                                                          YMM4r8_t(1.8434070063353677359298e+5_dp),   & 
                                                          YMM4r8_t(2.5666493484897117319268e+5_dp)]
      type(YMM4r8_t), dimension(0:8), save :: calcei_c = [YMM4r8_t(3.828573121022477169108e-1_dp),    & 
                                                          YMM4r8_t(1.107326627786831743809e+1_dp),    &
                                                          YMM4r8_t(7.246689782858597021199e+1_dp),    & 
                                                          YMM4r8_t(1.700632978311516129328e+2_dp),    & 
                                                          YMM4r8_t(1.698106763764238382705e+2_dp),    &
                                                          YMM4r8_t(7.633628843705946890896e+1_dp),    & 
                                                          YMM4r8_t(1.487967702840464066613e+1_dp),    & 
                                                          YMM4r8_t(9.999989642347613068437e-1_dp),    & 
                                                          YMM4r8_t(1.737331760720576030932e-8_dp)]
      type(YMM4r8_t), dimension(0:8), save :: calcei_d = [YMM4r8_t(8.258160008564488034698e-2_dp),    & 
	                                                  YMM4r8_t(4.344836335509282083360e+0_dp),    & 
                                                          YMM4r8_t(4.662179610356861756812e+1_dp),    & 
                                                          YMM4r8_t(1.775728186717289799677e+2_dp),    & 
                                                          YMM4r8_t(2.953136335677908517423e+2_dp),    & 
                                                          YMM4r8_t(2.342573504717625153053e+2_dp),    & 
                                                          YMM4r8_t(9.021658450529372642314e+1_dp),    & 
                                                          YMM4r8_t(1.587964570758947927903e+1_dp),    & 
                                                          YMM4r8_t(1.000000000000000000000e+0_dp)]
      type(YMM4r8_t), dimension(0:9), save :: calcei_e = [YMM4r8_t(1.3276881505637444622987e+2_dp),   &
                                                          YMM4r8_t(3.5846198743996904308695e+4_dp),   &
                                                          YMM4r8_t(1.7283375773777593926828e+5_dp),   &
                                                          YMM4r8_t(2.6181454937205639647381e+5_dp),   &
                                                          YMM4r8_t(1.7503273087497081314708e+5_dp),   & 
                                                          YMM4r8_t(5.9346841538837119172356e+4_dp),   &
                                                          YMM4r8_t(1.0816852399095915622498e+4_dp),   &
                                                          YMM4r8_t(1.0611777263550331766871e+03_dp),  &
                                                          YMM4r8_t(5.2199632588522572481039e+1_dp),   &
                                                          YMM4r8_t(9.9999999999999999087819e-1_dp)]
      type(YMM4r8_t), dimension(0:9), save :: calcei_f = [YMM4r8_t(3.9147856245556345627078e+4_dp),   &
                                                          YMM4r8_t(2.5989762083608489777411e+5_dp),   &
                                                          YMM4r8_t(5.5903756210022864003380e+5_dp),   &
                                                          YMM4r8_t(5.4616842050691155735758e+5_dp),   &
                                                          YMM4r8_t(2.7858134710520842139357e+5_dp),   &
                                                          YMM4r8_t(7.9231787945279043698718e+4_dp),   &
                                                          YMM4r8_t(1.2842808586627297365998e+4_dp),   &
                                                          YMM4r8_t(1.1635769915320848035459e+3_dp),   &
                                                          YMM4r8_t(5.4199632588522559414924e+1_dp),   &
                                                          YMM4r8_t(1.0000000000000000000000e+0_dp)]
      type(YMM4r8_t), dimension(0:3), save :: calcei_plg=[YMM4r8_t(-2.4562334077563243311e+01_dp),    &
                                                          YMM4r8_t(2.3642701335621505212e+02_dp),     &
                                                          YMM4r8_t(-5.4989956895857911039e+02_dp),    &
                                                          YMM4r8_t(3.5687548468071500413e+02_dp)]
      type(YMM4r8_t), dimension(0:3), save :: calcei_qlg=[YMM4r8_t(-3.5553900764052419184e+01_dp),    &
                                                          YMM4r8_t(1.9400230218539473193e+02_dp),     &
                                                          YMM4r8_t(-3.3442903192607538956e+02_dp),    &
                                                          YMM4r8_t(1.7843774234035750207e+02_dp)]
      type(YMM4r8_t), dimension(0:9), save :: calcei_p  =[YMM4r8_t(-1.2963702602474830028590e+01_dp), &
                                                          YMM4r8_t(-1.2831220659262000678155e+03_dp), &
                                                          YMM4r8_t(-1.4287072500197005777376e+04_dp), &
                                                          YMM4r8_t(-1.4299841572091610380064e+06_dp), &
                                                          YMM4r8_t(-3.1398660864247265862050e+05_dp), &
                                                          YMM4r8_t(-3.5377809694431133484800e+08_dp), &
                                                          YMM4r8_t(3.1984354235237738511048e+08_dp),  &
                                                          YMM4r8_t(-2.5301823984599019348858e+10_dp), &
                                                          YMM4r8_t(1.2177698136199594677580e+10_dp),  &
                                                          YMM4r8_t(-2.0829040666802497120940e+11_dp)]
      type(YMM4r8_t), dimension(0:9), save :: calcei_q  =[YMM4r8_t(7.6886718750000000000000e+01_dp),  &
                                                          YMM4r8_t(-5.5648470543369082846819e+03_dp), &
                                                          YMM4r8_t(1.9418469440759880361415e+05_dp),  &
                                                          YMM4r8_t(-4.2648434812177161405483e+06_dp), &
                                                          YMM4r8_t(6.4698830956576428587653e+07_dp),  &
                                                          YMM4r8_t(-7.0108568774215954065376e+08_dp), &
                                                          YMM4r8_t(5.4229617984472955011862e+09_dp),  &
                                                          YMM4r8_t(-2.8986272696554495342658e+10_dp), &
                                                          YMM4r8_t(9.8900934262481749439886e+10_dp),  &
                                                          YMM4r8_t(-8.9673749185755048616855e+10_dp)]
      type(YMM4r8_t), dimension(0:9), save :: calcei_r  =[YMM4r8_t(-2.645677793077147237806e+00_dp),  &
                                                          YMM4r8_t(-2.378372882815725244124e+00_dp),  &
                                                          YMM4r8_t(-2.421106956980653511550e+01_dp),  & 
                                                          YMM4r8_t(1.052976392459015155422e+01_dp),   &
                                                          YMM4r8_t(1.945603779539281810439e+01_dp),   &
                                                          YMM4r8_t(-3.015761863840593359165e+01_dp),  &
                                                          YMM4r8_t(1.120011024227297451523e+01_dp),   &
                                                          YMM4r8_t(-3.988850730390541057912e+00_dp),  &
                                                          YMM4r8_t(9.565134591978630774217e+00_dp),   & 
                                                          YMM4r8_t(9.981193787537396413219e-1_dp)]
      type(YMM4r8_t), dimension(0:8), save :: calcei_s  =[YMM4r8_t(1.598517957704779356479e-4_dp),    &
                                                          YMM4r8_t(4.644185932583286942650e+00_dp),   &
                                                          YMM4r8_t(3.697412299772985940785e+02_dp),   &
                                                          YMM4r8_t(-8.791401054875438925029e+00_dp),  &
                                                          YMM4r8_t(7.608194509086645763123e+02_dp),   &
                                                          YMM4r8_t(2.852397548119248700147e+01_dp),   &
                                                          YMM4r8_t(4.731097187816050252967e+02_dp),   &
                                                          YMM4r8_t(-2.369210235636181001661e+02_dp),  &
                                                          YMM4r8_t(1.249884822712447891440e+00_dp)]
      type(YMM4r8_t), dimension(0:9), save :: calcei_p1 =[YMM4r8_t(-1.647721172463463140042e+00_dp),  &
                                                          YMM4r8_t(-1.860092121726437582253e+01_dp),  &
                                                          YMM4r8_t(-1.000641913989284829961e+01_dp),  &
                                                          YMM4r8_t(-2.105740799548040450394e+01_dp),  &
                                                          YMM4r8_t(-9.13483569999874255243e-1_dp),    &
                                                          YMM4r8_t(-3.323612579343962284333e+01_dp),  &
                                                          YMM4r8_t(2.495487730402059440626e+01_dp),   &
                                                          YMM4r8_t(2.652575818452799819855e+01_dp),   &
                                                          YMM4r8_t(-1.845086232391278674524e+00_dp),  &
                                                          YMM4r8_t(9.999933106160568739091e-1_dp)]
      type(YMM4r8_t), dimension(0:8), save :: calcei_q1 =[YMM4r8_t(9.792403599217290296840e+01_dp),   &
                                                          YMM4r8_t(6.403800405352415551324e+01_dp),   &
                                                          YMM4r8_t(5.994932325667407355255e+01_dp),   &
                                                          YMM4r8_t(2.538819315630708031713e+02_dp),   &
                                                          YMM4r8_t(4.429413178337928401161e+01_dp),   &
                                                          YMM4r8_t(1.192832423968601006985e+03_dp),   &
                                                          YMM4r8_t(1.991004470817742470726e+02_dp),   &
                                                          YMM4r8_t(-1.093556195391091143924e+01_dp),  &
                                                          YMM4r8_t(1.001533852045342697818e+00_dp)]
      type(YMM4r8_t), dimension(0:9), save :: calcei_p2 =[YMM4r8_t(1.75338801265465972390e+02_dp),    &
                                                          YMM4r8_t(-2.23127670777632409550e+02_dp),   &
                                                          YMM4r8_t(-1.81949664929868906455e+01_dp),   &
                                                          YMM4r8_t(-2.79798528624305389340e+01_dp),   &
                                                          YMM4r8_t(-7.63147701620253630855e+00_dp),   &
                                                          YMM4r8_t(-1.52856623636929636839e+01_dp),   &
                                                          YMM4r8_t(-7.06810977895029358836e+00_dp),   &
                                                          YMM4r8_t(-5.00006640413131002475e+00_dp),   &
                                                          YMM4r8_t(-3.00000000320981265753e+00_dp),   &
                                                          YMM4r8_t(1.00000000000000485503e+00_dp)]
      type(YMM4r8_t), dimension(0:8), save :: calcei_q2 =[YMM4r8_t(3.97845977167414720840e+04_dp),    &
                                                          YMM4r8_t(3.97277109100414518365e+00_dp),    &
                                                          YMM4r8_t(1.37790390235747998793e+02_dp),    &
                                                          YMM4r8_t(1.17179220502086455287e+02_dp),    &
                                                          YMM4r8_t(7.04831847180424675988e+01_dp),    &
                                                          YMM4r8_t(-1.20187763547154743238e+01_dp),   &
                                                          YMM4r8_t(-7.99243595776339741065e+00_dp),   &
                                                          YMM4r8_t(-2.99999894040324959612e+00_dp),   &
                                                          YMM4r8_t(1.99999999999048104167e+00_dp)]
                                                          
                                                          
     !!
     !! calci0_ymm4r8 constant arrays (saved)
     !!
     !dir$ attributes align : 64 :: calci0_p
     !dir$ attributes align : 64 :: calci0_q
     !dir$ attributes align : 64 :: calci0_pp
     !dir$ attributes align : 64 :: calci0_qq
     type(YMM4r8_t), dimension(0:14), save :: calci0_p  =[YMM4r8_t(-5.2487866627945699800e-18_dp),        &
                                                          YMM4r8_t(-1.5982226675653184646e-14_dp),        &
                                                          YMM4r8_t(-2.6843448573468483278e-11_dp),        &
                                                          YMM4r8_t(-3.0517226450451067446e-08_dp),        &
                                                          YMM4r8_t(-2.5172644670688975051e-05_dp),        &
                                                          YMM4r8_t(-1.5453977791786851041e-02_dp),        &
                                                          YMM4r8_t(-7.0935347449210549190e+00_dp),        &
                                                          YMM4r8_t(-2.4125195876041896775e+03_dp),        &
                                                          YMM4r8_t(-5.9545626019847898221e+05_dp),        &
                                                          YMM4r8_t(-1.0313066708737980747e+08_dp),        &
                                                          YMM4r8_t(-1.1912746104985237192e+10_dp),        &
                                                          YMM4r8_t(-8.4925101247114157499e+11_dp),        &
                                                          YMM4r8_t(-3.2940087627407749166e+13_dp),        &
                                                          YMM4r8_t(-5.5050369673018427753e+14_dp),        &
                                                          YMM4r8_t(-2.2335582639474375249e+15_dp)]
      type(YMM4r8_t), dimension(0:4), save ::  calci0_q =[YMM4r8_t(3.7277560179962773046e+03_dp),         &
                                                          YMM4r8_t(6.5158506418655165707e+06_dp),         &
                                                          YMM4r8_t(-6.5626560740833869295e+09_dp),        &
                                                          YMM4r8_t(3.7604188704092954661e+12_dp),         &
                                                          YMM4r8_t(-9.7087946179594019126d+14_dp)]
      type(YMM4r8_t), dimension(0:7), save ::  calci0_pp=[YMM4r8_t(3.9843750000000000000e-01_dp),         & 
                                                          YMM4r8_t(2.9205384596336793945e+00_dp),         &
                                                          YMM4r8_t(-2.4708469169133954315e+00_dp),        &
                                                          YMM4r8_t(4.7914889422856814203e-01_dp),         &
                                                          YMM4r8_t(-3.7384991926068969150e-03_dp),        &
                                                          YMM4r8_t(-2.6801520353328635310e-03_dp),        &
                                                          YMM4r8_t(9.9168777670983678974e-05_dp),         &
                                                          YMM4r8_t(-2.1877128189032726730e-06_dp)]
      type(YMM4r8_t), dimension(0:6), save ::  calci0_qq=[YMM4r8_t(-3.1446690275135491500e+01_dp),        & 
                                                          YMM4r8_t(8.5539563258012929600e+01_dp),         &
                                                          YMM4r8_t(-6.0228002066743340583e+01_dp),        &
                                                          YMM4r8_t(1.3982595353892851542e+01_dp),         &
                                                          YMM4r8_t(-1.1151759188741312645e+00_dp),        &
                                                          YMM4r8_t(3.2547697594819615062e-02_dp),         &
                                                          YMM4r8_t(-5.5194330231005480228e-04_dp)] 
    
    !!
    !! calci1_ymm4r8 constant arrays (saved)
    !!
     !dir$ attributes align : 64 :: calci1_p
     !dir$ attributes align : 64 :: calci1_q
     !dir$ attributes align : 64 :: calci1_pp
     !dir$ attributes align : 64 :: calci1_qq
     type(YMM4r8_t), dimension(0:14), save :: calci1_p = [YMM4r8_t(-1.9705291802535139930e-19_dp),           &
                                                          YMM4r8_t(-6.5245515583151902910e-16_dp),           &
                                                          YMM4r8_t(-1.1928788903603238754e-12_dp),           &
                                                          YMM4r8_t(-1.4831904935994647675e-09_dp),           &
                                                          YMM4r8_t(-1.3466829827635152875e-06_dp),           &
                                                          YMM4r8_t(-9.1746443287817501309e-04_dp),           &
                                                          YMM4r8_t(-4.7207090827310162436e-01_dp),           &
                                                          YMM4r8_t(-1.8225946631657315931e+02_dp),           &
                                                          YMM4r8_t(-5.1894091982308017540e+04_dp),           &
                                                          YMM4r8_t(-1.0588550724769347106e+07_dp),           &
                                                          YMM4r8_t(-1.4828267606612366099e+09_dp),           &
                                                          YMM4r8_t(-1.3357437682275493024e+11_dp),           &
                                                          YMM4r8_t(-6.9876779648010090070e+12_dp),           &
                                                          YMM4r8_t(-1.7732037840791591320e+14_dp),           &
                                                          YMM4r8_t(-1.4577180278143463643e+15_dp)]
     type(YMM4r8_t), dimension(0:4), save ::  calci1_q = [YMM4r8_t(-4.0076864679904189921e+03_dp),           &
                                                          YMM4r8_t(7.4810580356655069138e+06_dp),            &
                                                          YMM4r8_t(-8.0059518998619764991e+09_dp),           &
                                                          YMM4r8_t(4.8544714258273622913e+12_dp),            &
                                                          YMM4r8_t(-1.3218168307321442305e+15_dp)]
     type(YMM4r8_t), dimension(0:7), save ::  calci1_pp =[YMM4r8_t(-6.0437159056137600000e-02_dp),           &
                                                          YMM4r8_t(4.5748122901933459000e-01_dp),            &
                                                          YMM4r8_t(-4.2843766903304806403e-01_dp),           &
                                                          YMM4r8_t(9.7356000150886612134e-02_dp),            &
                                                          YMM4r8_t(-3.2457723974465568321e-03_dp),           &
                                                          YMM4r8_t(-3.6395264712121795296e-04_dp),           &
                                                          YMM4r8_t( 1.6258661867440836395e-05_dp),           &
                                                          YMM4r8_t(-3.6347578404608223492e-07_dp)]
     type(YMM4r8_t), dimension(0:5), save ::  calci1_qq =[YMM4r8_t(-3.8806586721556593450e+00_dp),           &
                                                          YMM4r8_t(3.2593714889036996297e+00_dp),            &
                                                          YMM4r8_t(-8.5017476463217924408e-01_dp),           &
                                                          YMM4r8_t(7.4212010813186530069e-02_dp),            &
                                                          YMM4r8_t(-2.2835624489492512649e-03_dp),           &
                                                          YMM4r8_t(3.7510433111922824643e-05_dp)]
    !!
    !! calck0_ymm4r8 constant arrays (saved)
    !!
    !dir$ attributes align : 64 :: calck0_p
    !dir$ attributes align : 64 :: calck0_q
    !dir$ attributes align : 64 :: calck0_f
    !dir$ attributes align : 64 :: calck0_g
    !dir$ attributes align : 64 :: calck0_pp
    !dir$ attributes align : 64 :: calck0_qq
     type(YMM4r8_t), dimension(0:5), save ::  calck0_p = [YMM4r8_t(5.8599221412826100000e-04_dp),            & 
	                                                  YMM4r8_t(1.3166052564989571850e-01_dp),            &
                                                          YMM4r8_t(1.1999463724910714109e+01_dp),              & 
                                                          YMM4r8_t(4.6850901201934832188e+02_dp),              & 
                                                          YMM4r8_t(5.9169059852270512312e+03_dp),              &
                                                          YMM4r8_t(2.4708152720399552679e+03_dp)]
     type(YMM4r8_t), dimension(0:1), save ::  calck0_q = [YMM4r8_t(-2.4994418972832303646e+02_dp),           & 
	                                                  YMM4r8_t(2.1312714303849120380e+04_dp)] 
     type(YMM4r8_t), dimension(0:3), save ::  calck0_f = [YMM4r8_t(-1.6414452837299064100e+00_dp),           &
	                                                  YMM4r8_t(-2.9601657892958843866e+02_dp),           &
                                                          YMM4r8_t(-1.7733784684952985886e+04_dp),           &           
                                                          YMM4r8_t(-4.0320340761145482298e+05_dp)]
     type(YMM4r8_t), dimension(0:2), save ::  calck0_g = [YMM4r8_t(2.5064972445877992730e+02_dp),            &
	                                                 [YMM4r8_t(2.9865713163054025489e+04_dp),            & 
                                                         [YMM4r8_t(-1.6128136304458193998e+06_dp)]
     type(YMM4r8_t), dimension(0:9), save ::  calck0_pp= [YMM4r8_t(1.1394980557384778174e+02_dp),            & 
	                                                  YMM4r8_t(3.6832589957340267940e+03_dp),            & 
                                                          YMM4r8_t(3.1075408980684392399e+04_dp),            & 
                                                          YMM4r8_t(1.0577068948034021957e+05_dp),            & 
                                                          YMM4r8_t(1.7398867902565686251e+05_dp),            &
                                                          YMM4r8_t(1.5097646353289914539e+05_dp),            & 
                                                          YMM4r8_t(7.1557062783764037541e+04_dp),            & 
                                                          YMM4r8_t(1.8321525870183537725e+04_dp),            & 
                                                          YMM4r8_t(2.3444738764199315021e+03_dp),            & 
                                                          YMM4r8_t(1.1600249425076035558e+02_dp)]
     type(YMM4r8_t), dimension(0:9), save ::  calck0_qq= [YMM4r8_t(2.0013443064949242491e+02_dp),            & 
	                                                  YMM4r8_t(4.4329628889746408858e+03_dp),            & 
                                                          YMM4r8_t(3.1474655750295278825e+04_dp),            &
                                                          YMM4r8_t(9.7418829762268075784e+04_dp),            & 
                                                          YMM4r8_t(1.5144644673520157801e+05_dp),            & 
                                                          YMM4r8_t(1.2689839587977598727e+05_dp),            & 
                                                          YMM4r8_t(5.8824616785857027752e+04_dp),            & 
                                                          YMM4r8_t(1.4847228371802360957e+04_dp),            & 
                                                          YMM4r8_t(1.8821890840982713696e+03_dp),            & 
                                                          YMM4r8_t(9.2556599177304839811e+01_dp)]
    !!
    !! calck0_ymm4r8 constant arrays (saved)
    !!
    !dir$ attributes align : 64 :: calck1_p
    !dir$ attributes align : 64 :: calck1_q
    !dir$ attributes align : 64 :: calck1_f
    !dir$ attributes align : 64 :: calck1_g
    !dir$ attributes align : 64 :: calck1_pp
    !dir$ attributes align : 64 :: calck1_qq
     type(YMM4r8_t), dimension(0:4), save ::  calck1_p = [YMM4r8_t(4.8127070456878442310e-1_dp), & 
                                                          YMM4r8_t(9.9991373567429309922e+1_dp), &
                                                          YMM4r8_t(7.1885382604084798576e+3_dp), &
                                                          YMM4r8_t(1.7733324035147015630e+5_dp), &
                                                          YMM4r8_t(7.1938920065420586101e+5_dp)]
     type(YMM4r8_t), dimension(0:2), save ::  calck1_q = [YMM4r8_t(2.8143915754538725829e+2_dp), &
                                                          YMM4r8_t(3.7264298672067697862e+4_dp), &
                                                          YMM4r8_t(-2.2149374878243304548e+6_dp)]
     type(YMM4r8_t), dimension(0:4), save ::  calck1_f = [YMM4r8_t(-2.2795590826955002390e-1_dp), &
                                                          YMM4r8_t(-5.3103913335180275253e+1_dp), &
                                                          YMM4r8_t(-4.5051623763436087023e+3_dp), &
                                                          YMM4r8_t(-1.4758069205414222471e+5_dp), &
                                                          YMM4r8_t(-1.3531161492785421328e+6_dp)]
     type(YMM4r8_t), dimension(0:2), save ::  calck1_g = [YMM4r8_t(3.0507151578787595807e+2_dp),  &
                                                          YMM4r8_t(4.3117653211351080007e+4_dp),  &
                                                          YMM4r8_t(-2.7062322985570842656e+6_dp)]
     type(YMM4r8_t), dimension(0:10), save :: calck1_pp =[YMM4r8_t(6.4257745859173138767e-2_dp),  &
                                                          YMM4r8_t(7.5584584631176030810e+0_dp),  &
                                                          YMM4r8_t(1.3182609918569941308e+2_dp),  &
                                                          YMM4r8_t(8.1094256146537402173e+2_dp),  &
                                                          YMM4r8_t(2.3123742209168871550e+3_dp),  &
                                                          YMM4r8_t(3.4540675585544584407e+3_dp),  &
                                                          YMM4r8_t(2.8590657697910288226e+3_dp),  &
                                                          YMM4r8_t(1.3319486433183221990e+3_dp),  &
                                                          YMM4r8_t(3.4122953486801312910e+2_dp),  &
                                                          YMM4r8_t(4.4137176114230414036e+1_dp),  &
                                                          YMM4r8_t(2.2196792496874548962e+0_dp)] 
     type(YMM4r8_t), dimension(0:9), save :: calck1_qq  =[YMM4r8_t(3.6001069306861518855e+1_dp),  &
                                                          YMM4r8_t(3.3031020088765390854e+2_dp),  &
                                                          YMM4r8_t(1.2082692316002348638e+3_dp),  &
                                                          YMM4r8_t(2.1181000487171943810e+3_dp),  &
                                                          YMM4r8_t(1.9448440788918006154e+3_dp),  &
                                                          YMM4r8_t(9.6929165726802648634e+2_dp),  &
                                                          YMM4r8_t(2.5951223655579051357e+2_dp),  &
                                                          YMM4r8_t(3.4552228452758912848e+1_dp),  &
                                                          YMM4r8_t(1.7710478032601086579e+0_dp)]    
    !!
    !! caljy0_ymm4r8 constant arrays (saved)
    !!  
     !dir$ attributes align : 64 :: caljy0_plg
     !dir$ attributes align : 64 :: caljy0_qlg
     !dir$ attributes align : 64 :: caljy0_pj0
     !dir$ attributes align : 64 :: caljy0_qj0
     !dir$ attributes align : 64 :: caljy0_pj1
     !dir$ attributes align : 64 :: caljy0_qj1
     !dir$ attributes align : 64 :: caljy0_py0
     !dir$ attributes align : 64 :: caljy0_qy0
     !dir$ attributes align : 64 :: caljy0_py1
     !dir$ attributes align : 64 :: caljy0_qy1
     !dir$ attributes align : 64 :: caljy0_py2
     !dir$ attributes align : 64 :: caljy0_qy2
     !dir$ attributes align : 64 :: caljy0_p0
     !dir$ attributes align : 64 :: caljy0_q0
     !dir$ attributes align : 64 :: caljy0_p1
     !dir$ attributes align : 64 :: caljy0_q1
     type(YMM4r8_t), dimension(0:3), save ::  caljy0_plg=[YMM4r8_t(2.4562334077563243311e+1_dp), &
	                                                  YMM4r8_t(2.3642701335621505212e+2_dp), &
                                                          YMM4r8_t(-5.4989956895857911039e+2_dp),&
                                                          YMM4r8_t(3.5687548468071500413e+2_dp)]
     type(YMM4r8_t), dimension(0:3), save ::  caljy0_qlg=[YMM4r8_t(3.5553900764052419184e+1_dp), &
	                                                  YMM4r8_t(1.9400230218539473193e+2_dp), &
                                                          YMM4r8_t(-3.3442903192607538956e+2_dp),&
                                                          YMM4r8_t(1.7843774234035750207e+2_dp)]  
     type(YMM4r8_t), dimension(0:6), save ::  caljy0_pj0=[YMM4r8_t(6.6302997904833794242e+6_dp), &
	                                                  YMM4r8_t(-6.2140700423540120665e+8_dp),& 
                                                          YMM4r8_t(2.7282507878605942706e+10_dp),&
                                                          YMM4r8_t(-4.1298668500990866786e+11_dp),& 
                                                          YMM4r8_t(-1.2117036164593528341e-1_dp), &
                                                          YMM4r8_t(1.0344222815443188943e+2_dp),  &
                                                          YMM4r8_t(-3.6629814655107086448e+4_dp)]
     type(YMM4r8_t), dimension(0:4), save ::  caljy0_qj0=[YMM4r8_t(4.5612696224219938200e+5_dp),  & 
	                                                  YMM4r8_t(1.3985097372263433271e+8_dp),  & 
                                                          YMM4r8_t(2.6328198300859648632e+10_dp), & 
                                                          YMM4r8_t(2.3883787996332290397e+12_dp), & 
                                                          YMM4r8_t(9.3614022392337710626e+2_dp)]
     type(YMM4r8_t), dimension(0:7), save ::  caljy0_pj1=[YMM4r8_t(4.4176707025325087628e+3_dp),  & 
	                                                  YMM4r8_t(1.1725046279757103576e+4_dp),  &
                                                          YMM4r8_t(1.0341910641583726701e+4_dp),  &
                                                          YMM4r8_t(-7.2879702464464618998e+3_dp), & 
                                                          YMM4r8_t(-1.2254078161378989535e+4_dp), &
                                                          YMM4r8_t(-1.8319397969392084011e+3_dp), & 
                                                          YMM4r8_t(4.8591703355916499363e+1_dp),  & 
                                                          YMM4r8_t(7.4321196680624245801e+2_dp)]
     type(YMM4r8_t), dimension(0:6), save ::  caljy0_qj1=[YMM4r8_t(3.3307310774649071172e+2_dp),  &
	                                                  YMM4r8_t(-2.9458766545509337327e+3_dp), &
                                                          YMM4r8_t(1.8680990008359188352e+4_dp),  &
                                                          YMM4r8_t(-8.4055062591169562211e+4_dp), & 
                                                          YMM4r8_t(2.4599102262586308984e+5_dp),  &
                                                          YMM4r8_t(-3.5783478026152301072e+5_dp), & 
                                                          YMM4r8_t(-2.5258076240801555057e+1_dp)]
     type(YMM4r8_t), dimension(0:5), save ::  caljy0_py0=[YMM4r8_t(1.0102532948020907590e+4_dp),  &
	                                                  YMM4r8_t(-2.1287548474401797963e+6_dp), & 
                                                          YMM4r8_t(2.0422274357376619816e+8_dp),  &
                                                          YMM4r8_t(-8.3716255451260504098e+9_dp), & 
                                                          YMM4r8_t(1.0723538782003176831e+11_dp), &
                                                          YMM4r8_t(-1.8402381979244993524e+1_dp)]
     type(YMM4r8_t), dimension(0:4), save ::  caljy0_qy0=[YMM4r8_t(6.6475986689240190091e+2_dp),  & 
	                                                  YMM4r8_t(2.3889393209447253406e+5_dp),  & 
                                                          YMM4r8_t(5.5662956624278251596e+7_dp),  & 
                                                          YMM4r8_t(8.1617187777290363573e+9_dp),  & 
                                                          YMM4r8_t(5.8873865738997033405e+11_dp)]
     type(YMM4r8_t), dimension(0:6), save ::  caljy0_py1=[YMM4r8_t(1.4566865832663635920e+4_dp),  & 
	                                                  YMM4r8_t(4.6905288611678631510e+6_dp),  & 
                                                          YMM4r8_t(-6.9590439394619619534e+8_dp), &
                                                          YMM4r8_t(4.3600098638603061642e+10_dp), &
                                                          YMM4r8_t(-5.5107435206722644429e+11_dp),&
                                                          YMM4r8_t(-2.2213976967566192242e+13_dp),& 
                                                          YMM4r8_t(1.7427031242901594547e+1_dp)]
     type(YMM4r8_t), dimension(0:5), save ::  caljy0_qy1=[YMM4r8_t(8.3030857612070288823e+2_dp),  & 
	                                                  YMM4r8_t(4.0669982352539552018e+5_dp),  & 
                                                          YMM4r8_t(1.3960202770986831075e+8_dp),  & 
                                                          YMM4r8_t(3.4015103849971240096e+10_dp), &
                                                          YMM4r8_t(5.4266824419412347550e+12_dp), & 
                                                          YMM4r8_t(4.3386146580707264428e+14_dp)]
     type(YMM4r8_t), dimension(0:7), save ::  caljy0_py2=[YMM4r8_t(2.1363534169313901632e+4_dp),  &
	                                                  YMM4r8_t(-1.0085539923498211426e+7_dp), & 
                                                          YMM4r8_t(2.1958827170518100757e+9_dp),  &
                                                          YMM4r8_t(-1.9363051266772083678e+11_dp),& 
                                                          YMM4r8_t(-1.2829912364088687306e+11_dp),& 
                                                          YMM4r8_t(6.7016641869173237784e+14_dp), & 
                                                          YMM4r8_t(-8.0728726905150210443e+15_dp),&
                                                          YMM4r8_t(-1.7439661319197499338e+1_dp)]
     type(YMM4r8_t), dimension(0:6), save ::  caljy0_qy2=[YMM4r8_t(8.7903362168128450017e+2_dp),  & 
	                                                  YMM4r8_t(5.3924739209768057030e+5_dp),  &
                                                          YMM4r8_t(2.4727219475672302327e+8_dp),  & 
                                                          YMM4r8_t(8.6926121104209825246e+10_dp), &
                                                          YMM4r8_t(2.2598377924042897629e+13_dp), & 
                                                          YMM4r8_t(3.9272425569640309819e+15_dp), &
                                                          YMM4r8_t(3.4563724628846457519e+17_dp)]
     type(YMM4r8_t), dimension(0:5), save ::  caljy0_p0 =[YMM4r8_t(3.4806486443249270347e+3_dp),  & 
	                                                  YMM4r8_t(2.1170523380864944322e+4_dp),  & 
                                                          YMM4r8_t(4.1345386639580765797e+4_dp),  & 
                                                          YMM4r8_t(2.2779090197304684302e+4_dp),  &
                                                          YMM4r8_t(8.8961548424210455236e-1_dp),  & 
                                                          YMM4r8_t(1.5376201909008354296e+2_dp)]
     type(YMM4r8_t), dimension(0:4), save ::  caljy0_q0 =[YMM4r8_t(3.5028735138235608207e+3_dp),  & 
	                                                  YMM4r8_t(2.1215350561880115730e+4_dp),  & 
                                                          YMM4r8_t(4.1370412495510416640e+4_dp),  & 
                                                          YMM4r8_t(2.2779090197304684318e+4_dp),  &
                                                          YMM4r8_t(1.5711159858080893649e+2_dp)]
     type(YMM4r8_t), dimension(0:5), save ::  caljy0_p1 =[YMM4r8_t(-2.2300261666214198472e+1_dp), &
	                                                  YMM4r8_t(-1.1183429920482737611e+2_dp), & 
                                                          YMM4r8_t(-1.8591953644342993800e+2_dp), &
                                                          YMM4r8_t(-8.9226600200800094098e+1_dp), &
                                                          YMM4r8_t(-8.8033303048680751817e+3_dp), &
                                                          YMM4r8_t(-1.2441026745835638459e+00_dp)]
     type(YMM4r8_t), dimension(0:4), save ::  caljy0_q1 =[YMM4r8_t(1.4887231232283756582e+3_dp),  & 
	                                                  YMM4r8_t(7.2642780169211018836e+3_dp),  &
                                                          YMM4r8_t(1.1951131543434613647e+4_dp),  & 
                                                          YMM4r8_t(5.7105024128512061905e+3_dp),  &
                                                          YMM4r8_t(9.0593769594993125859e+1_dp)]
    !!
    !! caljy1_ymm4r8 constant arrays (saved)
    !!  
     !dir$ attributes align : 64 :: caljy1_plg
     !dir$ attributes align : 64 :: caljy1_qlg
     !dir$ attributes align : 64 :: caljy1_pj0
     !dir$ attributes align : 64 :: caljy1_qj0
     !dir$ attributes align : 64 :: caljy1_pj1
     !dir$ attributes align : 64 :: caljy1_qj1
     !dir$ attributes align : 64 :: caljy1_py0
     !dir$ attributes align : 64 :: caljy1_qy0
     !dir$ attributes align : 64 :: caljy1_py1
     !dir$ attributes align : 64 :: caljy1_qy1
     !dir$ attributes align : 64 :: caljy1_py2
     !dir$ attributes align : 64 :: caljy1_qy2
     !dir$ attributes align : 64 :: caljy1_p0
     !dir$ attributes align : 64 :: caljy1_q0
     !dir$ attributes align : 64 :: caljy1_p1
     !dir$ attributes align : 64 :: caljy1_q1
     type(YMM4r8_t), dimension(0:3), save ::  caljy1_plg =[YMM4r8_t(-2.4562334077563243311e+1_dp), &
	                                                   YMM4r8_t(2.3642701335621505212e+2_dp),  &
                                                           YMM4r8_t(-5.4989956895857911039e+2_dp), &
                                                           YMM4r8_t(3.5687548468071500413e+2_dp)]
     type(YMM4r8_t), dimension(0:3), save ::  caljy1_qlg =[YMM4r8_t(-3.5553900764052419184e+1_dp), &
	                                                   YMM4r8_t(1.9400230218539473193e+2_dp),  &
                                                           YMM4r8_t(-3.3442903192607538956e+2_dp), &
                                                           YMM4r8_t(1.7843774234035750207e+2_dp)]
     type(YMM4r8_t), dimension(0:6), save ::  caljy1_pj0 =[YMM4r8_t(9.8062904098958257677e+5_dp),  &
	                                                   YMM4r8_t(-1.1548696764841276794e+8_dp), & 
                                                           YMM4r8_t(6.6781041261492395835e+9_dp),  &
                                                           YMM4r8_t(-1.4258509801366645672e+11_dp),& 
                                                           YMM4r8_t(-4.4615792982775076130e+3_dp), & 
                                                           YMM4r8_t(1.0650724020080236441e+1_dp),  &
                                                           YMM4r8_t(-1.0767857011487300348e-2_dp)]
     type(YMM4r8_t), dimension(0:4), save ::  caljy1_qj0 =[YMM4r8_t(5.9117614494174794095e+5_dp),  & 
	                                                   YMM4r8_t(2.0228375140097033958e+8_dp),  & 
                                                           YMM4r8_t(4.2091902282580133541e+10_dp), & 
                                                           YMM4r8_t(4.1868604460820175290e+12_dp), & 
                                                           YMM4r8_t(1.0742272239517380498e+03_dp)]
     type(YMM4r8_t), dimension(0:7), save ::  caljy1_pj1 =[YMM4r8_t(4.6179191852758252280e+00_dp), &
	                                                   YMM4r8_t(-7.1329006872560947377e+3_dp), &
                                                           YMM4r8_t(4.5039658105749078904e+6_dp),  &
                                                           YMM4r8_t(-1.4437717718363239107e+9_dp), &
                                                           YMM4r8_t(2.3569285397217157313e+11_dp), &
                                                           YMM4r8_t(-1.6324168293282543629e+13_dp),&
                                                           YMM4r8_t(1.1357022719979468624e+14_dp), & 
                                                           YMM4r8_t(1.0051899717115285432e+15_dp)]
     type(YMM4r8_t), dimension(0:6), save ::  caljy1_qj1 =[YMM4r8_t(1.1267125065029138050e+6_dp),  & 
	                                                   YMM4r8_t(6.4872502899596389593e+8_dp),  &
                                                           YMM4r8_t(2.7622777286244082666e+11_dp), & 
                                                           YMM4r8_t(8.4899346165481429307e+13_dp), &
                                                           YMM4r8_t(1.7128800897135812012e+16_dp), & 
                                                           YMM4r8_t(1.7253905888447681194e+18_dp), & 
                                                           YMM4r8_t(1.3886978985861357615e+3_dp)]
     type(YMM4r8_t), dimension(0:6), save ::  caljy1_py0 =[YMM4r8_t(2.2157953222280260820e+5_dp),  &
	                                                   YMM4r8_t(-5.9157479997408395984e+7_dp), & 
                                                           YMM4r8_t(7.2144548214502560419e+9_dp),  &
                                                           YMM4r8_t(-3.7595974497819597599e+11_dp),&
                                                           YMM4r8_t(5.4708611716525426053e+12_dp), & 
                                                           YMM4r8_t(4.0535726612579544093e+13_dp), & 
                                                           YMM4r8_t(-3.1714424660046133456e+2_dp)]
     type(YMM4r8_t), dimension(0:5), save ::  caljy1_qy0 =[YMM4r8_t(8.2079908168393867438e+2_dp),  & 
	                                                   YMM4r8_t(3.8136470753052572164e+5_dp),  &
                                                           YMM4r8_t(1.2250435122182963220e+8_dp),  & 
                                                           YMM4r8_t(2.7800352738690585613e+10_dp), &
                                                           YMM4r8_t(4.1272286200406461981e+12_dp), & 
                                                           YMM4r8_t(3.0737873921079286084e+14_dp)]
     type(YMM4r8_t), dimension(0:8), save ::  caljy1_py1 =[YMM4r8_t(1.9153806858264202986e+6_dp),  &
	                                                   YMM4r8_t(-1.1957961912070617006e+9_dp), & 
                                                           YMM4r8_t(3.7453673962438488783e+11_dp), &
                                                           YMM4r8_t(-5.9530713129741981618e+13_dp),& 
                                                           YMM4r8_t(4.0686275289804744814e+15_dp), &
                                                           YMM4r8_t(-2.3638408497043134724e+16_dp),&
                                                           YMM4r8_t(-5.6808094574724204577e+18_dp),& 
                                                           YMM4r8_t(1.1514276357909013326e+19_dp), & 
                                                           YMM4r8_t(-1.2337180442012953128e+3_dp)]
     type(YMM4r8_t), dimension(0:7), save ::  caljy1_qy1 =[YMM4r8_t(1.2855164849321609336e+3_dp),  & 
	                                                   YMM4r8_t(1.0453748201934079734e+6_dp),  & 
                                                           YMM4r8_t(6.3550318087088919566e+8_dp),  & 
                                                           YMM4r8_t(3.0221766852960403645e+11_dp), & 
                                                           YMM4r8_t(1.1187010065856971027e+14_dp), & 
                                                           YMM4r8_t(3.0837179548112881950e+16_dp), &
                                                           YMM4r8_t(5.6968198822857178911e+18_dp), & 
                                                           YMM4r8_t(5.3321844313316185697e+20_dp)]
     type(YMM4r8_t), dimension(0:5), save ::  caljy1_p0  =[YMM4r8_t(1.0982405543459346727e+5_dp),  &
	                                                   YMM4r8_t(-1.5235293511811373833e+6_dp), &
                                                           YMM4r8_t(-6.6033732483649391093e+06_dp),&
                                                           YMM4r8_t(-9.9422465050776411957e+6_dp), &
                                                           YMM4r8_t(-4.4357578167941278571e+6_dp), &
                                                           YMM4r8_t(-1.6116166443246101165e+3_dp)]
     type(YMM4r8_t), dimension(0:5), save ::  caljy1_q0  =[YMM4r8_t(-1.0726385991103820119e+5_dp), &
	                                                   YMM4r8_t(-1.5118095066341608816e+6_dp), &
                                                           YMM4r8_t(-6.5853394797230870728e+6_dp), &
                                                           YMM4r8_t(-9.9341243899345856590e+6_dp), & 
                                                           YMM4r8_t(-4.4357578167941278568e+6_dp), &
                                                           YMM4r8_t(-1.4550094401904961825e+3_dp)]
     type(YMM4r8_t), dimension(0:5), save ::  caljy1_p1  =[YMM4r8_t(1.7063754290207680021e+3_dp),  & 
	                                                   YMM4r8_t(1.8494262873223866797e+4_dp),  & 
                                                           YMM4r8_t(6.6178836581270835179e+4_dp),  & 
                                                           YMM4r8_t(8.5145160675335701966e+4_dp),  &
                                                           YMM4r8_t(3.3220913409857223519e+4_dp),  & 
                                                           YMM4r8_t(3.5265133846636032186e+1_dp)]
     type(YMM4r8_t), dimension(0:5), save ::  caljy1_q1  =[YMM4r8_t(3.7890229745772202641e+4_dp),  & 
	                                                   YMM4r8_t(4.0029443582266975117e+5_dp),  &
                                                           YMM4r8_t(1.4194606696037208929e+6_dp),  & 
                                                           YMM4r8_t(1.8194580422439972989e+6_dp),  &
                                                           YMM4r8_t(7.0871281941028743574e+5_dp),  & 
                                                           YMM4r8_t(8.6383677696049909675e+2_dp)]
     
    !!
    !! calck0_ymm4r8 constant arrays (saved)
    !!  
     !dir$ attributes align : 64 :: calck0_p
     !dir$ attributes align : 64 :: calck0_q
     !dir$ attributes align : 64 :: calck0_f
     !dir$ attributes align : 64 :: calck0_g
     !dir$ attributes align : 64 :: calck0_pp
     !dir$ attributes align : 64 :: calck0_qq
     type(YMM4r8_t), dimension(0:5), save :: calck0_p    =[YMM4r8_t(5.8599221412826100000e-04_dp), & 
	                                                   YMM4r8_t(1.3166052564989571850e-01_dp), & 
                                                           YMM4r8_t(1.1999463724910714109e+01_dp), & 
                                                           YMM4r8_t(4.6850901201934832188e+02_dp), & 
                                                           YMM4r8_t(5.9169059852270512312e+03_dp), & 
                                                           YMM4r8_t(2.4708152720399552679e+03_dp)]
     type(YMM4r8_t), dimension(0:1), save :: calck0_q    =[YMM4r8_t(-2.4994418972832303646e+02_dp),& 
	                                                   YMM4r8_t(2.1312714303849120380e+04_dp)]
     type(YMM4r8_t), dimension(0:3), save :: calck0_f    =[YMM4r8_t(-1.6414452837299064100e+00_dp),&
	                                                   YMM4r8_t(-2.9601657892958843866e+02_dp),& 
                                                           YMM4r8_t(-1.7733784684952985886e+04_dp),&
                                                           YMM4r8_t(-4.0320340761145482298e+05_dp)]
     type(YMM4r8_t), dimension(0:2), save :: calck0_g    =[YMM4r8_t(-2.5064972445877992730e+02_dp),& 
	                                                   YMM4r8_t(2.9865713163054025489e+04_dp), &
                                                           YMM4r8_t(-1.6128136304458193998e+06_dp)]
     type(YMM4r8_t), dimension(0:9), save :: calck0_pp   =[YMM4r8_t(1.1394980557384778174e+02_dp), & 
	                                                   YMM4r8_t(3.6832589957340267940e+03_dp), & 
                                                           YMM4r8_t(3.1075408980684392399e+04_dp), & 
                                                           YMM4r8_t(1.0577068948034021957e+05_dp), & 
                                                           YMM4r8_t(1.7398867902565686251e+05_dp), &
                                                           YMM4r8_t(1.5097646353289914539e+05_dp), & 
                                                           YMM4r8_t(7.1557062783764037541e+04_dp), & 
                                                           YMM4r8_t(1.8321525870183537725e+04_dp), & 
                                                           YMM4r8_t(2.3444738764199315021e+03_dp), & 
                                                           YMM4r8_t(1.1600249425076035558e+02_dp)]
     type(YMM4r8_t), dimension(0:9), save :: calck0_qq   =[YMM4r8_t(2.0013443064949242491e+02_dp), & 
	                                                   YMM4r8_t(4.4329628889746408858e+03_dp), & 
                                                           YMM4r8_t(3.1474655750295278825e+04_dp), & 
                                                           YMM4r8_t(9.7418829762268075784e+04_dp), & 
                                                           YMM4r8_t(1.5144644673520157801e+05_dp), & 
                                                           YMM4r8_t(1.2689839587977598727e+05_dp), & 
                                                           YMM4r8_t(5.8824616785857027752e+04_dp), & 
                                                           YMM4r8_t(1.4847228371802360957e+04_dp), &
                                                           YMM4r8_t(1.8821890840982713696e+03_dp), & 
                                                           YMM4r8_t(9.2556599177304839811e+01_dp)] 
    !!
    !! calck1_ymm4r8 constant arrays (saved)
    !!   
     !dir$ attributes align : 64 :: calck0_p
     !dir$ attributes align : 64 :: calck0_q
     !dir$ attributes align : 64 :: calck0_f
     !dir$ attributes align : 64 :: calck0_g
     !dir$ attributes align : 64 :: calck0_pp
     !dir$ attributes align : 64 :: calck0_qq
     type(YMM4r8_t), dimension(0:4), save :: calck1_p  =[YMM4r8_t(4.8127070456878442310e-1_dp),   & 
	                                                 YMM4r8_t(9.9991373567429309922e+1_dp),   & 
                                                         YMM4r8_t(7.1885382604084798576e+3_dp),   & 
                                                         YMM4r8_t(1.7733324035147015630e+5_dp),   & 
                                                         YMM4r8_t(7.1938920065420586101e+5_dp)]
     type(YMM4r8_t), dimension(0:2), save :: calck1_q  =[YMM4r8_t(-2.8143915754538725829e+2_dp),  & 
	                                                 YMM4r8_t(3.7264298672067697862e+4_dp),   & 
                                                         YMM4r8_t(-2.2149374878243304548e+6_dp)]
     type(YMM4r8_t), dimension(0:4), save :: calck1_f  =[YMM4r8_t(-2.2795590826955002390e-1_dp),  &
	                                                 YMM4r8_t(-5.3103913335180275253e+1_dp),  & 
                                                         YMM4r8_t(-4.5051623763436087023e+3_dp),  &
                                                         YMM4r8_t(-1.4758069205414222471e+5_dp),  &
                                                         YMM4r8_t(-1.3531161492785421328e+6_dp)]
     type(YMM4r8_t), dimension(0:2), save :: calck1_g  =[YMM4r8_t(-3.0507151578787595807e+2_dp),  & 
	                                                 YMM4r8_t(4.3117653211351080007e+4_dp),   & 
                                                         YMM4r8_t(-2.7062322985570842656e+6_dp)]
     type(YMM4r8_t), dimension(0:10), save :: calck1_pp=[YMM4r8_t(6.4257745859173138767e-2_dp),   &
	                                                 YMM4r8_t(7.5584584631176030810e+0_dp),   & 
                                                         YMM4r8_t(1.3182609918569941308e+2_dp),   & 
                                                         YMM4r8_t(8.1094256146537402173e+2_dp),   &
                                                         YMM4r8_t(2.3123742209168871550e+3_dp),   & 
                                                         YMM4r8_t(3.4540675585544584407e+3_dp),   & 
                                                         YMM4r8_t(2.8590657697910288226e+3_dp),   & 
                                                         YMM4r8_t(1.3319486433183221990e+3_dp),   & 
                                                         YMM4r8_t(3.4122953486801312910e+2_dp),   & 
                                                         YMM4r8_t(4.4137176114230414036e+1_dp),   & 
                                                         YMM4r8_t(2.2196792496874548962e+0_dp)]
     type(YMM4r8_t), dimension(0:8), save :: calck1_qq =[YMM4r8_t(3.6001069306861518855e+1_dp),   & 
	                                                 YMM4r8_t(3.3031020088765390854e+2_dp),   & 
                                                         YMM4r8_t(1.2082692316002348638e+3_dp),   & 
                                                         YMM4r8_t(2.1181000487171943810e+3_dp),   &
                                                         YMM4r8_t(1.9448440788918006154e+3_dp),   & 
                                                         YMM4r8_t(9.6929165726802648634e+2_dp),   & 
                                                         YMM4r8_t(2.5951223655579051357e+2_dp),   & 
                                                         YMM4r8_t(3.4552228452758912848e+1_dp),   & 
                                                         YMM4r8_t(1.7710478032601086579e+0_dp)]
      
     contains
     
!! =============================================================================================================== //
!!                                  'Saved' arrays preload_calcei routines.
!!================================================================================================================ //  


      pure function preload_calcei_a() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_a
              !dir$ attributes forceinline :: preload_calcei_a
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_a
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v    = calcei_a(0).v+calcei_a(1).v
              t1.v    = calcei_a(2).v+calcei_a(3).v
              t2.v    = calcei_a(4).v+calcei_a(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_calcei_a 
       
       
       pure function preload_calcei_b() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_b
              !dir$ attributes forceinline :: preload_calcei_b
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_b
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v    = calcei_b(0).v+calcei_b(1).v
              t1.v    = calcei_b(2).v+calcei_b(3).v
              t2.v    = calcei_b(4).v+calcei_b(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_calcei_b 
       
       
      pure function preload_calcei_c() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_c
              !dir$ attributes forceinline :: preload_calcei_c
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_c
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v    = calcei_c(0).v+calcei_c(1).v
              t1.v    = calcei_c(2).v+calcei_c(3).v
              t2.v    = calcei_c(4).v+calcei_c(5).v
              t3.v    = calcei_c(6).v+calcei_c(7).v+ &
                        calcei_c(8).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_calcei_c 
       
       
       pure function preload_calcei_d() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_d
              !dir$ attributes forceinline :: preload_calcei_d
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_d
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v    = calcei_d(0).v+calcei_d(1).v
              t1.v    = calcei_d(2).v+calcei_d(3).v
              t2.v    = calcei_d(4).v+calcei_d(5).v
              t3.v    = calcei_d(6).v+calcei_d(7).v+ &
                        calcei_d(8).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_calcei_d 
       
       
       pure function preload_calcei_e() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_e
              !dir$ attributes forceinline :: preload_calcei_e
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_e
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3,t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3,t4
              t0.v    = calcei_e(0).v+calcei_e(1).v
              t1.v    = calcei_e(2).v+calcei_e(3).v
              t2.v    = calcei_e(4).v+calcei_e(5).v
              t3.v    = calcei_e(6).v+calcei_e(7).v
              t4.v    = calcei_e(8).v+calcei_e(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+t4.v
       end function preload_calcei_e 
       
       
       pure function preload_calcei_f() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_f
              !dir$ attributes forceinline :: preload_calcei_f
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_f
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3,t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3,t4
              t0.v    = calcei_f(0).v+calcei_f(1).v
              t1.v    = calcei_f(2).v+calcei_f(3).v
              t2.v    = calcei_f(4).v+calcei_f(5).v
              t3.v    = calcei_f(6).v+calcei_f(7).v
              t4.v    = calcei_f(8).v+calcei_f(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+t4.v
       end function preload_calcei_f 
       
       
       pure function preload_calcei_plg() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_plg
              !dir$ attributes forceinline :: preload_calcei_plg
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_plg
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v    = calcei_plg(0).v+calcei_plg(1).v
              t1.v    = calcei_plg(2).v+calcei_plg(3).v
              summa.v = t0.v+t1.v
       end function preload_calcei_plg 
       
       
       pure function preload_calcei_qlg() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_qlg
              !dir$ attributes forceinline :: preload_calcei_qlg
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_qlg
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v    = calcei_qlg(0).v+calcei_qlg(1).v
              t1.v    = calcei_qlg(2).v+calcei_qlg(3).v
              summa.v = t0.v+t1.v
       end function preload_calcei_qlg 
       
       
       pure function preload_calcei_p() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_p
              !dir$ attributes forceinline :: preload_calcei_p
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_p
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3,t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3,t4
              t0.v    = calcei_p(0).v+calcei_p(1).v
              t1.v    = calcei_p(2).v+calcei_p(3).v
              t2.v    = calcei_p(4).v+calcei_p(5).v
              t3.v    = calcei_p(6).v+calcei_p(7).v
              t4.v    = calcei_p(8).v+calcei_p(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+t4.v
       end function preload_calcei_p 
       
       
       pure function preload_calcei_q() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_q
              !dir$ attributes forceinline :: preload_calcei_q
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_q
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3,t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3,t4
              t0.v    = calcei_q(0).v+calcei_q(1).v
              t1.v    = calcei_q(2).v+calcei_q(3).v
              t2.v    = calcei_q(4).v+calcei_q(5).v
              t3.v    = calcei_q(6).v+calcei_q(7).v
              t4.v    = calcei_q(8).v+calcei_q(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+t4.v
       end function preload_calcei_q 
       
       
       pure function preload_calcei_r() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_r
              !dir$ attributes forceinline :: preload_calcei_r
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_r
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3,t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3,t4
              t0.v    = calcei_r(0).v+calcei_r(1).v
              t1.v    = calcei_r(2).v+calcei_r(3).v
              t2.v    = calcei_r(4).v+calcei_r(5).v
              t3.v    = calcei_r(6).v+calcei_r(7).v
              t4.v    = calcei_r(8).v+calcei_r(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+t4.v
       end function preload_calcei_r 
       
       
       pure function preload_calcei_s() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_s
              !dir$ attributes forceinline :: preload_calcei_s
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_s
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v    = calcei_r(0).v+calcei_r(1).v
              t1.v    = calcei_r(2).v+calcei_r(3).v
              t2.v    = calcei_r(4).v+calcei_r(5).v
              t3.v    = calcei_r(6).v+calcei_r(7).v+ &
                        calcei_r(8).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_calcei_s 
       
       
       pure function preload_calcei_p1() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_p1
              !dir$ attributes forceinline :: preload_calcei_p1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_p1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3,t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3,t4
              t0.v    = calcei_p1(0).v+calcei_p1(1).v
              t1.v    = calcei_p1(2).v+calcei_p1(3).v
              t2.v    = calcei_p1(4).v+calcei_p1(5).v
              t3.v    = calcei_p1(6).v+calcei_p1(7).v
              t4.v    = calcei_p1(8).v+calcei_p1(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+t4.v
       end function preload_calcei_p1 
       
       
       pure function preload_calcei_q1() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_q1
              !dir$ attributes forceinline :: preload_calcei_q1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_q1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v    = calcei_q1(0).v+calcei_q1(1).v
              t1.v    = calcei_q1(2).v+calcei_q1(3).v
              t2.v    = calcei_q1(4).v+calcei_q1(5).v
              t3.v    = calcei_q1(6).v+calcei_q1(7).v+ &
                        calcei_q1(8).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_calcei_q1 
       
       
        pure function preload_calcei_p2() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_p2
              !dir$ attributes forceinline :: preload_calcei_p2
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_p2
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3,t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3,t4
              t0.v    = calcei_p2(0).v+calcei_p2(1).v
              t1.v    = calcei_p2(2).v+calcei_p2(3).v
              t2.v    = calcei_p2(4).v+calcei_p2(5).v
              t3.v    = calcei_p2(6).v+calcei_p2(7).v
              t4.v    = calcei_p2(8).v+calcei_p2(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+t4.v
       end function preload_calcei_p2 
       
       
       pure function preload_calcei_q2() result(summa)
            
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calcei_q2
              !dir$ attributes forceinline :: preload_calcei_q2
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calcei_q2
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v    = calcei_q2(0).v+calcei_q2(1).v
              t1.v    = calcei_q2(2).v+calcei_q2(3).v
              t2.v    = calcei_q2(4).v+calcei_q2(5).v
              t3.v    = calcei_q2(6).v+calcei_q2(7).v+ &
                        calcei_q2(8).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_calcei_q2 
       
!! =============================================================================================================== //
!!                                  'Saved' arrays preload_calci0 routines.
!!================================================================================================================ //        
       
       
       pure function preload_calci0_p() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calci0_p
              !dir$ attributes forceinline :: preload_calci0_p
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calci0_p
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              !dir$ attributes align : 64 :: t4,t5,t6
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t4,t5,t6
              t0.v = calci0_p(0).v+calci0_p(1).v
              t1.v = calci0_p(2).v+calci0_p(3).v
              t2.v = calci0_p(4).v+calci0_p(5).v
              t3.v = calci0_p(6).v+calci0_p(7).v
              t4.v = calci0_p(8).v+calci0_p(9).v
              t5.v = calci0_p(10).v+calci0_p(11).v
              t6.v = calci0_p(12).v+calci0_p(13).v+ &
                     calc0_p(14).v
              summa.v = t0.v+t1.v+t2.v+t3.v+ &
                        t4.v+t5.v+t6.v
       end function preload_calci0_p
       
       
       pure function preload_calci0_q() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calci0_q
              !dir$ attributes forceinline :: preload_calci0_q
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calci0_q
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = calci0_p(0).v+calci0_p(1).v
              t1.v = calci0_p(2).v+calci0_p(3).v+ &
                     calci0_p(4)
              summa.v = t0.v+t1.v
       end function preload_calci0_q
       
       
       pure function preload_calci0_pp() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calci0_pp
              !dir$ attributes forceinline :: preload_calci0_pp
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calci0_pp
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v = calci0_pp(0).v+calci0_pp(1).v
              t1.v = calci0_pp(2).v+calci0_pp(3).v
              t2.v = calci0_pp(4).v+calci0_pp(5).v
              t3.v = calci0_pp(6).v+calci0_pp(7).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_calci0_pp
       
       
       pure function preload_calci0_qq() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calci0_qq
              !dir$ attributes forceinline :: preload_calci0_qq
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calci0_qq
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = calci0_qq(0).v+calci0_qq(1).v
              t1.v = calci0_qq(2).v+calci0_qq(3).v
              t2.v = calci0_qq(4).v+calci0_qq(5).v+ &
                     calci0_qq(6).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_calci0_qq
       
!! =============================================================================================================== //
!!                                  'Saved' arrays preload_calci1 routines.
!!================================================================================================================ //        
       
       
      pure function preload_calci1_p() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calc1_p
              !dir$ attributes forceinline :: preload_calci1_p
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calci1_p
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              !dir$ attributes align : 64 :: t4,t5,t6
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t4,t5,t6
              t0.v = calci1_p(0).v+calci1_p(1).v
              t1.v = calci1_p(2).v+calci1_p(3).v
              t2.v = calci1_p(4).v+calci1_p(5).v
              t3.v = calci1_p(6).v+calci1_p(7).v
              t4.v = calci1_p(8).v+calci1_p(9).v
              t5.v = calci1_p(10).v+calci1_p(11).v
              t6.v = calci1_p(12).v+calci1_p(13).v+ &
                     calc1_p(14).v
              summa.v = t0.v+t1.v+t2.v+t3.v+ &
                        t4.v+t5.v+t6.v
       end function preload_calci1_p       
       
       
       pure function preload_calci1_q() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calci1_q
              !dir$ attributes forceinline :: preload_calci1_q
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calci1_q
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = calci1_p(0).v+calci1_p(1).v
              t1.v = calci1_p(2).v+calci1_p(3).v+ &
                     calci1_p(4).v
              summa.v = t0.v+t1.v
       end function preload_calci1_q
       
       
       pure function preload_calci1_pp() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calci1_pp
              !dir$ attributes forceinline :: preload_calci1_pp
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calci1_pp
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v = calci1_pp(0).v+calci1_pp(1).v
              t1.v = calci1_pp(2).v+calci1_pp(3).v
              t2.v = calci1_pp(4).v+calci1_pp(5).v
              t3.v = calci1_pp(6).v+calci1_pp(7).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_calci1_pp
       
       
       pure function preload_calci1_qq() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calci1_qq
              !dir$ attributes forceinline :: preload_calci1_qq
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calci0_qq
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = calci1_qq(0).v+calci1_qq(1).v
              t1.v = calci1_qq(2).v+calci1_qq(3).v
              t2.v = calci1_qq(4).v+calci1_qq(5).v
              
              summa.v = t0.v+t1.v+t2.v
       end function preload_calci1_qq

!! =============================================================================================================== //
!!                                  'Saved' arrays preload_calck0 routines.
!!================================================================================================================ //         
     
        
      pure function preload_calck0_p() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_p
              !dir$ attributes forceinline :: preload_calck0_p
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_p
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = calck0_p(0).v+calck0_p(1).v
              t1.v = calck0_p(2).v+calck0_p(3).v
              t2.v = calck0_p(4).v+calck0_p(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_calck0_p 
       
       
       pure function preload_calck0_q() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_q
              !dir$ attributes forceinline :: preload_calck0_q
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_q
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0
              type(YMM4r8_t), automatic :: t0
              t0.v = calck0_q(0).v+calck0_q(1).v
              summa.v = t0.v
       end function preload_calck0_q 
       
       
       pure function preload_calck0_f() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_f
              !dir$ attributes forceinline :: preload_calck0_f
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_f
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = calck0_f(0).v+calck0_f(1).v
              t1.v = calck0_f(2).v+calck0_f(3).v
              summa.v = t0.v+t1.v
       end function preload_calck0_f 
       
       
       pure function preload_calck0_g() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_g
              !dir$ attributes forceinline :: preload_calck0_g
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_g
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0
              type(YMM4r8_t), automatic :: t0
              t0.v = calck0_g(0).v+calck0_g(1).v+ &
                     calck0_g(2).v
              summa.v = t0.v
       end function preload_calck0_g
       
       
         
      pure function preload_calck0_pp() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_pp
              !dir$ attributes forceinline :: preload_calck0_pp
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_pp
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              !dir$ attributes align : 64 :: t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t4
              t0.v = calck0_pp(0).v+calck0_pp(1).v
              t1.v = calck0_pp(2).v+calck0_pp(3).v
              t2.v = calck0_pp(4).v+calck0_pp(5).v
              t3.v = calck0_pp(6).v+calck0_pp(7).v
              t4.v = calck0_pp(8).v+calck0_pp(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+ &
                        t4.v
       end function preload_calck0_pp       
       
     
       pure function preload_calck0_qq() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_qq
              !dir$ attributes forceinline :: preload_calck0_qq
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_qq
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              !dir$ attributes align : 64 :: t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t4
              t0.v = calck0_qq(0).v+calck0_qq(1).v
              t1.v = calck0_qq(2).v+calck0_qq(3).v
              t2.v = calck0_qq(4).v+calck0_qq(5).v
              t3.v = calck0_qq(6).v+calck0_qq(7).v
              t4.v = calck0_qq(8).v+calck0_qq(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+ &
                        t4.v
       end function preload_calck0_qq   
       
!! =============================================================================================================== //
!!                                  'Saved' arrays preload_calck1 routines.
!!================================================================================================================ //  


       pure function preload_calck1_p() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_p
              !dir$ attributes forceinline :: preload_calck1_p
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_p
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = calck1_p(0).v+calck1_p(1).v
              t1.v = calck1_p(2).v+calck1_p(3).v+ &
                     calck1_p(4).v
              summa.v = t0.v+t1.v
       end function preload_calck1_p  
       
       
       pure function preload_calck1_q() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_q
              !dir$ attributes forceinline :: preload_calck1_q
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_q
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0
              type(YMM4r8_t), automatic :: t0
              t0.v = calck1_q(0).v+calck1_q(1).v+ &
                     calck1_q(2).v
              summa.v = t0.v
       end function preload_calck1_q   
       
       
       pure function preload_calck1_f() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_f
              !dir$ attributes forceinline :: preload_calck1_f
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_f
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = calck1_f(0).v+calck1_f(1).v
              t1.v = calck1_f(2).v+calck1_f(3).v+ &
              t2.v = calck1_f(4).v
              summa.v = t0.v+t1.v
       end function preload_calck1_f  
       
       
       pure function preload_calck1_g() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_g
              !dir$ attributes forceinline :: preload_calck1_g
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_g
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0
              type(YMM4r8_t), automatic :: t0
              t0.v = calck1_g(0).v+calck1_g(1).v+ &
                     calck1_g(2).v
              summa.v = t0.v
       end function preload_calck1_g  
       
       
       pure function preload_calck1_pp() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_pp
              !dir$ attributes forceinline :: preload_calck1_pp
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_pp
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              !dir$ attributes align : 64 :: t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t4
              t0.v = calck1_pp(0).v+calck1_pp(1).v
              t1.v = calck1_pp(2).v+calck1_pp(3).v
              t2.v = calck1_pp(4).v+calck1_pp(5).v
              t3.v = calck1_pp(6).v+calck1_pp(7).v
              t4.v = calck1_pp(8).v+calck2_pp(9).v+
                     calck1_pp(10).v
              summa.v = t0.v+t1.v+t2.v+t3.v+ &
                        t4.v
       end function preload_calck1_pp  
       
       
       pure function preload_calck1_qq() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_qq
              !dir$ attributes forceinline :: preload_calck1_qq
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_qq
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v = calck1_qq(0).v+calck1_qq(1).v
              t1.v = calck1_qq(2).v+calck1_qq(3).v
              t2.v = calck1_qq(4).v+calck1_qq(5).v
              t3.v = calck1_qq(6).v+calck1_qq(7).v+ &
                     calck1_qq(8).v
              summa.v = t0.v+t1.v+t2.v+t3.v
         
         
       end function preload_calck1_qq   
       
!! =============================================================================================================== //
!!                                  'Saved' arrays preload_caljy0 routines.
!!================================================================================================================ //   


      pure function preload_caljy0_plg() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_plg
              !dir$ attributes forceinline :: preload_calcjy0_plg
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_plg
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = caljy0_plg(0).v+caljy0_plg(1).v
              t1.v = caljy0_plg(2).v+caljy0_plg(3).v
              summa.v = t0.v+t1.v
       end function preload_caljy0_plg      
       
       
       pure function preload_caljy0_qlg() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_qlg
              !dir$ attributes forceinline :: preload_calcjy0_qlg
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_qlg
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = caljy0_qlg(0).v+caljy0_qlg(1).v
              t1.v = caljy0_qlg(2).v+caljy0_qlg(3).v
              summa.v = t0.v+t1.v
       end function preload_caljy0_qlg  
       
       
       pure function preload_caljy0_pj0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_pj0
              !dir$ attributes forceinline :: preload_calcjy0_pj0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_pj0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v = caljy0_pj0(0).v+caljy0_pj0(1).v
              t1.v = caljy0_pj0(2).v+caljy0_pj0(3).v
              t2.v = caljy0_pj0(4).v+caljy0_pj0(5).v
              t3.v = caljy0_pj0(6).v+caljy0_pj0(7).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_caljy0_pj0   
       
       
       pure function preload_caljy0_qj0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_qj0
              !dir$ attributes forceinline :: preload_calcjy0_qj0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_qj0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy0_qj0(0).v+caljy0_qj0(1).v
              t1.v = caljy0_qj0(2).v+caljy0_qj0(3).v
              t2.v = caljy0_qj0(4).v+caljy0_qj0(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy0_qj0 
       
       
       pure function preload_caljy0_pj1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_pj1
              !dir$ attributes forceinline :: preload_calcjy0_pj1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_pj1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v = caljy0_pj1(0).v+caljy0_pj1(1).v
              t1.v = caljy0_pj1(2).v+caljy0_pj1(3).v
              t2.v = caljy0_pj1(4).v+caljy0_pj1(5).v
              t3.v = caljy0_pj1(6).v+caljy0_pj1(7).v+ &
                     caljy0_pj1(8).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_caljy0_pj1   
       
       
       pure function preload_caljy0_qj1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_qj1
              !dir$ attributes forceinline :: preload_calcjy0_qj1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_qj1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy0_qj1(0).v+caljy0_qj1(1).v
              t1.v = caljy0_qj1(2).v+caljy0_qj1(3).v
              t2.v = caljy0_qj1(4).v+caljy0_qj1(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy0_qj1 
       
       
       pure function preload_caljy0_py0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_py0
              !dir$ attributes forceinline :: preload_calcjy0_py0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_py0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy0_py0(0).v+caljy0_py0(1).v
              t1.v = caljy0_py0(2).v+caljy0_py0(3).v
              t2.v = caljy0_py0(4).v+caljy0_py0(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy0_py0 
       
       
       pure function preload_caljy0_qy0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_qy0
              !dir$ attributes forceinline :: preload_calcjy0_qy0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_qy0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy0_qy0(0).v+caljy0_qy0(1).v
              t1.v = caljy0_qy0(2).v+caljy0_qy0(3).v+ &
                     caljy0_qy0(4).v
              summa.v = t0.v+t1.v
       end function preload_caljy0_qy0 
       
       
       pure function preload_caljy0_py1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_py1
              !dir$ attributes forceinline :: preload_calcjy0_py1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_py1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy0_py1(0).v+caljy0_py1(1).v
              t1.v = caljy0_py1(2).v+caljy0_py1(3).v
              t2.v = caljy0_py1(4).v+caljy0_py1(5).v+ &
                     caljy0_py1(6).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy0_py1 
       
       
       pure function preload_caljy0_qy1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_qy1
              !dir$ attributes forceinline :: preload_calcjy0_qy1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_qy1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy0_qy1(0).v+caljy0_qy1(1).v
              t1.v = caljy0_qy1(2).v+caljy0_qy1(3).v
              t2.v = caljy0_qy1(4).v+caljy0_qy1(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy0_qy1 
       
       
       pure function preload_caljy0_py2() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_py2
              !dir$ attributes forceinline :: preload_calcjy0_py2
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_py2
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v = caljy0_py1(0).v+caljy0_py1(1).v
              t1.v = caljy0_py1(2).v+caljy0_py1(3).v
              t2.v = caljy0_py1(4).v+caljy0_py1(5).v
              t3.v = caljy0_py2(6).v+caljy0_py1(7).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_caljy0_py2 
       
       
       pure function preload_caljy0_qy2() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_qy2
              !dir$ attributes forceinline :: preload_calcjy0_qy2
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_qy2
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy0_qy2(0).v+caljy0_qy2(1).v
              t1.v = caljy0_qy2(2).v+caljy0_qy2(3).v
              t2.v = caljy0_qy2(4).v+caljy0_qy2(5).v+ &
                     caljy0_qy2(6).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy0_qy2
       
       
       pure function preload_caljy0_p0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_p0
              !dir$ attributes forceinline :: preload_calcjy0_p0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_p0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy0_p0(0).v+caljy0_p0(1).v
              t1.v = caljy0_p0(2).v+caljy0_p0(3).v
              t2.v = caljy0_p0(4).v+caljy0_p0(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy0_p0
       
       
       pure function preload_caljy0_q0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_q0
              !dir$ attributes forceinline :: preload_calcjy0_q0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_q0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = caljy0_q0(0).v+caljy0_q0(1).v
              t1.v = caljy0_q0(2).v+caljy0_q0(3).v+ &
                     caljy0_q0(4).v
              summa.v = t0.v+t1.v
       end function preload_caljy0_q0
       
       
       pure function preload_caljy0_p1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_p1
              !dir$ attributes forceinline :: preload_calcjy0_p1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_p1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy0_p1(0).v+caljy0_p1(1).v
              t1.v = caljy0_p1(2).v+caljy0_p1(3).v
              t2.v = caljy0_p1(4).v+caljy0_p1(5).v+ &
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy0_p1
       
       
       pure function preload_caljy0_q1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy0_q1
              !dir$ attributes forceinline :: preload_calcjy0_q1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy0_q1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = caljy0_q1(0).v+caljy0_q1(1).v
              t1.v = caljy0_q1(2).v+caljy0_q1(3).v+ &
                     caljy0_q1(4).v
              summa.v = t0.v+t1.v
       end function preload_caljy0_q1
       
!! =============================================================================================================== //
!!                                  'Saved' arrays preload_caljy1 routines.
!!================================================================================================================ //       
                              
       
       pure function preload_caljy1_plg() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_plg
              !dir$ attributes forceinline :: preload_calcjy1_plg
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_plg
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = caljy1_plg(0).v+caljy1_plg(1).v
              t1.v = caljy1_plg(2).v+caljy1_plg(3).v
              summa.v = t0.v+t1.v
       end function preload_caljy1_plg 
       
       
       pure function preload_caljy1_qlg() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_qlg
              !dir$ attributes forceinline :: preload_calcjy1_qlg
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_qlg
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = caljy1_qlg(0).v+caljy1_qlg(1).v
              t1.v = caljy1_qlg(2).v+caljy1_qlg(3).v
              summa.v = t0.v+t1.v
       end function preload_caljy1_qlg  
       
       
       pure function preload_caljy1_pj0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_pj0
              !dir$ attributes forceinline :: preload_calcjy1_pj0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_pj0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy1_pj0(0).v+caljy1_pj0(1).v
              t1.v = caljy1_pj0(2).v+caljy1_pj0(3).v
              t2.v = caljy1_pj0(4).v+caljy1_pj0(5).v+ &
                     caljy1_pj0(6).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy1_pj0  
       
       
       pure function preload_caljy1_qj0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_qj0
              !dir$ attributes forceinline :: preload_calcjy1_qj0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_qj0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = caljy1_qj0(0).v+caljy1_qj0(1).v
              t1.v = caljy1_qj0(2).v+caljy1_qj0(3).v+ &
                     caljy1_qj0(4).v
              summa.v = t0.v+t1.v
       end function preload_caljy1_qj0 
       
       
       pure function preload_caljy1_pj1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_pj1
              !dir$ attributes forceinline :: preload_calcjy1_pj1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_pj1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v = caljy1_pj1(0).v+caljy1_pj1(1).v
              t1.v = caljy1_pj1(2).v+caljy1_pj1(3).v
              t2.v = caljy1_pj1(4).v+caljy1_pj1(5).v
              t3.v = caljy1_pj1(6).v+caljy1_pj1(7).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_caljy1_pj1   
       
       
       pure function preload_caljy1_qj1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_qj1
              !dir$ attributes forceinline :: preload_calcjy1_qj1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_qj1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy1_qj1(0).v+caljy1_qj1(1).v
              t1.v = caljy1_qj1(2).v+caljy1_qj1(3).v
              t2.v = caljy1_qj1(4).v+caljy1_qj1(5).v+ &
                     caljy1_qj1(6).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy1_qj1   
       
       
       pure function preload_caljy1_py0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_py0
              !dir$ attributes forceinline :: preload_calcjy1_py0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_py0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy1_py0(0).v+caljy1_py0(1).v
              t1.v = caljy1_py0(2).v+caljy1_py0(3).v
              t2.v = caljy1_py0(4).v+caljy1_py0(5).v+ &
                     caljy1_py0(6).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy1_py0  
       
       
       pure function preload_caljy1_qy0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_qy0
              !dir$ attributes forceinline :: preload_calcjy1_qy0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_qy0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy1_qy0(0).v+caljy1_qy0(1).v
              t1.v = caljy1_qy0(2).v+caljy1_qy0(3).v
              t2.v = caljy1_qy0(4).v+caljy1_qy0(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy1_qy0   
       
       
       pure function preload_caljy1_py1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_py1
              !dir$ attributes forceinline :: preload_calcjy1_py1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_py1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v = caljy1_py1(0).v+caljy1_py1(1).v
              t1.v = caljy1_py1(2).v+caljy1_py1(3).v
              t2.v = caljy1_py1(4).v+caljy1_py1(5).v
              t3.v = caljy1_py1(6).v+caljy1_py1(7).v+ &
                     caljy1_py1(8).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_caljy1_py1  
       
       
       pure function preload_caljy1_qy1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_qy1
              !dir$ attributes forceinline :: preload_calcjy1_qy1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_qy1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v = caljy1_qy1(0).v+caljy1_qy1(1).v
              t1.v = caljy1_qy1(2).v+caljy1_qy1(3).v
              t2.v = caljy1_qy1(4).v+caljy1_qy1(5).v
              t3.v = caljy1_qy1(6).v+caljy1_qy1(7).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_caljy1_qy1 
       
       
       pure function preload_caljy1_p0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_p0
              !dir$ attributes forceinline :: preload_calcjy1_p0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_p0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy1_p0(0).v+caljy1_p0(1).v
              t1.v = caljy1_p0(2).v+caljy1_p0(3).v
              t2.v = caljy1_p0(4).v+caljy1_p0(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy1_p0
       
       
       pure function preload_caljy1_q0() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_q0
              !dir$ attributes forceinline :: preload_calcjy1_q0
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_q0
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy1_q0(0).v+caljy1_q0(1).v
              t1.v = caljy1_q0(2).v+caljy1_q0(3).v
              t2.v = caljy1_q0(4).v+caljy1_q0(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy1_q0
       
       
       pure function preload_caljy1_p1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_p1
              !dir$ attributes forceinline :: preload_calcjy1_p1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_p1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy1_p1(0).v+caljy1_p1(1).v
              t1.v = caljy1_p1(2).v+caljy1_p1(3).v
              t2.v = caljy1_p1(4).v+caljy1_p1(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy1_p1
       
       
       pure function preload_caljy1_q1() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_caljy1_q1
              !dir$ attributes forceinline :: preload_calcjy1_q1
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_caljy1_q1
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = caljy1_q1(0).v+caljy1_q1(1).v
              t1.v = caljy1_q1(2).v+caljy1_q1(3).v
              t2.v = caljy1_q1(4).v+caljy1_q1(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_caljy1_q1
       
       
!! =============================================================================================================== //
!!                                  'Saved' arrays preload_calck0 routines.
!!================================================================================================================ //  


      pure function preload_calck0_p() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_p
              !dir$ attributes forceinline :: preload_calck0_p
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_p
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2
              type(YMM4r8_t), automatic :: t0,t1,t2
              t0.v = calck0_p(0).v+calck0_p(1).v
              t1.v = calck0_p(2).v+calck0_p(3).v
              t2.v = calck0_p(4).v+calck0_p(5).v
              summa.v = t0.v+t1.v+t2.v
       end function preload_calck0_p
       
       
       pure function preload_calck0_q() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_q
              !dir$ attributes forceinline :: preload_calck0_q
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_q
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0
              type(YMM4r8_t), automatic :: t0
              t0.v = calck0_p(0).v+calck0_p(1).v
              summa.v = t0.v
       end function preload_calck0_q
       
       
       pure function preload_calck0_f() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_f
              !dir$ attributes forceinline :: preload_calck0_f
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_f
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = calck0_f(0).v+calck0_f(1).v
              t1.v = calck0_f(2).v+calck0_f(3).v
              summa.v = t0.v+t1.v
       end function preload_calck0_f
       
       
       pure function preload_calck0_g() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_g
              !dir$ attributes forceinline :: preload_calck0_g
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_g
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0
              type(YMM4r8_t), automatic :: t0
              t0.v = calck0_g(0).v+calck0_g(1).v+ &
                     calck0_g(2).v
              summa.v = t0.v
       end function preload_calck0_g
       
       
       pure function preload_calck0_pp() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_pp
              !dir$ attributes forceinline :: preload_calck0_pp
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_pp
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3,t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3,t4
              t0.v = calck0_pp(0).v+calck0_pp(1).v
              t1.v = calck0_pp(2).v+calck0_pp(3).v
              t2.v = calck0_pp(4).v+calck0_pp(5).v
              t3.v = calck0_pp(6).v+calck0_pp(7).v
              t4.v = calck0_pp(8).v+calck0_pp(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+t4.v
       end function preload_calck0_pp
       
       
       pure function preload_calck0_qq() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck0_qq
              !dir$ attributes forceinline :: preload_calck0_qq
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck0_qq
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3,t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3,t4
              t0.v = calck0_qq(0).v+calck0_qq(1).v
              t1.v = calck0_qq(2).v+calck0_qq(3).v
              t2.v = calck0_qq(4).v+calck0_qq(5).v
              t3.v = calck0_qq(6).v+calck0_qq(7).v
              t4.v = calck0_qq(8).v+calck0_qq(9).v
              summa.v = t0.v+t1.v+t2.v+t3.v+t4.v
       end function preload_calck0_qq
       
       
!! =============================================================================================================== //
!!                                  'Saved' arrays preload_calck1 routines.
!!================================================================================================================ //   


       pure function preload_calck1_p() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_p
              !dir$ attributes forceinline :: preload_calck1_p
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_p
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = calck1_p(0).v+calck1_p(1).v
              t1.v = calck1_p(2).v+calck1_p(3).v+ &
                     calck1_p(4).v
              summa.v = t0.v+t1.v
       end function preload_calck1_p    
       
       
       pure function preload_calck1_q() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_q
              !dir$ attributes forceinline :: preload_calck1_q
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_q
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0
              type(YMM4r8_t), automatic :: t0
              t0.v = calck1_p(0).v+calck1_p(1).v+ &
                     calck1_p(2).v
              summa.v = t0.v
       end function preload_calck1_q 
       
       
       pure function preload_calck1_f() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_f
              !dir$ attributes forceinline :: preload_calck1_f
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_f
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1
              type(YMM4r8_t), automatic :: t0,t1
              t0.v = calck1_f(0).v+calck1_f(1).v
              t1.v = calck1_f(2).v+calck1_f(3).v+ &
                     calck1_f(4).v
              summa.v = t0.v+t1.v
       end function preload_calck1_f
       
       
       pure function preload_calck1_g() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_g
              !dir$ attributes forceinline :: preload_calck1_g
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_g
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0
              type(YMM4r8_t), automatic :: t0
              t0.v = calck1_f(0).v+calck1_f(1).v+ &
                     calck1_f(2).v
              summa.v = t0.v
       end function preload_calck1_g
       
       
       pure function preload_calck1_pp() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_pp
              !dir$ attributes forceinline :: preload_calck1_pp
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_pp
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3,t4
              type(YMM4r8_t), automatic :: t0,t1,t2,t3,t4
              t0.v = calck1_pp(0).v+calck1_pp(1).v
              t1.v = calck1_pp(2).v+calck1_pp(3).v
              t2.v = calck1_pp(4).v+calck1_pp(5).v
              t3.v = calck1_pp(6).v+calck1_pp(7).v
              t4.v = calck1_pp(8).v+calck1_pp(9).v+ &
                     calck1_pp(10).v
              summa.v = t0.v+t1.v+t2.v+t3.v+t4.v
       end function preload_calck1_pp
       
       
       pure function preload_calck1_qq() result(summa)
       
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: preload_calck1_qq
              !dir$ attributes forceinline :: preload_calck1_qq
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: preload_calck1_qq
              type(YMM4r8_t) :: summa
              !dir$ attributes align : 64 :: t0,t1,t2,t3
              type(YMM4r8_t), automatic :: t0,t1,t2,t3
              t0.v = calck1_qq(0).v+calck1_qq(1).v
              t1.v = calck1_qq(2).v+calck1_qq(3).v
              t2.v = calck1_qq(4).v+calck1_qq(5).v
              t3.v = calck1_qq(6).v+calck1_qq(7).v+ &
                     calck1_qq(8).v
              summa.v = t0.v+t1.v+t2.v+t3.v
       end function preload_calck1_qq
       
       
#if 0
/*
               !*****************************************************************************80
!
!! BESEI0 evaluates the exponentially scaled Bessel I0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the modified Bessel
!    function of the first kind of order zero multiplied by EXP(-ABS(X)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESEI0, the value of the function.
!
               
*/

#endif


          pure function besei0_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besei0_ymm4r8
              !dir$ attributes forceinline :: besei0_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besei0_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 2
              val  = calci0_ymm4r8(x,jint)
          end function besei0_ymm4r8
          

#if 0
 /*
!*****************************************************************************80
!
!! BESEI1 evaluates the exponentially scaled Bessel I1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the first kind of order one
!    multiplied by EXP(-ABS(X)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESEI1, the value of the function.
*/	         

#endif      


          pure function besei1_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besei1_ymm4r8
              !dir$ attributes forceinline :: besei1_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besei1_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 2
              val  = calci1_ymm4r8(x,jint) 
          end function besei1_ymm4r8
          

#if 0         
/*
   !*****************************************************************************80
!
!! BESEK0 evaluates the exponentially scaled Bessel K0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order zero
!    multiplied by the exponential function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!    0 < X.
!
!    Output, real ( kind = 8 ) BESK0, the value of the function.
*/   
#endif


           pure function besek0_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besek0_ymm4r8
              !dir$ attributes forceinline :: besek0_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besek0_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 2
              val  = calck0_ymm4r8(x,jint) 
          end function besek0_ymm4r8    
          

#if 0          
/*
!*****************************************************************************80
!
!! BESEK1 evaluates the exponentially scaled Bessel K1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order one
!    multiplied by the exponential function, for arguments
!    XLEAST <= ARG <= XMAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESEK1, the value of the function.
*/	             
#endif


          pure function besek1_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besek1_ymm4r8
              !dir$ attributes forceinline :: besek1_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besek1_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 2
              val  = calck1_ymm4r8(x,jint) 
          end function besek1_ymm4r8   
          
#if 0          
/*
     !*****************************************************************************80
!
!! BESI0 evaluates the Bessel I0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for
!    modified Bessel functions of the first kind of order zero for
!    arguments ABS(ARG) <= XMAX.
!
!    See comments heading CALCI0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESI0, the value of the function.
*/ 
#endif


          pure function besi0_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besi0_ymm4r8
              !dir$ attributes forceinline :: besi0_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besi0_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 1
              val  = calci0_ymm4r8(x,jint) 
          end function besi0_ymm4r8  
          
          
#if 0
/*
    !*****************************************************************************80
!
!! BESI1 evaluates the Bessel I1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for
!    modified Bessel functions of the first kind of order one for
!    arguments ABS(ARG) <= XMAX.
!
!    See comments heading CALCI1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESI1, the value of the function.        
*/
#endif  


           pure function besi1_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besi1_ymm4r8
              !dir$ attributes forceinline :: besi1_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besi1_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 1
              val  = calci1_ymm4r8(x,jint) 
          end function besi1_ymm4r8   
          
          
#if 0
/*
    *****************************************************************************80
!
!! BESJ0 evaluates the Bessel J0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for Bessel functions
!    of the first kind of order zero for arguments  |X| <= XMAX
!
!    See comments heading CALJY0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESJ0, the value of the function.       
*/
#endif    


         pure function besj0_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besj0_ymm4r8
              !dir$ attributes forceinline :: besj0_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besj0_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 0
              val  = calcjy0_ymm4r8(x,jint) 
          end function besj0_ymm4r8   
          
          
#if 0
   /*
*****************************************************************************80
!
!! BESJ1 evaluates the Bessel J1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for Bessel functions
!    of the first kind of order zero for arguments  |X| <= XMAX
!
!    See comments heading CALJY1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESJ1, the value of the function.
*/	 

#endif

   
          pure function besj1_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besj1_ymm4r8
              !dir$ attributes forceinline :: besj1_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besj1_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 0
              val  = calcjy1_ymm4r8(x,jint) 
          end function besj1_ymm4r8   
          
          
#if 0
/*
    !*****************************************************************************80
!
!! BESK0 evaluates the Bessel K0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order zero
!    for arguments 0.0 < ARG <= XMAX.
!
!    See comments heading CALCK0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESK0, the value of the function.
*/
#endif   


         pure function besk0_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besk0_ymm4r8
              !dir$ attributes forceinline :: besk0_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besk0_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 1
              val  = calck0_ymm4r8(x,jint) 
          end function besk0_ymm4r8   
          
          
#if 0
/*
      *****************************************************************************80
!
!! BESK1 evaluates the Bessel K1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order one
!    for arguments XLEAST <= ARG <= XMAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESK1, the value of the function.        
*/  
#endif


         pure function besk1_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besk1_ymm4r8
              !dir$ attributes forceinline :: besk1_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besk1_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 1
              val  = calck1_ymm4r8(x,jint) 
          end function besk1_ymm4r8   
          
          
#if 0
/*
  !*****************************************************************************80
!
!! BESY0 evaluates the Bessel Y0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for Bessel functions
!    of the second kind of order zero for arguments 0 < X <= XMAX.
!
!    See comments heading CALJY0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESY0, the value of the function.
*/	      
#endif


          pure function besy0_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besy0_ymm4r8
              !dir$ attributes forceinline :: besy0_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besy0_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 1
              val  = caly0_ymm4r8(x,jint) 
          end function besy0_ymm4r8   
          
          
#if 0
/*
    !*****************************************************************************80
!
!! BESY1 evaluates the Bessel Y1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for Bessel functions
!    of the second kind of order zero for arguments 0 < X <= XMAX.
!
!    See comments heading CALJY1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESY1, the value of the function.     

*/    
#endif


          pure function besy1_ymm4r8(x) result(val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: besy1_ymm4r8
              !dir$ attributes forceinline :: besy1_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: besy1_ymm4r8
              type(YMM4r8_t),   intent(in) :: x
              type(YMM4r8_t)  :: val
              integer(kind=i4), automatic :: jint
              jint = 1
              val  = caly1_ymm4r8(x,jint) 
          end function besy1_ymm4r8   

!! /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !!

#if 0
        !*****************************************************************************80
!
!! CALCEI computes various exponential integrals.
!
!  Discussion:
!
!    This routine computes the exponential integrals Ei(x),
!    E1(x), and  exp(-x)*Ei(x) for real arguments x where
!
!           integral (from t=-oo to t=x) (exp(t)/t),  x > 0,
!    Ei(x) =
!          -integral (from t=-x to t=+oo) (exp(t)/t),  x < 0,
!
!    and where the first integral is a principal value integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Henry Thacher,
!    Rational Chebyshev Approximations for the Exponential
!    Integral E1(x),
!    Mathematics of Computation,
!    Volume 22, Number 103, July 1968, pages 641-649.
!
!    William Cody, Henry Thacher,
!    Chebyshev Approximations for the Exponential
!    Integral Ei(x),
!    Mathematics of Computation,
!    Volume 23, Number 106, April 1969, pages 289-303.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  The argument must not
!    be zero.  If JINT = 2, then the argument must be strictly positive.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = EI ( ARG );
!    2, RESULT = EONE ( ARG );
!    3, RESULT = EXPEI ( ARG ).
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, Ei(x);
!    2, -Ei(-x);
!    3, exp(-x)*Ei(x).
!
#endif   


         subroutine calcei_ymm4r8(arg,val,jint) 
           
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: calcei_ymm4r8
              !dir$ attributes forceinline :: calcei_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: calcei_ymm4r8  
               type(YMM4r8_t),   intent(in)   :: arg
               type(YMM4r8_t),   intent(out)  :: val
               integer(kind=i4), intent(in)   :: jint
               
               !dir$ attributes align : 64 :: q
               !dir$ attributes align : 64 :: qlq
               !dir$ attributes align : 64 :: qx
               !dir$ attributes align : 64 :: px
               type(YMM4r8_t), dimension(0:9), automatic  :: q
               type(YMM4r8_t), dimension(0:9), automatic  :: qlq
               type(YMM4r8_t), dimension(0:9), automatic  :: qx
               type(YMM4r8_t), dimension(0:9), automatic  :: px
               !dir$ attributes align : 64 :: zero
               !dir$ attributes align : 64 :: p037
               !dir$ attributes align : 64 :: half
               !dir$ attributes align : 64 :: one
               !dir$ attributes align : 64 :: two
               !dir$ attributes align : 64 :: three
               !dir$ attributes align : 64 :: four
               !dir$ attributes align : 64 :: six
               !dir$ attributes align : 64 :: twlve
               !dir$ attributes align : 64 :: two4
               !dir$ attributes align : 64 :: frty
               !dir$ attributes align : 64 :: exp40
               !dir$ attributes align : 64 :: x01
               !dir$ attributes align : 64 :: x11
               !dir$ attributes align : 64 :: x02
               !dir$ attributes align : 64 :: x0
               !dir$ attributes align : 64 :: xinf
               !dir$ attributes align : 64 :: xmax
               !dir$ attributes align : 64 :: xbig
               type(YMM4r8_t), parameter :: zero = YMM4r8_t(0.0e+0_dp)
               type(YMM4r8_t), parameter :: p037 = YMM4r8_t(0.037e+0_dp)
               type(YMM4r8_t), parameter :: half = YMM4r8_t(0.5e+0_dp)
               type(YMM4r8_t), parameter :: one  = YMM4r8_t(1.0e+0_dp)
               type(YMM4r8_t), parameter :: two  = YMM4r8_t(2.0e+0_dp)
               type(YMM4r8_t), parameter :: three= YMM4r8_t(3.0e+0_dp)
               type(YMM4r8_t), parameter :: four = YMM4r8_t(4.0e+0_dp)
               type(YMM4r8_t), parameter :: six  = YMM4r8_t(6.0e+0_dp)
               type(YMM4r8_t), parameter :: twlve= YMM4r8_t(12.0e+0_dp)
               type(YMM4r8_t), parameter :: two4 = YMM4r8_t(24.0e+0_dp)
               type(YMM4r8_t), parameter :: frty = YMM4r8_t(40.0e+0_dp)
               type(YMM4r8_t), parameter :: exp40= YMM4r8_t(2.3538526683701998541e+17)
               type(YMM4r8_t), parameter :: x01  = YMM4r8_t(381.5e+0_dp)
               type(YMM4r8_t), parameter :: x11  = YMM4r8_t(1024.0e+0_dp)
               type(YMM4r8_t), parameter :: x02  = YMM4r8_t(-5.1182968633365538008e-5_dp)
               type(YMM4r8_t), parameter :: x0   = YMM4r8_t(3.7250741078136663466e-1_dp)
               type(YMM4r8_t), parameter :: xinf = YMM4r8_t(1.79e+308_dp)
               type(YMM4r8_t), parameter :: xmax = YMM4r8_t(716.351e+0_dp)
               type(YMM4r8_t), parameter :: xbig = YMM4r8_t(701.84e+0_dp)
               !dir$ attributes align : 64 :: ei
               !dir$ attributes align : 64 :: frac
               !dir$ attributes align : 64 :: res
               !dir$ attributes align : 64 :: sump
               !dir$ attributes align : 64 :: sumq
               !dir$ attributes align : 64 :: t
               !dir$ attributes align : 64 :: w
               !dir$ attributes align : 64 :: x
               !dir$ attributes align : 64 :: mx0
               !dir$ attributes align : 64 :: y
               !dir$ attributes align : 64 :: ysq
               !dir$ attributes align : 64 :: t0
               !dir$ attributes align : 64 :: t1
               type(YMM4r8_t), automatic :: ei,frac
               type(YMM4r8_t), automatic :: res
               type(YMM4r8_t), automatic :: sump,sumq
               type(YMM4r8_t), automatic :: t,w
               type(YMM4r8_t), automatic :: x,mx0
               type(YMM4r8_t), automatic :: y,ysq
               type(YMM4r8_t), automatic :: t0,t1
               type(Mask4_t),  automatic :: msk1
               type(Mask4_t),  automatic :: msk2
               type(Mask4_t),  automatic :: msk3
               type(Mask4_t),  automatic :: msk4
               type(Mask4_t),  automatic :: msk5
               type(Mask4_t),  automatic :: msk6
               type(Mask4_t),  automatic :: msk7
               type(Mask4_t),  automatic :: msk8
               type(Mask4_t),  automatic :: msk9
               type(Mask4_t),  automatic :: msk10
#if (GMS_EXPLICIT_VECTORIZE) == 1
               integer(kind=i4) :: j
#endif   
#if (GMS_EXPLICIT_VECTORIZE) == 1
           
               !dir$ loop_count(4)
               !dir$ vector aligned
               !dir$ vector vectorlength(8)
               !dir$ vector always
               do j=0,3
                   x.v(j)    = arg.v(j)
                   msk1.m(j) = (x.v(j)==zero.v(j))
                   msk2.m(j) = (x.v(j)<zero.v(j))
                   msk6.m(j) = (x.v(j)<six.v(j))
                   msk8.m(j) = (x.v(j)<twlve.v(j))
                   msk9.m(j) = (x.v(j)<=two4.v(j))
                   if(all(msk1.m(j))) then
               
                       ei.v(j) = -xinf.v(j)
                       if(jint==2) ei.v(j) = -ei.v(j)
                   ! /*
	           !             !
                   !              !  Calculate EI for negative argument or for E1.
                   !             !   
	           !          */
	           else if(all(msk2.m(j)).or.jint==2) then
	       
	                  y.v(j)    = abs(x.v(j))
	                  msk3.m(j) = (y.v(j)<one.v(j))
	                  msk4.m(j) = (y.v(j)<=four.v(j))
	                  msk5.m(j) = (xbig.v(j)<y.v(j))
	                  if(all(msk3.m(j))) then
	           
	                     sump.v(j) = calcei_a(6).v(j)*y.v(j)+calcei_a(0).v(j))
	                     sumq.v(j) = y.v(j)+calcei_b(0).v(j)
	                     sump.v(j) = sump.v(j)*y.v(j)+calcei_a(1).v(j)
	                     sumq.v(j) = sumq.v(j)*y.v(j)+calcei_b(1).v(j)
	                     sump.v(j) = sump.v(j)*y.v(j)+calcei_a(2).v(j)
	                     sumq.v(j) = sumq.v(j)*y.v(j)+calcei_b(2).v(j)
	                     sump.v(j) = sump.v(j)*y.v(j)+calcei_a(3).v(j)
	                     sumq.v(j) = sumq.v(j)*y.v(j)+calcei_b(3).v(j)
	                     sump.v(j) = sump.v(j)*y.v(j)+calcei_a(4).v(j)
	                     sumq.v(j) = sumq.v(j)*y.v(j)+calcei_b(4).v(j)
	                     sump.v(j) = sump.v(j)*y.v(j)+calcei_a(5).v(j)
	                     sumq.v(j) = sumq.v(j)*y.v(j)+calcei_b(5).v(j)
	                     ei.v(j)   = log(y.v(j))-(sump.v(j)/sumq.v(j))
	                     if(jint==3) ei.v(j) = ei.v(j)*exp(y.v(j))
	              
	                 else if(all(msk4.m(j))) then
	              
	                     w.v(j)    = one.v(j)/y.v(j)
	                     sump.v(j) = calcei_c(0).v(j)
	                     sumq.v(j) = calcei_d(0).v(j)
	                     sump.v(j) = sump.v(j)*w.v(j)+calcei_c(1).v(j)
	                     sumq.v(j) = sumq.v(j)*w.v(j)+calcei_d(1).v(j)
	                     sump.v(j) = sump.v(j)*w.v(j)+calcei_c(2).v(j)
	                     sumq.v(j) = sumq.v(j)*w.v(j)+calcei_d(2).v(j)
	                     sump.v(j) = sump.v(j)*w.v(j)+calcei_c(3).v(j)
	                     sumq.v(j) = sumq.v(j)*w.v(j)+calcei_d(3).v(j)
	                     sump.v(j) = sump.v(j)*w.v(j)+calcei_c(4).v(j)
	                     sumq.v(j) = sumq.v(j)*w.v(j)+calcei_d(4).v(j)
	                     sump.v(j) = sump.v(j)*w.v(j)+calcei_c(5).v(j)
	                     sumq.v(j) = sumq.v(j)*w.v(j)+calcei_d(5).v(j)
	                     sump.v(j) = sump.v(j)*w.v(j)+calcei_c(6).v(j)
	                     sumq.v(j) = sumq.v(j)*w.v(j)+calcei_d(6).v(j)
	                     sump.v(j) = sump.v(j)*w.v(j)+calcei_c(7).v(j)
	                     sumq.v(j) = sumq.v(j)*w.v(j)+calcei_d(7).v(j)
	                     sump.v(j) = sump.v(j)*w.v(j)+calcei_c(8).v(j)
	                     sumq.v(j) = sumq.v(j)*w.v(j)+calcei_d(8).v(j)
	                     ei.v(j)   = -sump.v(j)/sumq.v(j)
	              
	                     if(jint/=3) ei.v(j) = ei.v(j)*exp(-y.v(j)) 
	              
	               else 
	               
	                     if(all(msk5.m(j)).and.jint<3) then
	              
	                         ei.v(j) = zero.v(j)
	               
	                    else
	               
	                         w.v(j)    = one.v(j)/y.v(j)
	                         sump.v(j) = calcei_e(0).v(j)
	                         sumq.v(j) = calcei_f(0).v(j)
	                         sump.v(j) = sump.v(j)*w.v(j)+calcei_e(1).v(j)
	                         sumq.v(j) = sumq.v(j)*w.v(j)+calcei_f(1).v(j)
	                         sump.v(j) = sump.v(j)*w.v(j)+calcei_e(2).v(j)
	                         sumq.v(j) = sumq.v(j)*w.v(j)+calcei_f(2).v(j)
	                         sump.v(j) = sump.v(j)*w.v(j)+calcei_e(3).v(j)
	                         sumq.v(j) = sumq.v(j)*w.v(j)+calcei_f(3).v(j) 
	                         sump.v(j) = sump.v(j)*w.v(j)+calcei_e(4).v(j)
	                         sumq.v(j) = sumq.v(j)*w.v(j)+calcei_f(4).v(j)
	                         sump.v(j) = sump.v(j)*w.v(j)+calcei_e(5).v(j)
	                         sumq.v(j) = sumq.v(j)*w.v(j)+calcei_f(5).v(j)
	                         sump.v(j) = sump.v(j)*w.v(j)+calcei_e(6).v(j)
	                         sumq.v(j) = sumq.v(j)*w.v(j)+calcei_f(6).v(j)
	                         sump.v(j) = sump.v(j)*w.v(j)+calcei_e(7).v(j)
	                         sumq.v(j) = sumq.v(j)*w.v(j)+calcei_f(7).v(j)
	                         sump.v(j) = sump.v(j)*w.v(j)+calcei_e(8).v(j)
	                         sumq.v(j) = sumq.v(j)*w.v(j)+calcei_f(8).v(j)
	                         sump.v(j) = sump.v(j)*w.v(j)+calcei_e(9).v(j)
	                         sumq.v(j) = sumq.v(j)*w.v(j)+calcei_f(9).v(j)
	                         t0.v(j)   = sump.v(j)/sumq.v(j)
	                         t1.v(j)   = one.v(j)-w.v(j)
	                         ei.v(j)   = -w.v(j)*t0.v(j)*t1.v(j)
	               
	                         if(jint/=3) ei.v(j) = -y.v(j)*ei.v(j)
	               
	                 end if
	           
	             end if
	       
	             if(jint==2) ei.v(j) = -ei.v(j)
	              !    /*
	              !                  !
                      !                  !  To improve conditioning, rational approximations are expressed
                      !                  !  in terms of Chebyshev polynomials for 0 <= X < 6, and in
                      !                  !  continued fraction form for larger X.
                      !                  !
	              !               */
	       else if(all(msk6.m(j))) then
	             
	             t.v(j)     = x.v(j)+x.v(j)
	             t.v(j)     = (t.v(j)/three.v(j))-two.v(j)
	             px(0).v(j) = zero.v(j)
	             qx(0).v(j) = zero.v(j)
	             px(1).v(j) = p(0).v(j)
	             qx(1).v(j) = q(0).v(j)
	             px(2).v(j) = t.v(j)*px(1).v(j)-px(0).v(j)+p(1).v(j)
	             qx(2).v(j) = t.v(j)*qx(1).v(j)-qx(0).v(j)+q(1).v(j)
	             px(3).v(j) = t.v(j)*px(2).v(j)-px(1).v(j)+p(2).v(j)
	             qx(3).v(j) = t.v(j)*qx(2).v(j)-qx(1).v(j)+q(2).v(j)
	             px(4).v(j) = t.v(j)*px(3).v(j)-px(2).v(j)+p(3).v(j)
	             qx(4).v(j) = t.v(j)*qx(3).v(j)-qx(2).v(j)+q(3).v(j)
	             px(5).v(j) = t.v(j)*px(4).v(j)-px(3).v(j)+p(4).v(j)
	             qx(5).v(j) = t.v(j)*qx(4).v(j)-qx(3).v(j)+q(4).v(j)
	             px(6).v(j) = t.v(j)*px(5).v(j)-px(4).v(j)+p(5).v(j)
	             qx(6).v(j) = t.v(j)*qx(5).v(j)-qx(4).v(j)+q(5).v(j)
	             px(7).vv(j) = t.v(j)*px(6).v(j)-px(5).v(j)+p(6).v(j)
	             qx(7).v(j) = t.v*qx(6).v-qx(5).v+q(6).v(j)
	             px(8).v(j) = t.v(j)*px(7).v(j)-px(6).v(j)+p(7).v(j)
	             qx(8).v(j) = t.v(j)*qx(7).v(j)-qx(6).v(j)+q(7).v(j)
	             t0.v(j)    = half.v(j)*t.v(j)
	             sump.v(j)  = t0.v(j)*px(9).v(j)-px(8).v(j)+p(9).v(j)
	             sumq.v(j)  = t0.v(j)*qx(9).v(j)-qx(8).v(j)+q(9).v(j)
	             frac.v(j)  = sump.v(j)/sumq.v(j)
	             t0.v(j)    = x.v(j)-x01.v(j)/x11.v(j)
	             xmx0.v(j)  = t0.v(j)-x02.v(j)
	             msk7.m(j)  = (p037.v(j)<=abs(xmx0.v(j))
	             
	             if(all(msk7.m(j))) then
	                   t0.v(j) = x.v(j)/x0.v(j)
	                   ei.v(j) = frac.v(j)*xmx0.v(j)+log(t0.v(j))
	                   if(jint==3) ei.v(j) = exp(-x.v(j))*ei.v(j)
	             else
	                 !  //Special approximation to ln(X/X0) for X close to X0. 
	                 y.v(j)    = xmx0.v(j)/(x.v(j)+x0.v(j))
	                 ysq.v(j)  = y.v(j)*y.v(j)
	                 sump.v(j) = calcei_plg(0).v(j)
	                 sumq.v(j) = ysq.v(j)+calcei_qlg(0).v(j)
	                 sump.v(j) = sump.v(j)*ysq.v(j)+calcei_plg(1).v(j)
	                 sumq.v(j) = sumq.v(j)*ysq.v(j)+calcei_qlg(1).v(j)
	                 sump.v(j) = sump.v(j)*ysq.v(j)+calcei_plg(2).v(j)
	                 sumq.v(j) = sumq.v(j)*ysq.v(j)+calcei_qlg(2).v(j)
	                 sump.v(j) = sump.v(j)*ysq.v(j)+calcei_plg(3).v(j)
	                 sumq.v(j) = sumq.v(j)*ysq.v(j)+calcei_qlg(3).v(j)
	                 t0.v(j)   = sumq.v(j)*(x.v(j)+x0.v(j))+frac.v(j)
	                 ei.v(j)   = (sump.v(j)/t0.v(j))*xmx0.v(j)
	                 
	                 if(jint==3) ei.v(j) = exp(-x.v(j))*ei.v(j) 
	             end if
	             
	       else if(all(msk8.m(j))) then
	       
	              frac.v(j)= zero.v(j)
	              frac.v(j)= calcei_s(0).v(j)/(calcei_r(0).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_s(1).v(j)/(calcei_r(1).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_s(2).v(j)/(calcei_r(2).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_s(3).v(j)/(calcei_r(3).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_s(4).v(j)/(calcei_r(4).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_s(5).v(j)/(calcei_r(5).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_s(6).v(j)/(calcei_r(6).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_s(7).v(j)/(calcei_r(7).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_s(8).v(j)/(calcei_r(8).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              ei.v(j)  = (calcei_r(9).v(j)+frac.v(j))/x.v(j)
	              
	              if(jint==3) ei.v(j) = ei.v(j)*exp(x.v(j))
	              
	        else if(all(msk9.m(j))) then
	        
	              frac.v(j) = zero.v(j)
	              frac.v(j)= calcei_q1(0).v(j)/(calcei_p1(0).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_q1(1).v(j)/(calcei_p1(1).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_q1(2).v(j)/(calcei_p1(2).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_q1(3).v(j)/(calcei_p1(3).v(j)+ &
	                                          x.v+frac.v(j))
	              frac.v(j)= calcei_q1(4).v(j)/(calcei_p1(4).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_q1(5).v(j)/(calcei_p1(5).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_q1(6).v(j)/(calcei_p1(6).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_q1(7).v(j)/(calcei_p1(7).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              frac.v(j)= calcei_q1(8).v(j)/(calcei_p1(8).v(j)+ &
	                                          x.v(j)+frac.v(j))
	              ei.v(j)  = (calcei_p1(9).v(j)+frac.v(j))/x.v(j)
	              
	              if(jint/=3) ei.v(j) = ei.v(j)*exp(x.v(j))
	                                     
	        else
	             
	               msk10.m(j) = (xmax.v(j)<=x.v(j))
	               if(all(msk10.m).and.jint==3) then
	                   ei.v(j) = xinf.v(j)
	               else
	                   y.v(j)    = one.v(j)/x.v(j)
	                   frac.v(j) = zero.v(j)
	                   frac.v(j)= calcei_q2(0).v(j)/(calcei_p2(0).v(j)+ &
	                                               x.v+frac.v(j))
	                   frac.v(j)= calcei_q2(1).v(j)/(calcei_p2(1).v(j)+ &
	                                          x.v(j)+frac.v(j))
	                   frac.v(j)= calcei_q2(2).v(j)/(calcei_p2(2).v(j)+ &
	                                          x.v(j)+frac.v(j))
	                   frac.v(j)= calcei_q2(3).v(j)/(calcei_p2(3).v(j)+ &
	                                          x.v(j)+frac.v(j))
	                   frac.v(j)= calcei_q2(4).v(j)/(calcei_p2(4).v(j)+ &
	                                          x.v+frac.v(j))
	                   frac.v(j)= calcei_q2(5).v(j)/(calcei_p2(5).v(j)+ &
	                                          x.v+frac.v(j))
	                   frac.v(j)= calcei_q2(6).v(j)/(calcei_p2(6).v(j)+ &
	                                          x.v+frac.v(j))
	                   frac.v(j)= calcei_q2(7).v(j)/(calcei_p2(7).v(j)+ &
	                                          x.v(j)+frac.v(j))
	                   frac.v(j)= calcei_q2(8).v(j)/(calcei_p2(8).v(j)+ &
	                                          x.v(j)+frac.v(j))   
	                   frac.v(j)= calcei_p2(9).v(j)+frac.v(j)
	                   ei.v(j)  = frac.v(j)*y.v(j)*y.v(j)+y.v(j)
	                   
	                   if(jint/=3) then
	                      msk10.m(j) = (x.v(j)<=(xmax.v(j)-two4.v(j)))
	                      if(all(msk10.m(j))) then
	                          ei.v(j) = ei.v(j)*exp(x.v(j))
	                      else
	                          t0.v(j) = x.v(j)-frty.v(j)
	                          ei.v(j) = ei.v(j)*t0.v(j)*exp40.v(j)
	                      end if
	                   end if
	               end if
	           end if
	           val.v(j) = ei.v(j) 
               end do
#else
               x.v    = arg.v
               msk1.m = (x.v==zero.v)
               msk2.m = (x.v<zero.v)
               msk6.m = (x.v<six.v)
               msk8.m = (x.v<twlve.v)
               msk9.m = (x.v<=two4.v)
               if(all(msk1.m)) then
               
                   ei.v = -xinf.v
                   if(jint==2) ei.v = -ei.v
                   ! /*
	           !             !
                   !              !  Calculate EI for negative argument or for E1.
                   !             !   
	           !          */
	       else if(all(msk2.m).or.jint==2) then
	       
	             y.v    = abs(x.v)
	             msk3.m = (y.v<one.v)
	             msk4.m = (y.v<=four.v)
	             msk5.m = (xbig.v<y.v)
	             if(all(msk3.m)) then
	           
	                sump.v = calcei_a(6).v*y.v+calcei_a(0).v
	                sumq.v = y.v+calcei_b(0).v
	                sump.v = sump.v*y.v+calcei_a(1).v
	                sumq.v = sumq.v*y.v+calcei_b(1).v
	                sump.v = sump.v*y.v+calcei_a(2).v
	                sumq.v = sumq.v*y.v+calcei_b(2).v
	                sump.v = sump.v*y.v+calcei_a(3).v
	                sumq.v = sumq.v*y.v+calcei_b(3).v
	                sump.v = sump.v*y.v+calcei_a(4).v
	                sumq.v = sumq.v*y.v+calcei_b(4).v
	                sump.v = sump.v*y.v+calcei_a(5).v
	                sumq.v = sumq.v*y.v+calcei_b(5).v
	                ei.v   = log(y.v)-(sump.v/sumq.v)
	                if(jint==3) ei.v = ei.v*exp(y.v)
	              
	             else if(all(msk4.m)) then
	              
	                w.v    = one.v/y.v
	                sump.v = calcei_c(0).v
	                sumq.v = calcei_d(0).v
	                sump.v = sump.v*w.v+calcei_c(1).v
	                sumq.v = sumq.v*w.v+calcei_d(1).v
	                sump.v = sump.v*w.v+calcei_c(2).v
	                sumq.v = sumq.v*w.v+calcei_d(2).v
	                sump.v = sump.v*w.v+calcei_c(3).v
	                sumq.v = sumq.v*w.v+calcei_d(3).v
	                sump.v = sump.v*w.v+calcei_c(4).v
	                sumq.v = sumq.v*w.v+calcei_d(4).v
	                sump.v = sump.v*w.v+calcei_c(5).v
	                sumq.v = sumq.v*w.v+calcei_d(5).v
	                sump.v = sump.v*w.v+calcei_c(6).v
	                sumq.v = sumq.v*w.v+calcei_d(6).v
	                sump.v = sump.v*w.v+calcei_c(7).v
	                sumq.v = sumq.v*w.v+calcei_d(7).v
	                sump.v = sump.v*w.v+calcei_c(8).v
	                sumq.v = sumq.v*w.v+calcei_d(8).v
	                ei.v   = -sump.v/sumq.v
	              
	                if(jint/=3) ei.v = ei.v*exp(-y.v) 
	              
	             else 
	               
	                if(all(msk5.m).and.jint<3) then
	              
	                    ei.v = zero.v
	               
	                else
	               
	                    w.v    = one.v/y.v
	                    sump.v = calcei_e(0).v
	                    sumq.v = calcei_f(0).v
	                    sump.v = sump.v*w.v+calcei_e(1).v
	                    sumq.v = sumq.v*w.v+calcei_f(1).v
	                    sump.v = sump.v*w.v+calcei_e(2).v
	                    sumq.v = sumq.v*w.v+calcei_f(2).v
	                    sump.v = sump.v*w.v+calcei_e(3).v
	                    sumq.v = sumq.v*w.v+calcei_f(3).v 
	                    sump.v = sump.v*w.v+calcei_e(4).v
	                    sumq.v = sumq.v*w.v+calcei_f(4).v
	                    sump.v = sump.v*w.v+calcei_e(5).v
	                    sumq.v = sumq.v*w.v+calcei_f(5).v
	                    sump.v = sump.v*w.v+calcei_e(6).v
	                    sumq.v = sumq.v*w.v+calcei_f(6).v
	                    sump.v = sump.v*w.v+calcei_e(7).v
	                    sumq.v = sumq.v*w.v+calcei_f(7).v
	                    sump.v = sump.v*w.v+calcei_e(8).v
	                    sumq.v = sumq.v*w.v+calcei_f(8).v
	                    sump.v = sump.v*w.v+calcei_e(9).v
	                    sumq.v = sumq.v*w.v+calcei_f(9).v
	                    t0.v   = sump.v/sumq.v
	                    t1.v   = one.v-w.v
	                    ei.v   = -w.v*t0.v*t1.v
	               
	                    if(jint/=3) ei.v = -y.v*ei.v
	               
	                end if
	           
	             end if
	       
	             if(jint==2) ei.v = -ei.v
	              !    /*
	              !                  !
                      !                  !  To improve conditioning, rational approximations are expressed
                      !                  !  in terms of Chebyshev polynomials for 0 <= X < 6, and in
                      !                  !  continued fraction form for larger X.
                      !                  !
	              !               */
	       else if(all(msk6.m)) then
	             
	             t.v     = x.v+x.v
	             t.v     = (t.v/three.v)-two.v
	             px(0).v = zero.v
	             qx(0).v = zero.v
	             px(1).v = p(0).v
	             qx(1).v = q(0).v
	             px(2).v = t.v*px(1).v-px(0).v+p(1).v
	             qx(2).v = t.v*qx(1).v-qx(0).v+q(1).v
	             px(3).v = t.v*px(2).v-px(1).v+p(2).v
	             qx(3).v = t.v*qx(2).v-qx(1).v+q(2).v
	             px(4).v = t.v*px(3).v-px(2).v+p(3).v
	             qx(4).v = t.v*qx(3).v-qx(2).v+q(3).v
	             px(5).v = t.v*px(4).v-px(3).v+p(4).v
	             qx(5).v = t.v*qx(4).v-qx(3).v+q(4).v
	             px(6).v = t.v*px(5).v-px(4).v+p(5).v
	             qx(6).v = t.v*qx(5).v-qx(4).v+q(5).v
	             px(7).v = t.v*px(6).v-px(5).v+p(6).v
	             qx(7).v = t.v*qx(6).v-qx(5).v+q(6).v
	             px(8).v = t.v*px(7).v-px(6).v+p(7).v
	             qx(8).v = t.v*qx(7).v-qx(6).v+q(7).v
	             t0.v    = half.v*t.v
	             sump.v  = t0.v*px(9).v-px(8).v+p(9).v
	             sumq.v  = t0.v*qx(9).v-qx(8).v+q(9).v
	             frac.v  = sump.v/sumq.v
	             t0.v    = x.v-x01.v/x11.v
	             xmx0.v  = t0.v-x02.v
	             msk7.m  = (p037.v<=abs(xmx0.v)
	             
	             if(all(msk7.m)) then
	                   t0.v = x.v/x0.v
	                   ei.v = frac.v*xmx0.v+log(t0.v)
	                   if(jint==3) ei.v = exp(-x.v)*ei.v
	             else
	                 !  //Special approximation to ln(X/X0) for X close to X0. 
	                 y.v    = xmx0.v/(x.v+x0.v)
	                 ysq.v  = y.v*y.v
	                 sump.v = calcei_plg(0).v
	                 sumq.v = ysq.v+calcei_qlg(0).v
	                 sump.v = sump.v*ysq.v+calcei_plg(1).v
	                 sumq.v = sumq.v*ysq.v+calcei_qlg(1).v
	                 sump.v = sump.v*ysq.v+calcei_plg(2).v
	                 sumq.v = sumq.v*ysq.v+calcei_qlg(2).v
	                 sump.v = sump.v*ysq.v+calcei_plg(3).v
	                 sumq.v = sumq.v*ysq.v+calcei_qlg(3).v
	                 t0.v   = sumq.v*(x.v+x0.v)+frac.v
	                 ei.v   = (sump.v/t0.v)*xmx0.v
	                 
	                 if(jint==3) ei.v = exp(-x.v)*ei.v 
	             end if
	             
	       else if(all(msk8.m)) then
	       
	              frac.v = zero.v
	              frac.v= calcei_s(0).v/(calcei_r(0).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_s(1).v/(calcei_r(1).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_s(2).v/(calcei_r(2).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_s(3).v/(calcei_r(3).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_s(4).v/(calcei_r(4).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_s(5).v/(calcei_r(5).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_s(6).v/(calcei_r(6).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_s(7).v/(calcei_r(7).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_s(8).v/(calcei_r(8).v+ &
	                                          x.v+frac.v)
	              ei.v  = (calcei_r(9).v+frac.v)/x.v
	              
	              if(jint==3) ei.v = ei.v*exp(x.v)
	              
	        else if(all(msk9.m)) then
	        
	              frac.v = zero.v
	              frac.v= calcei_q1(0).v/(calcei_p1(0).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_q1(1).v/(calcei_p1(1).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_q1(2).v/(calcei_p1(2).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_q1(3).v/(calcei_p1(3).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_q1(4).v/(calcei_p1(4).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_q1(5).v/(calcei_p1(5).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_q1(6).v/(calcei_p1(6).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_q1(7).v/(calcei_p1(7).v+ &
	                                          x.v+frac.v)
	              frac.v= calcei_q1(8).v/(calcei_p1(8).v+ &
	                                          x.v+frac.v)
	              ei.v  = (calcei_p1(9).v+frac.v)/x.v
	              
	              if(jint/=3) ei.v = ei.v*exp(x.v)
	                                     
	        else
	             
	               msk10.m = (xmax.v<=x.v)
	               if(all(msk10.m).and.jint==3) then
	                   ei.v = xinf.v
	               else
	                   y.v    = one.v/x.v
	                   frac.v = zero.v
	                   frac.v= calcei_q2(0).v/(calcei_p2(0).v+ &
	                                               x.v+frac.v)
	                   frac.v= calcei_q2(1).v/(calcei_p2(1).v+ &
	                                          x.v+frac.v)
	                   frac.v= calcei_q2(2).v/(calcei_p2(2).v+ &
	                                          x.v+frac.v)
	                   frac.v= calcei_q2(3).v/(calcei_p2(3).v+ &
	                                          x.v+frac.v)
	                   frac.v= calcei_q2(4).v/(calcei_p2(4).v+ &
	                                          x.v+frac.v)
	                   frac.v= calcei_q2(5).v/(calcei_p2(5).v+ &
	                                          x.v+frac.v)
	                   frac.v= calcei_q2(6).v/(calcei_p2(6).v+ &
	                                          x.v+frac.v)
	                   frac.v= calcei_q2(7).v/(calcei_p2(7).v+ &
	                                          x.v+frac.v)
	                   frac.v= calcei_q2(8).v/(calcei_p2(8).v+ &
	                                          x.v+frac.v)   
	                   frac.v= calcei_p2(9).v+frac.v
	                   ei.v  = frac.v*y.v*y.v+y.v
	                   
	                   if(jint/=3) then
	                      msk10.m = (x.v<=(xmax.v-two4.v))
	                      if(all(msk10.m)) then
	                          ei.v = ei.v*exp(x.v)
	                      else
	                          t0.v = x.v-frty.v
	                          ei.v = ei.v*t0.v*exp40.v
	                      end if
	                   end if
	               end if
	           end if
	           val.v = ei.v
#endif
       end subroutine calcei_ymm4r8
       
       
#if 0
   *
    !*****************************************************************************80
!
!! CALCI0 computes various I0 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the first kind
!    and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards, Chalk
!    River (Atomic Energy of Canada Limited) Report AECL-4928,
!    October, 1974.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
!    the argument must be less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = I0(x);
!    2, RESULT = exp(-x) * I0(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, I0(x);
!    2, exp(-x) * I0(x);      
*/
#endif


        subroutine calci0_ymm4r8(arg,val,jint)
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: calci0_ymm4r8
              !dir$ attributes forceinline :: calci0_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: calci0_ymm4r8  
               type(YMM4r8_t),   intent(in)   :: arg
               type(YMM4r8_t),   intent(out)  :: val
               integer(kind=i4), intent(in)   :: jint 
               !dir$ attributes align : 64 :: one
               !dir$ attributes align : 64 :: one5
               !dir$ attributes align : 64 :: exp40
               !dir$ attributes align : 64 :: frty
               !dir$ attributes align : 64 :: rec15
               !dir$ attributes align : 64 :: two25
               !dir$ attributes align : 64 :: xsmall
               !dir$ attributes align : 64 :: xinf
               !dir$ attributes align : 64 :: xmax
               type(YMM4r8_t),   parameter    :: one   = YMM4r8_t(1.0e+0_dp)
               type(YMM4r8_t),   parameter    :: one5  = YMM4r8_t(15.0e+0_dp)
               type(YMM4r8_t),   parameter    :: exp40 = YMM4r8_t(2.353852668370199854e+17_dp)
               type(YMM4r8_t),   parameter    :: frty  = YMM4r8_t(40.0e+0_dp)
               type(YMM4r8_t),   parameter    :: rec15 = YMM4r8_t(6.6666666666666666666e-2_dp)
               type(YMM4r8_t),   parameter    :: two25 = YMM4r8_t(225.0e+0_dp)
               type(YMM4r8_t),   parameter    :: xsmall= YMM4r8_t(5.55e-17_dp)
               type(YMM4r8_t),   parameter    :: xinf  = YMM4r8_t(1.79e+308_dp)
               type(YMM4r8_t),   parameter    :: xmax  = YMM4r8_t(713.986e+0_dp)
               !dir$ attributes align : 64 :: a
               !dir$ attributes align : 64 :: b
               !dir$ attributes align : 64 :: sump
               !dir$ attributes align : 64 :: sumq
               !dir$ attributes align : 64 :: x
               !dir$ attributes align : 64 :: xx
               !dir$ attributes align : 64 :: t0
               !dir$ attributes align : 64 :: t1
               type(YMM4r8_t),   automatic    :: a,b
               type(YMM4r8_t),   automatic    :: sump,sumq
               type(YMM4r8_t),   automatic    :: x,xx
               type(YMM4r8_t),   automatic    :: t0,t1
               type(Mask4_t),    automatic    :: msk1,msk2
               type(Mask4_t),    automatic    :: msk3,msk4
#if (GMS_EXPLICIT_VECTORIZE) == 1
               integer(kind=i4) :: j
#endif   
#if (GMS_EXPLICIT_VECTORIZE) == 1
           
               !dir$ loop_count(4)
               !dir$ vector aligned
               !dir$ vector vectorlength(8)
               !dir$ vector always
               do j=0,3
                   x.v(j)    = abs(arg.v(j))
                   msk1.m(j) = (x.v(j)<xsmall.v(j))
                   msk2.m(j) = (x.v(j)<one5.v(j))
                   msk3.m(j) = (one5.v(j)<=x.v(j))
                   if(all(msk1.m(j))) then
                       val.v(j) = one.v(j)
                   else if(all(msk2.m(j))) then
                       xx.v(j)   = x.v(j)*x.v(j)
                       sump.v(j) = calci0_p(0).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(1).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(2).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(3).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(4).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(5).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(6).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(7).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(8).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(9).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(10).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(11).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(12).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(13).v(j)
                       sump.v(j) = sump.v(j)*xx.v(j)+calci0_p(14).v(j)
                       xx.v(j)   = xx.v(j)-two25.v(j)
                       sumq.v(j) = (((( &
                             xx.v(j)+calcei0_q(0).v(j)) &
                           * xx.v(j)+calcei0_q(1).v(j)) &
                           * xx.v(j)+calcei0_q(2).v(j)) &
                           * xx.v(j)+calcei0_q(3).v(j)) &
                           * xx.v(j)+calcei0_q(4).v(j)
                       val.v(j)  = sump.v(j)/sumq.v(j)
                       if(jint==2) val.v(j) = val.v(j)*exp(-x.v(j))
                   else if(all(msk3.m(j))) then
                           msk4.m(j) = (xmax.v(j)<=x.v(j))
                           if(jint==1.and.all(msk4.m(j))) then
                               val.v = xinf.v
                           else
                               xx.v(j)    = one.v(j)/(x.v(j)-rec15.v(j))
                               sump.v(j)  = ((((((   &
                                      calci0_pp(0).v(j)      &
                                   * xx.v(j)+calci0_pp(1).v(j))   &
                                   * xx.v(j)+calci0_pp(2).v(j))   &
                                   * xx.v(j)+calci0_pp(3).v(j))   &
                                   * xx.v(j)+calci0_pp(4).v(j))   &
                                   * xx.v(j)+calci0_pp(5).v(j))   &
                                   * xx.v(j)+calci0_pp(6).v(j))   &
                                   * xx.v(j)+calci0_pp(7).v(j)
                                sumq.v(j)  = ((((((    &
                                     xx.v(j)+calci0_qq(0).v(j)) &
                                   * xx.v(j)+calci0_qq(1).v(j)) &
                                   * xx.v(j)+calci0_qq(2).v(j)) &
                                   * xx.v(j)+calci0_qq(3).v(j)) &
                                   * xx.v(j)+calci0_qq(4).v(j)) &
                                   * xx.v(j)+calci0_qq(5).v(j)) &
                                   * xx.v(j)+calci0_qq(6).v(j)
                                val.v(j) = sump.v(j)/sumq.v(j)
                                if(jint==2) val.v(j) = (val.v(j)-calci0_pp(0).v(j)/sqrt(x.v(j)))
                   else
                        msk4.m(j) = (x.v(j)<=(xmamx.v(j)-one5.v(j)))
                        if(all(msk4.m(j))) then
                            a.v(j) = exp(x.v(j))
                            b.v(j) = one.v(j)
                        else
                            a.v(j) = exp(x.v(j)-frty.v(j))
                            b.v(j) = exp40.v(j)
                        end if
                        t0.v(j)  = calci0_pp(1).v(j)*a.v(j)
                        t1.v(j)  = sqrt(x.v(j))
                        val.v(j) = ((val.v(j)*a.v(j)-t0.v(j))/t1.v(j))*b.v(j)
                    end if
                  end if
               end if
           end do
#else
               x.v    = abs(arg.v)
               msk1.m = (x.v<xsmall.v)
               msk2.m = (x.v<one5.v)
               msk3.m = (one5.v<=x.v)
               if(all(msk1.m)) then
                  val.v = one.v
               else if(all(msk2.m)) then
                  xx.v   = x.v*x.v
                  sump.v = calci0_p(0).v
                  sump.v = sump.v*xx.v+calci0_p(1).v
                  sump.v = sump.v*xx.v+calci0_p(2).v
                  sump.v = sump.v*xx.v+calci0_p(3).v
                  sump.v = sump.v*xx.v+calci0_p(4).v
                  sump.v = sump.v*xx.v+calci0_p(5).v
                  sump.v = sump.v*xx.v+calci0_p(6).v
                  sump.v = sump.v*xx.v+calci0_p(7).v
                  sump.v = sump.v*xx.v+calci0_p(8).v
                  sump.v = sump.v*xx.v+calci0_p(9).v
                  sump.v = sump.v*xx.v+calci0_p(10).v
                  sump.v = sump.v*xx.v+calci0_p(11).v
                  sump.v = sump.v*xx.v+calci0_p(12).v
                  sump.v = sump.v*xx.v+calci0_p(13).v
                  sump.v = sump.v*xx.v+calci0_p(14).v
                  xx.v   = xx.v-two25.v
                  sumq.v = (((( &
                             xx.v+calcei0_q(0).v) &
                           * xx.v+calcei0_q(1).v) &
                           * xx.v+calcei0_q(2).v) &
                           * xx.v+calcei0_q(3).v) &
                           * xx.v+calcei0_q(4)
                  val.v  = sump.v/sumq.v
                  if(jint==2) val.v = val.v*exp(-x.v)
             else if(all(msk3.m)) then
                     msk4.m = (xmax.v<=x.v)
                     if(jint==1.and.all(msk4.m)) then
                         val.v = xinf.v
                     else
                         xx.v    = one.v/(x.v-rec15.v)
                         sump.v  = ((((((   &
                                      calci0_pp(0).v      &
                                   * xx.v+calci0_pp(1).v)   &
                                   * xx.v+calci0_pp(2).v)   &
                                   * xx.v+calci0_pp(3).v)   &
                                   * xx.v+calci0_pp(4).v)   &
                                   * xx.v+calci0_pp(5).v)   &
                                   * xx.v+calci0_pp(6).v)   &
                                   * xx.v+calci0_pp(7).v
                         sumq.v  = ((((((    &
                                     xx.v+calci0_qq(0).v) &
                                   * xx.v+calci0_qq(1).v) &
                                   * xx.v+calci0_qq(2).v) &
                                   * xx.v+calci0_qq(3).v) &
                                   * xx.v+calci0_qq(4).v) &
                                   * xx.v+calci0_qq(5).v) &
                                   * xx.v+calci0_qq(6).v
                        val.v    = sump.v/sumq.v
                        if(jint==2) val.v = (val.v-calci0_pp(0).v/sqrt(x.v))
                    else
                        msk4.m = (x.v<=(xmamx.v-one5.v))
                        if(all(msk4.m)) then
                            a.v = exp(x.v)
                            b.v = one.v
                        else
                            a.v = exp(x.v-frty.v)
                            b.v = exp40.v
                        end if
                        t0.v  = calci0_pp(1).v*a.v
                        t1.v  = sqrt(x.v)
                        val.v = ((val.v*a.v-t0.v)/t1.v)*b.v
                    end if
                  end if
               end if
#endif
        end subroutine calci0_ymm4r8
        
        
#if 0
  !*****************************************************************************80
!
!! CALCI1 computes various I1 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functioons of the first kind
!    and order one, I1(X) and EXP(-ABS(X))*I1(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards, Chalk
!    River (Atomic Energy of Canada Limited) Report AECL-4928,
!    October, 1974.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
!    the argument must be less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = I1(x);
!    2, RESULT = exp(-x) * I1(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, I1(x);
!    2, exp(-x) * I1(x);   
#endif


        subroutine calci1_ymm4r8(arg,val,jint)
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: calci1_ymm4r8
              !dir$ attributes forceinline :: calci1_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: calci1_ymm4r8  
              type(YMM4r8_t),   intent(in)   :: arg
              type(YMM4r8_t),   intent(out)  :: val
              integer(kind=i4), intent(in)   :: jint 
              !dir$ attributes align : 64 :: one
              !dir$ attributes align : 64 :: one5
              !dir$ attributes align : 64 :: exp40
              !dir$ attributes align : 64 :: frty
              !dir$ attributes align : 64 :: rec15
              !dir$ attributes align : 64 :: two25
              !dir$ attributes align : 64 :: xsmall
              !dir$ attributes align : 64 :: xinf
              !dir$ attributes align : 64 :: xmax
              !dir$ attributes align : 64 :: half
              !dir$ attributes align : 64 :: zero
              !dir$ attributes align : 64 :: pbar
              type(YMM4r8_t),   parameter    :: one   = YMM4r8_t(1.0e+0_dp)
              type(YMM4r8_t),   parameter    :: one5  = YMM4r8_t(15.0e+0_dp)
              type(YMM4r8_t),   parameter    :: exp40 = YMM4r8_t(2.353852668370199854e+17_dp)
              type(YMM4r8_t),   parameter    :: frty  = YMM4r8_t(40.0e+0_dp)
              type(YMM4r8_t),   parameter    :: rec15 = YMM4r8_t(6.6666666666666666666e-2_dp)
              type(YMM4r8_t),   parameter    :: two25 = YMM4r8_t(225.0e+0_dp)
              type(YMM4r8_t),   parameter    :: xsmall= YMM4r8_t(5.55e-17_dp)
              type(YMM4r8_t),   parameter    :: xinf  = YMM4r8_t(1.79e+308_dp)
              type(YMM4r8_t),   parameter    :: xmax  = YMM4r8_t(713.986e+0_dp)
              type(YMM4r8_t),   parameter    :: half  = YMM4r8_t(0.5e+00_dp)
              type(YMM4r8_t),   parameter    :: zero  = YMM4r8_t(0.0e+00_dp)
              type(YMM4r8_t),   parameter    :: pbar  = YMM4r8_t(3.98437500e-01_dp)
              !dir$ attributes align : 64 :: sump
              !dir$ attributes align : 64 :: sumq
              !dir$ attributes align : 64 :: x
              !dir$ attributes align : 64 :: a
              !dir$ attributes align : 64 :: b
              !dir$ attributes align : 64 :: t0
              !dir$ attributes align : 64 :: xx 
              type(YMM4r8_t),   automatic    :: sump,sumq
              type(YMM4r8_t),   automatic    :: x,a
              type(YMM4r8_t),   automatic    :: b,t0
              type(YMM4r8_t),   automatic    :: xx
              type(Mask4_t),    automatic    :: msk1,msk2
              type(Mask4_t),    automatic    :: msk3,msk4
              type(Mask4_t),    automatic    :: msk5
#if (GMS_EXPLICIT_VECTORIZE) == 1
               integer(kind=i4) :: j
#endif   
#if (GMS_EXPLICIT_VECTORIZE) == 1
           
               !dir$ loop_count(4)
               !dir$ vector aligned
               !dir$ vector vectorlength(8)
               !dir$ vector always
               do j=0,3
                     x.v(j)    = abs(arg.v(j))
                     msk1.m(j) = (x.v(j)<small.v(j))
                     msk2.m(j) = (x.v(j)<one5.v(j))
                     msk3.m(j) = (xmax.v(j)<x.v(j))
                     if(all(msk1.m(j))) then
                         val.v(j) = half.v(j)*x.v(j)
                     else if(all(msk2.m)) then
                         xx.v(j)  = x.v(j)*x.v(j)
                         sump.v(j)= calci1.p(0).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(1).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(2).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(3).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(4).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(6).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(7).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(8).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(9).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(10).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(11).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(12).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(13).v(j)
                         sump.v(j)= sump.v(j)*xx.v(j)+calci1.p(14).v(j)
                         xx.v(j)  = xx.v(j)-two25.v(j)
                         sumq.v(j)= (((((  &
                            xx.v(j)+calci1.q(0).v(j))  &
                          * xx.v(j)+calci1.q(1).v(j))  &
                          * xx.v(j)+calci1.q(2).v(j))  &
                          * xx.v(j)+calci1.q(3).v(j))  &
                          * xx.v(j)+calci1.q(4).v(j)
                         val.v(j) = (sump.v(j)/sumq.v(j))*x.v(j)
                         if(jint==2) val.v(j) = val.v(j)*exp(-x.v(j))
                    else if(jint==1.and.all(msk3.m(j))) then
                         val.v(j) = xinf.v(j)        
                    else
                         xx.v(j)   = one.v(j)/x.v(j)-rec15.v(j)
                         sump.v(j) = ((((((   &
                                      calci1.pp(0).v(j)     &
                               * xx.v(j)+calci1.pp(1).v(j))  &
                               * xx.v(j)+calci1.pp(2).v(j))  &
                               * xx.v(j)+calci1.pp(3).v(j))  &
                               * xx.v(j)+calci1.pp(4).v(j))  &
                               * xx.v(j)+calci1.pp(5).v(j))  &
                               * xx.v(j)+calci1.pp(6).v(j))  &
                               * xx.v(j)+calci1.pp(7).v(j)
                         sumq.v(j) = (((((    &
                                 xx.v(j)+calci1.qq(0).v(j))  &
                               * xx.v(j)+calci1.qq(1).v(j))  &
                               * xx.v(j)+calci1.qq(2).v(j))  &
                               * xx.v(j)+calci1.qq(3).v(j))  &
                               * xx.v(j)+calci1.qq(4).v(j))  &
                               * xx.v(j)+calci1.qq(5).v(j)
                         val.v(j)  = sump.v(j)/sumq.v(j)
                         msk4.m(j) = (xmax.v(j)-one.v(j)<x.v(j))
                         if(jint/=1) then 
                             val.v(j) = val.v(j)+pbar.v(j)/sqrt(x.v(j))
                         else
                         
                             if(all(msk4.m(j))) then
                                a.v(j) = exp(x.v(j)-frty40.v(j))
                                b.v(j) = exp40.v(j)
                             else
                                a.v(j) = exp(x.v(j))
                                b.v(j) = one.v(j)
                         end if
                         t0.v(j)   = val.v(j)*a.v(j)+pbar.v(j)*a.v(j)
                         val.v  = (t0.v(j)/sqrt(x.v(j)))*b.v(j)
                     end if
                 end if
                 msk5.m(j) = (arg.v(j)<zero.v(j))
                 if(all(msk5.m(j))) val.v(j) = -val.v(j)
            end do
#else              
              x.v    = abs(arg.v)
              msk1.m = (x.v<small.v)
              msk2.m = (x.v<one5.v)
              msk3.m = (xmax.v<x.v)
              if(all(msk1.m)) then
                  val.v = half.v*x.v
              else if(all(msk2.m)) then
                  xx.v  = x.v*x.v
                  sump.v= calci1.p(0).v
                  sump.v= sump.v*xx.v+calci1.p(1).v
                  sump.v= sump.v*xx.v+calci1.p(2).v
                  sump.v= sump.v*xx.v+calci1.p(3).v
                  sump.v= sump.v*xx.v+calci1.p(4).v
                  sump.v= sump.v*xx.v+calci1.p(6).v
                  sump.v= sump.v*xx.v+calci1.p(7).v
                  sump.v= sump.v*xx.v+calci1.p(8).v
                  sump.v= sump.v*xx.v+calci1.p(9).v
                  sump.v= sump.v*xx.v+calci1.p(10).v
                  sump.v= sump.v*xx.v+calci1.p(11).v
                  sump.v= sump.v*xx.v+calci1.p(12).v
                  sump.v= sump.v*xx.v+calci1.p(13).v
                  sump.v= sump.v*xx.v+calci1.p(14).v
                  xx.v  = xx.v-two25.v
                  sumq.v= (((((  &
                            xx.v+calci1.q(0).v)  &
                          * xx.v+calci1.q(1).v)  &
                          * xx.v+calci1.q(2).v)  &
                          * xx.v+calci1.q(3).v)  &
                          * xx.v+calci1.q(4).v
                  val.v = (sump.v/sumq.v)*x.v
                  if(jint==2) val.v = val.v*exp(-x.v)
              else if(jint==1.and.all(msk3.m)) then
                      val.v = xinf.v        
              else
                      xx.v   = one.v/x.v-rec15.v
                      sump.v = ((((((   &
                                      calci1.pp(0).v     &
                               * xx.v+calci1.pp(1).v)  &
                               * xx.v+calci1.pp(2).v)  &
                               * xx.v+calci1.pp(3).v)  &
                               * xx.v+calci1.pp(4).v)  &
                               * xx.v+calci1.pp(5).v)  &
                               * xx.v+calci1.pp(6).v)  &
                               * xx.v+calci1.pp(7).v
                      sumq.v = (((((    &
                                 xx.v+calci1.qq(0).v)  &
                               * xx.v+calci1.qq(1).v)  &
                               * xx.v+calci1.qq(2).v)  &
                               * xx.v+calci1.qq(3).v)  &
                               * xx.v+calci1.qq(4).v)  &
                               * xx.v+calci1.qq(5).v
                      val.v  = sump.v/sumq.v
                      msk4.m = (xmax.v-one.v<x.v)
                      if(jint/=1) then 
                         val.v = val.v+pbar.v/sqrt(x.v)
                      else
                         
                         if(all(msk4.m)) then
                            a.v = exp(x.v-frty40.v)
                            b.v = exp40.v
                         else
                            a.v = exp(x.v)
                            b.v = one.v
                         end if
                      t0.v   = val.v*a.v+pbar.v*a.v
                      val.v  = (t0.v/sqrt(x.v))*b.v
                  end if
              end if
              msk5.m = (arg.v<zero.v)
              if(all(msk5.m)) val.v = -val.v
#endif
        end subroutine calci1_ymm4r8
        
        
#if 0
/*
*****************************************************************************80
!
!! CALCK0 computes various K0 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order zero, K0(X) and EXP(X)*K0(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  0 < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K0(x);
!    2, RESULT = exp(x) * K0(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K0(x);
!    2, exp(x) * K0(x);
*/
#endif


      
         subroutine calck0_ymm4r8(val,arg,jint)
              
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: calck0_ymm4r8
              !dir$ attributes forceinline :: calck0_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: calck0_ymm4r8  
              type(YMM4r8_t),   intent(in)   :: arg
              type(YMM4r8_t),   intent(out)  :: val
              integer(kind=i4), intent(in)   :: jint  
              !dir$ attributes align : 64 :: zero
              !dir$ attributes align : 64 :: one
              !dir$ attributes align : 64 :: xsmall
              !dir$ attributes align : 64 :: xinf
              !dir$ attributes align : 64 :: xmax
              !dir$ attributes align : 64 :: sumf
              !dir$ attributes align : 64 :: sumg
              !dir$ attributes align : 64 :: sump
              !dir$ attributes align : 64 :: sumq
              !dir$ attributes align : 64 :: temp
              !dir$ attributes align : 64 :: x
              !dir$ attributes align : 64 :: xx
              !dir$ attributes align : 64 :: t0
              !dir$ attributes align : 64 :: t1
              !dir$ attributes align : 64 :: t2
              type(YMM4r8_t),   parameter    :: zero   = YMM4r8_t(0.0_dp)
              type(YMM4r8_t),   parameter    :: one    = YMM4r8_t(1.0_dp)
              type(YMM4r8_t),   parameter    :: xsmall = YMM4r8_t(1.11e-16_dp)
              type(YMM4r8_t),   parameter    :: xinf   = YMM4r8_t(1.79e+308_dp)
              type(YMM4r8_t),   parameter    :: xmax   = YMM4r8_t(705.342e+00_dp)
              type(YMM4r8_t),   automatic    :: sumf,sumg
              type(YMM4r8_t),   automatic    :: sump,sumq
              type(YMM4r8_t),   automatic    :: temp,x
              type(YMM4r8_t),   automatic    :: xx,t0,t1,t2
              type(Mask4_t),    automatic    :: msk1,msk2
              type(Mask4_t),    automatic    :: msk3,msk4
#if (GMS_EXPLICIT_VECTORIZE) == 1
               integer(kind=i4) :: j
#endif   
#if (GMS_EXPLICIT_VECTORIZE) == 1
           
               !dir$ loop_count(4)
               !dir$ vector aligned
               !dir$ vector vectorlength(8)
               !dir$ vector always
               do j=0,3
                   x.v(j)    = arg.v(j) 
                   msk1.m(j)  = (zero.v(j)<x.v(j))
                   msk4.m(j) = (xmax.v(j)<x.v(j))
                   if(all(msk1.m(j))) then
                       msk2.m(j) = (x.v(j)<=one.v(j))
                       if(all(msk2.m(j))) then
                           temp.v(j) = log(x.v(j))
                           msk3.m(j) = (x.v(j)<=xsmall.v(j))
                           if(all(msk3.m(j))) then
                              val.v(j) = calck0_p(5).v(j)/calck0_q(1).v(j)- &
                                         temp.v(j)
                           else
                      
                               xx.v(j)   = x.v(j)*x.v(j)
                               sump.v(j) = ((((  &
                                        calck0_p(0).v(j)   &
                                 * xx.v(j)+calck0_p(1).v(j))  &
                                 * xx.v(j)+calck0_p(2).v(j))  &
                                 * xx.v(j)+calck0_p(3).v(j))  &
                                 * xx.v(j)+calck0_p(4).v(j))  &
                                 * xx.v(j)+calck0_p(5).v(j)
                               sumq.v(j) = (xx.v(j)+calck0_q(0).v(j)) * &
                                   xx.v(j)+calck0_q(1).v(j)
                               sumf.v(j) = ((  &
                                         calck0_f(0).v(j)) &
                                  * xx.v(j)+calck0_f(1).v(j)) &
                                  * xx.v(j)+calck0_f(2).v(j)) &
                                  * xx.v(j)+calck0_f(3).v(j)
                               sumg.v(j) = ((xx.v(j)+calck0_g(0).v(j)) * &
                                    xx.v(j)+calck0_g(1).v(j)) * &
                                    xx.v(j)+calck0_g(2).v(j)
                               t0.v(j)   = sump.v(j)/sumq.v(j)
                               t1.v(j)   = xx.v(j)*sumf.v(j)
                               t2.v(j)   = temp.v(j)/sumg.v(j)-temp.v(j)
                               val.v(j)  = t0.v(j)-t1.v(j)*t2.v(j)
                               if(jint==2) val.v(j) = val.v(j)*exp(x.v(j))
                          end if
                    else if(jint==1.and.all(msk4.m(j))) then
                          val.v(j) = zero.v(j)
                    else
                          xx.v(j)  = one.v(j)/x.v(j)
                          t0.v(j)  = sqrt(x.v(j))
                          sump.v(j)= calck0_pp(0).v(j)
                          sump.v(j)= sump.v(j)*xx.v(j)+calck0_pp(1).v(j)
                          sump.v(j)= sump.v(j)*xx.v(j)+calck0_pp(2).v(j)
                          sump.v(j)= sump.v(j)*xx.v(j)+calck0_pp(3).v(j)
                          sump.v(j)= sump.v(j)*xx.v(j)+calck0_pp(4).v(j)
                          sump.v(j)= sump.v(j)*xx.v(j)+calck0_pp(5).v(j)
                          sump.v(j)= sump.v(j)*xx.v(j)+calck0_pp(6).v(j)
                          sump.v(j)= sump.v(j)*xx.v(j)+calck0_pp(7).v(j)
                          sump.v(j)= sump.v(j)*xx.v(j)+calck0_pp(8).v(j)
                          sump.v(j)= sump.v(j)*xx.v(j)+calck0_pp(9).v(j)
                          sumq.v(j)= xx.v(j)
                          sumq.v(j)= (sumq.v(j)+calck0_qq(1).v(j))*xx.v(j)
                          sumq.v(j)= (sumq.v(j)+calck0_qq(2).v(j))*xx.v(j)
                          sumq.v(j)= (sumq.v(j)+calck0_qq(3).v(j))*xx.v(j)
                          sumq.v(j)= (sumq.v(j)+calck0_qq(4).v(j))*xx.v(j)
                          sumq.v(j)= (sumq.v(j)+calck0_qq(5).v(j))*xx.v(j)
                          sumq.v(j)= (sumq.v(j)+calck0_qq(6).v(j))*xx.v(j)
                          sumq.v(j)= (sumq.v(j)+calck0_qq(7).v(j))*xx.v(j)
                          sumq.v(j)= (sumq.v(j)+calck0_qq(8).v(j))*xx.v(j)
                          sumq.v(j)= sumq.v(j)+calck0_qq(9).v(j)
                          val.v(j) = sump.v(j)/sumq.v(j)/t0.v(j)
                          if(jint==1) val.v(j) = val.v(j)*exp(x.v(j))
                 end if
              else
                 val.v(j) = xinf.v(j)
              end if
          end do
#else              
              x.v    = arg.v
              msk1.m = (zero.v<x.v)
              msk4.m = (xmax.v<x.v)
              if(all(msk1.m)) then
                  msk2.m = (x.v<=one.v)
                  if(all(msk2.m)) then
                      temp.v = log(x.v)
                      msk3.m = (x.v<=xsmall.v)
                      if(all(msk3.m)) then
                         val.v = calck0_p(5).v/calck0_q(1).v- &
                                 temp.v
                      else
                      
                         xx.v   = x.v*x.v
                         sump.v = ((((  &
                                        calck0_p(0).v   &
                                 * xx.v+calck0_p(1).v)  &
                                 * xx.v+calck0_p(2).v)  &
                                 * xx.v+calck0_p(3).v)  &
                                 * xx.v+calck0_p(4).v)  &
                                 * xx.v+calck0_p(5).v
                         sumq.v = (xx.v+calck0_q(0).v) * &
                                   xx.v+calck0_q(1).v
                         sumf.v = ((  &
                                         calck0_f(0).v) &
                                  * xx.v+calck0_f(1).v) &
                                  * xx.v+calck0_f(2).v) &
                                  * xx.v+calck0_f(3).v
                         sumg.v = ((xx.v+calck0_g(0).v) * &
                                    xx.v+calck0_g(1).v) * &
                                    xx.v+calck0_g(2).v
                         t0.v   = sump.v/sumq.v
                         t1.v   = xx.v*sumf.v
                         t2.v   = temp.v/sumg.v-temp.v
                         val.v  = t0.v-t1.v*t2.v
                         if(jint==2) val.v = val.v*exp(x.v)
                     end if
                 else if(jint==1.and.all(msk4.m)) then
                      val.v = zero.v
                 else
                      xx.v  = one.v/x.v
                      t0.v  = sqrt(x.v)
                      sump.v= calck0_pp(0).v
                      sump.v= sump.v*xx.v+calck0_pp(1).v
                      sump.v= sump.v*xx.v+calck0_pp(2).v
                      sump.v= sump.v*xx.v+calck0_pp(3).v
                      sump.v= sump.v*xx.v+calck0_pp(4).v
                      sump.v= sump.v*xx.v+calck0_pp(5).v
                      sump.v= sump.v*xx.v+calck0_pp(6).v
                      sump.v= sump.v*xx.v+calck0_pp(7).v
                      sump.v= sump.v*xx.v+calck0_pp(8).v
                      sump.v= sump.v*xx.v+calck0_pp(9).v
                      sumq.v= xx.v
                      sumq.v= (sumq.v+calck0_qq(1).v)*xx.v
                      sumq.v= (sumq.v+calck0_qq(2).v)*xx.v
                      sumq.v= (sumq.v+calck0_qq(3).v)*xx.v
                      sumq.v= (sumq.v+calck0_qq(4).v)*xx.v
                      sumq.v= (sumq.v+calck0_qq(5).v)*xx.v
                      sumq.v= (sumq.v+calck0_qq(6).v)*xx.v
                      sumq.v= (sumq.v+calck0_qq(7).v)*xx.v
                      sumq.v= (sumq.v+calck0_qq(8).v)*xx.v
                      sumq.v= sumq.v+calck0_qq(9).v
                      val.v = sump.v/sumq.v/t0.v
                      if(jint==1) val.v = val.v*exp(x.v)
                 end if
              else
                 val.v = xinf.v
              end if
#endif
         end subroutine calck0_ymm4r8
         
         
#if 0
  !*****************************************************************************80
!
!! CALCK1 computes various K1 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order one, K1(X) and EXP(X)*K1(X), for real arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  XLEAST < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K1(x);
!    2, RESULT = exp(x) * K1(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K1(x);
!    2, exp(x) * K1(x);  
#endif


       subroutine calck1_ymm4r8(val,arg,jint)
              
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: calck1_ymm4r8
              !dir$ attributes forceinline :: calck1_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: calck1_ymm4r8  
              type(YMM4r8_t),   intent(in)   :: arg
              type(YMM4r8_t),   intent(out)  :: val
              integer(kind=i4), intent(in)   :: jint  
              !dir$ attributes align : 64 :: zero
              !dir$ attributes align : 64 :: one
              !dir$ attributes align : 64 :: xsmall
              !dir$ attributes align : 64 :: xinf
              !dir$ attributes align : 64 :: xmax
              !dir$ attributes align : 64 :: sumf
              !dir$ attributes align : 64 :: sumg
              !dir$ attributes align : 64 :: sump
              !dir$ attributes align : 64 :: sumq
              !dir$ attributes align : 64 :: temp
              !dir$ attributes align : 64 :: x
              !dir$ attributes align : 64 :: xx
              !dir$ attributes align : 64 :: t0
              !dir$ attributes align : 64 :: t1
              !dir$ attributes align : 64 :: t2
              type(YMM4r8_t),   parameter    :: zero   = YMM4r8_t(0.0_dp)
              type(YMM4r8_t),   parameter    :: one    = YMM4r8_t(1.0_dp)
              type(YMM4r8_t),   parameter    :: xsmall = YMM4r8_t(1.11e-16_dp)
              type(YMM4r8_t),   parameter    :: xinf   = YMM4r8_t(1.79e+308_dp)
              type(YMM4r8_t),   parameter    :: xmax   = YMM4r8_t(705.342e+00_dp) 
              type(YMM4r8_t),   parameter    :: xleast = YMM4r8_t(2.23e-308_dp)
              type(YMM4r8_t),   automatic    :: sumf,sumg
              type(YMM4r8_t),   automatic    :: sump,sumq
              type(YMM4r8_t),   automatic    :: x,xx
              type(YMM4r8_t),   automatic    :: t0,t1,t2
              type(Mask4_t),    automatic    :: msk1,msk2
              type(Mask4_t),    automatic    :: msk3,msk4   
#if (GMS_EXPLICIT_VECTORIZE) == 1
               integer(kind=i4) :: j
#endif   
#if (GMS_EXPLICIT_VECTORIZE) == 1
           
               !dir$ loop_count(4)
               !dir$ vector aligned
               !dir$ vector vectorlength(8)
               !dir$ vector always
               do j=0,3
                    x.v(j)     = arg.v(j) 
                    msk1.m(j)  = (x.v(j)<=least.v(j))
                    msk2.m(j) = (x.v(j)<=one.v(j))
                    msk4.m(j) = (xmax.v(j)<x.v(j))
                    if(all(msk1.m(j))) then
                        val.v(j) = xinf.v(j)
                    else if(all(msk2.m(j))) then
                            msk3.m(j) = (x.v(j)<=xsmall.v(j))
                            if(all(msk3.m(j))) then
                                  val.v(j) = one.v(j)/x.v(j)
                            else
                                  xx.v(j)    = x.v(j)*x.v(j)
                                  sump.v(j)  = ((((    &
                                               calck1_p(0).v(j)       &
                                          *    xx.v(j)+calck1_p(1).v(j)) & 
                                          *    xx.v(j)+calck1_p(2).v(j)) &
                                          *    xx.v(j)+calck1_p(3).v(j)) &
                                          *    xx.v(j)+calck1_p(4).v(j)) &
                                          *    xx.v(j)+calck1_q(2).v(j)
                                  sumq.v(j)   = ((     &
                                               xx.v(j)+calck1_q(0).v(j)) &
                                          * xx.v(j)+calck1_q(1).v(j)) &
                                          * xx.v(j)+calck1_q(2).v(j))
                                  t1.v(j)     = sump.v(j)/sumq.v(j)  
                                  sumf.v(j)   = (((    &
                                              calck1_f(0).v(j)       &
                                          *    xx.v(j)+calck1_f(1).v(j)) &
                                          *    xx.v(j)+calck1_f(2).v(j)) &
                                          *    xx.v(j)+calck1_f(3).v(j)) &
                                          *    xx.v(j)+calck1_f(4).v(j)
                                  t2.v(j)     = xx.v(j)*log(x.v(j))
                                  sumg.v(j)   = (((    &
                                               xx.v(j)+calck1_g(0).v(j)) &
                                          *    xx.v(j)+calck1_g(1).v(j)) &
                                          *    xx.v(j)+calck1_g(2).v(j)
                                  t0.v(j)     = sumf.v(j)/sumg.v(j)
                                  val.v(j)    = (t2.v(j)*t0.v(j)+t1.v(j))/x.v(j)
                                  if(jint==2) val.v(j) = val.v(j)*exp(x.v(j))
                            end if
                 else if(all(msk4.m(j)).and.jint==1) then
                            val.v(j) = zero.v(j)
                 else
                            xx.v(j)  = one.v(j)/x.v(j)
                            sump.v(j)= calck1_pp(0).v(j)
                            sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(1).v(j)
                            sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(2).v(j)
                            sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(3).v(j)
                            sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(4).v(j)
                            sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(5).v(j)
                            sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(6).v(j)
                            sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(7).v(j)
                            sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(8).v(j)
                            sump.v(j)= sump.v(j)*xx.v+calck1_pp(9).v(j)
                            sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(10).v(j)
                            t0.v(j)  = sqrt(x.v(j))
                            sumq.v(j)= xx.v(j)
                            sumq.v(j)= (sumq.v(j)+calck1_qq(0).v(j))*xx.v(j)
                            sumq.v(j)= (sumq.v(j)+calck1_qq(1).v(j))*xx.v(j)
                            sumq.v(j)= (sumq.v(j)+calck1_qq(2).v(j))*xx.v(j)
                            sumq.v(j)= (sumq.v(j)+calck1_qq(3).v(j))*xx.v(j)
                            sumq.v(j)= (sumq.v(j)+calck1_qq(4).v(j))*xx.v(j)
                            sumq.v(j)= (sumq.v(j)+calck1_qq(5).v(j))*xx.v(j)
                            sumq.v(j)= (sumq.v(j)+calck1_qq(6).v(j))*xx.v(j)
                            sumq.v(j)= (sumq.v(j)+calck1_qq(7).v(j))*xx.v(j)
                            sumq.v(j)= sumq.v(j)+calck1_qq(8).v(j)
                            val.v = sump.v(j)/sumq.v(j)/t0.v(j)
                            if(jint==1) val.v(j) = val.v(j)*exp(-x.v(j))
                end if
            end do
#else         
              x.v    = arg.v
              msk1.m = (x.v<=least.v)
              msk2.m = (x.v<=one.v)
              msk4.m = (xmax.v<x.v)
              if(all(msk1.m)) then
                  val.v = xinf.v
              else if(all(msk2.m)) then
                  msk3.m = (x.v<=xsmall.v)
                  if(all(msk3.m)) then
                      val.v = one.v/x.v
                  else
                      xx.v    = x.v*x.v
                      sump.v  = ((((    &
                                    calck1_p(0).v       &
                               *    xx.v+calck1_p(1).v) & 
                               *    xx.v+calck1_p(2).v) &
                               *    xx.v+calck1_p(3).v) &
                               *    xx.v+calck1_p(4).v) &
                               *    xx.v+calck1_q(2).v
                      sumq.v   = ((     &
                                    xx.v+calck1_q(0).v) &
                                  * xx.v+calck1_q(1).v) &
                                  * xx.v+calck1_q(2).v)
                      t1.v     = sump.v/sumq.v  
                      sumf.v   = (((    &
                                    calck1_f(0).v       &
                               *    xx.v+calck1_f(1).v) &
                               *    xx.v+calck1_f(2).v) &
                               *    xx.v+calck1_f(3).v) &
                               *    xx.v+calck1_f(4).v
                      t2.v     = xx.v*log(x.v)
                      sumg.v   = (((    &
                                    xx.v+calck1_g(0).v) &
                               *    xx.v+calck1_g(1).v) &
                               *    xx.v+calck1_g(2).v
                      t0.v     = sumf.v/sumg.v
                      val.v    = (t2.v*t0.v+t1.v)/x.v
                      if(jint==2) val.v = val.v*exp(x.v)
                  end if
             else if(all(msk4.m).and.jint==1) then
                      val.v = zero.v
             else
                      xx.v  = one.v/x.v
                      sump.v= calck1_pp(0).v
                      sump.v= sump.v*xx.v+calck1_pp(1).v
                      sump.v= sump.v*xx.v+calck1_pp(2).v
                      sump.v= sump.v*xx.v+calck1_pp(3).v
                      sump.v= sump.v*xx.v+calck1_pp(4).v
                      sump.v= sump.v*xx.v+calck1_pp(5).v
                      sump.v= sump.v*xx.v+calck1_pp(6).v
                      sump.v= sump.v*xx.v+calck1_pp(7).v
                      sump.v= sump.v*xx.v+calck1_pp(8).v
                      sump.v= sump.v*xx.v+calck1_pp(9).v
                      sump.v= sump.v*xx.v+calck1_pp(10).v
                      t0.v  = sqrt(x.v)
                      sumq.v= xx.v
                      sumq.v= (sumq.v+calck1_qq(0).v)*xx.v
                      sumq.v= (sumq.v+calck1_qq(1).v)*xx.v
                      sumq.v= (sumq.v+calck1_qq(2).v)*xx.v
                      sumq.v= (sumq.v+calck1_qq(3).v)*xx.v
                      sumq.v= (sumq.v+calck1_qq(4).v)*xx.v
                      sumq.v= (sumq.v+calck1_qq(5).v)*xx.v
                      sumq.v= (sumq.v+calck1_qq(6).v)*xx.v
                      sumq.v= (sumq.v+calck1_qq(7).v)*xx.v
                      sumq.v= sumq.v+calck1_qq(8).v
                      val.v = sump.v/sumq.v/t0.v
                      if(jint==1) val.v = val.v*exp(-x.v)
              end if
#endif
        end subroutine calck1_ymm4r8
        
        
#if 0
    /*
    !*****************************************************************************80
!
!! CALJY0 computes various J0 and Y0 Bessel functions.
!
!  Discussion:
!
!    This routine computes zero-order Bessel functions of the first and
!    second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
!    for Y0, and |X| <= XMAX for J0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 0, ARG
!    must satisfy
!     -XMAX < ARG < XMAX;
!    If JINT = 1, then ARG must satisfy
!      0 < ARG < XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = J0(x);
!    1, RESULT = Y0(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    0, J0(x);
!    1, Y0(x);  
#endif


        subroutine caljy0_ymm4r8(arg,jint,val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: caljy0_ymm4r8
              !dir$ attributes forceinline :: caljy0_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: caljy0_ymm4r8  
              use mod_vectypes, only : ZMM8i8_t
              type(YMM4r8_t),   intent(in)   :: arg
              type(YMM4r8_t),   intent(out)  :: val
              integer(kind=i4), intent(in)   :: jint 
              !dir$ attributes align : 64 :: zero
              !dir$ attributes align : 64 :: one
              !dir$ attributes align : 64 :: three
              !dir$ attributes align : 64 :: four
              !dir$ attributes align : 64 :: eight
              !dir$ attributes align : 64 :: five5
              !dir$ attributes align : 64 :: sixty4
              !dir$ attributes align : 64 :: oneov8
              !dir$ attributes align : 64 :: p17
              !dir$ attributes align : 64 :: two56
              !dir$ attributes align : 64 :: cons
              !dir$ attributes align : 64 :: pi2
              !dir$ attributes align : 64 :: twopi
              !dir$ attributes align : 64 :: twopi1
              !dir$ attributes align : 64 :: twopi2
              !dir$ attributes align : 64 :: xmax
              !dir$ attributes align : 64 :: xsmall
              !dir$ attributes align : 64 :: xinf
              !dir$ attributes align : 64 :: xj0
              !dir$ attributes align : 64 :: xj1
              !dir$ attributes align : 64 :: xy0
              !dir$ attributes align : 64 :: xy1
              !dir$ attributes align : 64 :: xy2
              !dir$ attributes align : 64 :: xj01
              !dir$ attributes align : 64 :: xj02
              !dir$ attributes align : 64 :: xj11
              !dir$ attributes align : 64 :: xj12
              !dir$ attributes align : 64 :: xy01
              !dir$ attributes align : 64 :: xy02
              !dir$ attributes align : 64 :: xy11
              !dir$ attributes align : 64 :: xy12
              !dir$ attributes align : 64 :: xy21
              !dir$ attributes align : 64 :: xy22
              type(YMM4r8_t),   parameter :: zero  = YMM4r8_t(0.0e+0_dp);
              type(YMM4r8_t),   parameter :: one   = YMM4r8_t(1.0e+0_dp);
              type(YMM4r8_t),   parameter :: three = YMM4r8_t(3.0e+0_dp);
              type(YMM4r8_t),   parameter :: four  = YMM4r8_t(4.0e+0_dp);
              type(YMM4r8_t),   parameter :: eight = YMM4r8_t(8.0e+0_dp);
              type(YMM4r8_t),   parameter :: five5 = YMM4r8_t(5.5e+0_dp); 
              type(YMM4r8_t),   parameter :: sixty4= YMM4r8_t(64.0e+0_dp);
              type(YMM4r8_t),   parameter :: oneov8= YMM4r8_t(0.125e+0_dp); 
              type(YMM4r8_t),   parameter :: p17   = YMM4r8_t(1.716e-1_dp);
              type(YMM4r8_t),   parameter :: two56 = YMM4r8_t(256.0e+0_dp);
              type(YMM4r8_t),   parameter :: cons  = YMM4r8_t(-1.1593151565841244881e-1_dp);
              type(YMM4r8_t),   parameter :: pi2   = YMM4r8_t(6.3661977236758134308e-1_dp);
              type(YMM4r8_t),   parameter :: twopi = YMM4r8_t(6.2831853071795864769e+0_dp);
              type(YMM4r8_t),   parameter :: twopi1= YMM4r8_t(6.28125e+0_dp);
              type(YMM4r8_t),   parameter :: twopi2= YMM4r8_t(1.9353071795864769253e-3_dp);
              type(YMM4r8_t),   parameter :: xmax  = YMM4r8_t(1.07e+09_dp);
              type(YMM4r8_t),   parameter :: xsmall= YMM4r8_t(9.31e-10_dp);
              type(YMM4r8_t),   parameter :: xinf  = YMM4r8_t(1.7e+38_dp);
              type(YMM4r8_t),   parameter :: xj0   = YMM4r8_t(2.4048255576957727686e+0_dp);
              type(YMM4r8_t),   parameter :: xj1   = YMM4r8_t(5.5200781102863106496e+0_dp);
              type(YMM4r8_t),   parameter :: xy0   = YMM4r8_t(8.9357696627916752158e-1_dp);
              type(YMM4r8_t),   parameter :: xy1   = YMM4r8_t(3.9576784193148578684e+0_dp);
              type(YMM4r8_t),   parameter :: xy2   = YMM4r8_t(7.0860510603017726976e+0_dp);
              type(YMM4r8_t),   parameter :: xj01  = YMM4r8_t(616.0e+0_dp);
              type(YMM4r8_t),   parameter :: xj02  = YMM4r8_t(-1.4244423042272313784e-3_dp);
              type(YMM4r8_t),   parameter :: xj11  = YMM4r8_t(1413.0e+0_dp);
              type(YMM4r8_t),   parameter :: xj12  = YMM4r8_t(5.4686028631064959660e-4_dp);
              type(YMM4r8_t),   parameter :: xy01  = YMM4r8_t(228.0e+0_dp);
              type(YMM4r8_t),   parameter :: xy02  = YMM4r8_t(2.9519662791675215849e-3_dp);
              type(YMM4r8_t),   parameter :: xy11  = YMM4r8_t(1013.0e+0_dp);
              type(YMM4r8_t),   parameter :: xy12  = YMM4r8_t(6.4716931485786837568e-4_dp);
              type(YMM4r8_t),   parameter :: xy21  = YMM4r8_t(1814.0e+0_dp);
              type(YMM4r8_t),   parameter :: xy22  = YMM4r8_t(1.1356030177269762362e-4_dp);  
              !dir$ attributes align : 64 :: ax
              !dir$ attributes align : 64 :: down
              !dir$ attributes align : 64 :: prod
              !dir$ attributes align : 64 :: resj
              !dir$ attributes align : 64 :: r0
              !dir$ attributes align : 64 :: r1
              !dir$ attributes align : 64 :: up
              !dir$ attributes align : 64 :: w
              !dir$ attributes align : 64 :: wsq
              !dir$ attributes align : 64 :: xden
              !dir$ attributes align : 64 :: xy
              !dir$ attributes align : 64 :: z
              !dir$ attributes align : 64 :: zsq
              !dir$ attributes align : 64 :: t0
              !dir$ attributes align : 64 :: t1
              !dir$ attributes align : 64 :: pi2ax
              !dir$ attributes align : 64 :: t2
              !dir$ attributes align : 64 :: t3
              !dir$ attributes align : 64 :: ti
              type(YMM4r8_t),   automatic :: ax,down
              type(YMM4r8_t),   automatic :: prod,resj
              type(YMM4r8_t),   automatic :: r0,r1
              type(YMM4r8_t),   automatic :: up,w
              type(YMM4r8_t),   automatic :: wsq,xden
              type(YMM4r8_t),   automatic :: xy,z,zsq
              type(YMM4r8_t),   automatic :: t0,t1
              type(YMM4r8_t),   automatic :: pi2ax,t2
              type(YMM4r8_t),   automatic :: t3
              type(ZMM8i8_t),   automatic :: ti
              type(Mask4_t),    automatic :: m0,m1
#if (GMS_EXPLICIT_VECTORIZE) == 1
               integer(kind=i4) :: j
#endif   
              ax.v   = abs(arg.v)
              m0.m   = (arg.v<=zero.v)
              pi2ax.v= pi2.v/ax.v
              m1.m   = (xmax.v<ax.v)
              if(all(m0.m).and.jint==1) then
                  val.v = -xinf.v
                  return
              else if(all(m1.m)) then
                  val.v = zero.v
                  return
              end if
              m0.m = (eight.v<ax.v)
              m1.m = (ax.v<=small.v)
              if(all(m0.m)) goto 800
              if(all(m1.m)) then
                 if(jint==0) then
                     val.v = one.v
                 else
                     val.v = pi2.v*log(ax.v)+cons.v
                     return
                 end if
              end if
#if (GMS_EXPLICIT_VECTORIZE) == 1
           
               !dir$ loop_count(4)
               !dir$ vector aligned
               !dir$ vector vectorlength(8)
               !dir$ vector always
               do j=0,3
                    m0.m(j)  = (ax.v(j)<=four.v(j))
                    zsq.v(j) = ax.v(j)*ax.v(j)
                    if(all(m0.m(j))) then
                        xnum.v(j) = caljy0_pj0(4).v(j)*zsq.v(j)+ &
                                 caljy0_pj0(5).v(j)*zsq.v(j)+ &
                                 caljy0_pj0(6).v(j)
                        xden.v(j) = zsq.v(j)+caljy0_qj0(4).v(j)
                        xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_pj0(0).v(j)
                        xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qj0(0).v(j)
                        xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_pj0(1).v(j)
                        xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qj0(1).v(j)
                        xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_pj0(2).v(j)
                        xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qj0(2).v(j)
                        xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_pj0(3).v(j)
                        xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qj0(3).v(j)
                        t0.v(j)   = ax.v(j)-(xj01.v(j)/two56.v(j))
                        t1.v(j)   = ax.v(j)+xj0.v(j)
                        prod.v(j) = (t0.v(j)-xj02.v(j))*t1.v(j)
                    else
                        wsq.v(j)  = one.v(j)-(zsq.v(j)/sixty4.v(j))
                        xnum.v(j) = caljy0_pj1(6).v(j)*wsq.v(j)+ &
                                    caljy0_pji(7).v(j)
                        xden.v(j) = wsq.v(j)+caljy0_qj1(6).v(j)
                        xnum.v(j) = xnum.v(j)*wsq.v(j)+caljy0_pj1(0).v(j)
                        xden.v(j) = xden.v(j)*wsq.v(j)+caljy0_qj1(0).v(j)
                        xnum.v(j) = xnum.v(j)*wsq.v(j)+caljy0_pj1(1).v(j)
                        xden.v(j) = xden.v(j)*wsq.v(j)+caljy0_qj1(1).v(j)
                        xnum.v(j) = xnum.v(j)*wsq.v(j)+caljy0_pj1(2).v(j)
                        xden.v(j) = xden.v(j)*wsq.v(j)+caljy0_qj1(2).v(j)
                        xnum.v(j) = xnum.v(j)*wsq.v(j)+caljy0_pj1(3).v(j)
                        xden.v(j) = xden.v(j)*wsq.v(j)+caljy0_qj1(3).v(j)
                        xnum.v(j) = xnum.v(j)*wsq.v(j)+caljy0_pj1(4).v(j)
                        xden.v(j) = xden.v(j)*wsq.v(j)+caljy0_qj1(4).v(j)
                        xnum.v(j) = xnum.v(j)*wsq.v(j)+caljy0_pj1(5).v(j)
                        xden.v(j) = xden.v(j)*wsq.v(j)+caljy0_qj1(5).v(j)
                        t0.v(j)   = ax.v(j)-(xj11.v(j)/two56.v(j))
                        t1.v(j)   = ax.v(j)+xj1.v(j)
                        prod.v(j) = (t0.v(j)*t1.v(j))-xj12.v(j)
                   end if
                   val.v(j)  = prod.v(j)*(xnum.v(j)/xden.v(j))
                   if(jint==0) return
!    /*
!                          Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
!                          !  where xn is a zero of Y0.
!                       */      
                   m0.m(j)  = (ax.v(j)<=three.v(j))       
                   m1.m(j)  = (ax.v(j)<=five.v(j))
                   if(all(m0.m(j))) then
                       up.v(j) = (ax.v(j)-xy01.v(j)/two56.v(j))-xy02.v(j)
                       xy.v(j) = xy0.v(j)
                   else if(all(m1.m(j))) then
                       up.v(j) = (ax.v(j)-xy11.v(j)/two56.v(j))-xy12.v(j)
                       xy.v(j) = xy1.v(j)
                   else
                       up.v(j) = (ax.v(j)-xy21.v(j)/two56.v(j))-xy22.v(j)
                       xy.v(j) = xy1.v(j)
                   end if
                   down.v(j)   = ax.v(j)*xy.v(j)
                   t0.v(j)     = abs(up.v(j))
                   t1.v(j)     = p17.v(j)*down.v(j)
                   m0.m(j)     = (t0.v(j)<t1.v(j))
                   if(all(m0.m(j))) then
                        w.v(j)   = up.v(j)/down.v(j)
                        wsq.v(j) = w.v(j)*w.v(j)
                        xnum.v(j)= wsq.v(j)+caljy0_qlg(0).v(j)
                        xnum.v(j)= caljy0_plg(0).v(j)
                        xden.v(j)= wsq.v(j)+caljy0_qlg(0).v(j)
                        xnum.v(j)= xnum.v(j)*wsq.v(j)+caljy0_plg(1).v(j)
                        xden.v(j)= xden.v(j)*wsq.v(j)+caljy0_qlg(1).v(j)
                        xnum.v(j)= xnum.v(j)*wsq.v(j)+caljy0_plg(2).v(j)
                        xden.v(j)= xden.v(j)*wsq.v(j)+caljy0_qlg(2).v(j)
                        xnum.v(j)= xnum.v(j)*wsq.v(j)+caljy0_plg(3).v(j)
                        xden.v(j)= xden.v(j)*wsq.v(j)+caljy0_qlg(3).v(j)
                        t0.v(j)  = xnum.v(j)/xden.v(j)
                        t1.v(j)  = pi.v(j)*val.v(j)
                        resj.v(j)= t1.v(j)*w.v(j)*t0.v(j)
                  else
                        t0.v(j)  = xnum.v(j)/xden.v(j)
                        t1.v(j)  = pi.v(j)*val.v(j)
                        resj.v(j)= t1.v(j)*log(t0.v(j))
                  end if
!               /*
!                           Now calculate Y0 for appropriate interval, preserving
!                           !  accuracy near the zero of Y0.
!                       */
                   m0.v(j) = (ax.v(j)<=three.v(j))
                   m1.v(j) = (ax.v(j)<=five5.v(j))
                   if(all(m0.m(j))) then
                      xnum.v(j) = caljy0_py0(5).v(j)*zsq.v(j)+ &
                                  caljy0_py0(0).v(j)
                      xden.v(j) = zsq.v(j)+caljy0_qy0(0).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py0(1).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qy0(1).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py0(2).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qy0(2).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py0(3).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qy0(3).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py0(4).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qy0(4).v(j)
                  else if(all(m1.m(j))) then
                      xnum.v(j) = caljy0_py1(6).v(j)*zsq.v(j)+ &
                                  caljy0_py1(0).v(j)
                      xden.v(j) = zsq.v(j)+caljy0_qy1(0).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py1(1).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qy1(1).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py1(2).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qy1(2).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py1(3).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qy1(3).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py1(4).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qy1(4).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py1(5).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy0_qy1(5).v(j)
                 else
                      xnum.v(j) = caljy0_py2(7).v(j)*zsq.v(j)+ &
                                  caljy0_py2(0).v(j)
                      xden.v(j) = zsq.v(j)+caljy0_qy2(0).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py2(1).v(j)
                      xden.v(j) = xnum.v(j)*zsq.v(j)+caljy0_qy2(1).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py2(2).v(j)
                      xden.v(j) = xnum.v(j)*zsq.v(j)+caljy0_qy2(2).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py2(3).v(j)
                      xden.v(j) = xnum.v(j)*zsq.v(j)+caljy0_qy2(3).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py2(4).v(j)
                      xden.v(j) = xnum.v(j)*zsq.v(j)+caljy0_qy2(4).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py2(5).v(j)
                      xden.v(j) = xnum.v(j)*zsq.v(j)+caljy0_qy2(5).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy0_py2(6).v(j)
                      xden.v(j) = xnum.v(j)*zsq.v(j)+caljy0_qy2(6).v(j)
                 end if
                 t0.v(j)  = xnum.v(j)/xden.v(j)
                 t1.v(j)  = up.v(j)*down.v(j)
                 val.v(j) = t0.v(j)*t1.v(j)+resj.v(j)
                 return
800              z.v(j)   = eight.v(j)/ax.v(j)
                 w.v(j)   = ax.v(j)/twopi.v(j)
                 t1.v(j)  = int(w.v(j),kind=i8)
                 w.v(j)   = real(t1.v(j),kind=dp)+oneov8.v(j)
                 t0.v(j)  = w.v(j)*twopi2.v(j)
                 t1.v(j)  = ax.v(j)-w.v(j)
                 w.v(j)   = t1.v(j)*twopi1.v(j)-t0.v(j)
                 zsq.v(j) = z.v(j)*z.v(j)
                 xnum.v(j)= caljy0_p0(4).v(j)*zsq.v(j)+ &
                      caljy0_p0(5).v(j)
                 xden.v(j)= zsq.v(j)+caljy0_q0(4).v(j)
                 up.v(j)  = caljy0_p1(4).v(j)*zsq.v(j)+ &
                      caljy0_p1(5).v(j)
                 down.v(j)= zsq.v(j)+caljy0_q1(4).v(j)
                 xnum.v(j)= xnum.v(j)*zsq.v(j)+caljy0_p0(0).v(j)
                 xden.v(j)= xden.v(j)*zsq.v(j)+caljy0_q0(0).v(j)
                 xnum.v(j)= xnum.v(j)*zsq.v(j)+caljy0_p0(1).v(j)
                 xden.v(j)= xden.v(j)*zsq.v(j)+caljy0_q0(1).v(j)
                 xnum.v(j)= xnum.v(j)*zsq.v(j)+caljy0_p0(2).v(j)
                 xden.v(j)= xden.v(j)*zsq.v(j)+caljy0_q0(2).v(j)
                 xnum.v(j)= xnum.v(j)*zsq.v(j)+caljy0_p0(3).v(j)
                 xden.v(j)= xden.v(j)*zsq.v(j)+caljy0_q0(3).v(j)
                 r0.v(j)  = xnum.v(j)/xden.v(j)
                 t1.v(j)  = cos(w.v(j))
                 r1.v(j)  = up.v(j)/down.v(j)
                 t0.v(j)  = sqrt(pi2ax.v(j))
                 t2.v(j)  = sin(w.v(j))
                 t3.v(j)  = z.v(j)*r1.v(j)
                 if(jint==0) then
                     val.v(j) = t0.v(j)*r0.v(j)*t1.v(j)- &
                          t3.v(j)*t2.v(j)
                 else
                     val.v(j) = t0.v(j)*r0.v(j)*t2.v(j)+ &
                             t3.v(j)*t1.v(j)
                 end if
           end do
#else
             
!               /*
!                            !  Calculate J0 for appropriate interval, preserving
!                            !  accuracy near the zero of J0.
!                        
! */                     
              m0.m  = (ax.v<=four.v)
              zsq.v = ax.v*ax.v
              if(all(m0.m)) then
                  xnum.v = caljy0_pj0(4).v*zsq.v+ &
                           caljy0_pj0(5).v*zsq.v+ &
                           caljy0_pj0(6).v
                  xden.v = zsq.v+caljy0_qj0(4).v
                  xnum.v = xnum.v*zsq.v+caljy0_pj0(0).v
                  xden   = xden.v*zsq.v+caljy0_qj0(0).v
                  xnum.v = xnum.v*zsq.v+caljy0_pj0(1).v
                  xden   = xden.v*zsq.v+caljy0_qj0(1).v
                  xnum.v = xnum.v*zsq.v+caljy0_pj0(2).v
                  xden   = xden.v*zsq.v+caljy0_qj0(2).v
                  xnum.v = xnum.v*zsq.v+caljy0_pj0(3).v
                  xden   = xden.v*zsq.v+caljy0_qj0(3).v
                  t0.v   = ax.v-(xj01.v/two56.v)
                  t1.v   = ax.v+xj0.v
                  prod.v = (t0.v-xj02.v)*t1.v
              else
                  wsq.v  = one.v-(zsq.v/sixty4.v)
                  xnum.v = caljy0_pj1(6).v*wsq.v+ &
                           caljy0_pji(7).v
                  xden.v = wsq.v+caljy0_qj1(6).v
                  xnum.v = xnum.v*wsq.v+caljy0_pj1(0).v
                  xden.v = xden.v*wsq.v+caljy0_qj1(0).v
                  xnum.v = xnum.v*wsq.v+caljy0_pj1(1).v
                  xden.v = xden.v*wsq.v+caljy0_qj1(1).v
                  xnum.v = xnum.v*wsq.v+caljy0_pj1(2).v
                  xden.v = xden.v*wsq.v+caljy0_qj1(2).v
                  xnum.v = xnum.v*wsq.v+caljy0_pj1(3).v
                  xden.v = xden.v*wsq.v+caljy0_qj1(3).v
                  xnum.v = xnum.v*wsq.v+caljy0_pj1(4).v
                  xden.v = xden.v*wsq.v+caljy0_qj1(4).v
                  xnum.v = xnum.v*wsq.v+caljy0_pj1(5).v
                  xden.v = xden.v*wsq.v+caljy0_qj1(5).v
                  t0.v   = ax.v-(xj11.v/two56.v)
                  t1.v   = ax.v+xj1.v
                  prod.v = (t0.v*t1.v)-xj12.v
              end if
              val.v  = prod.v*(xnum.v/xden.v)
              if(jint==0) return
!    /*
!                          Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
!                          !  where xn is a zero of Y0.
!                       */      
              m0.m  = (ax.v<=three.v)       
              m1.m  = (ax.v<=five.v)
              if(all(m0.m)) then
                  up.v = (ax.v-xy01.v/two56.v)-xy02.v
                  xy.v = xy0.v
              else if(all(m1.m)) then
                  up.v = (ax.v-xy11.v/two56.v)-xy12.v
                  xy.v = xy1.v
              else
                  up.v = (ax.v-xy21.v/two56.v)-xy22.v
                  xy.v = xy1.v
              end if
              down.v   = ax.v*xy.v
              t0.v     = abs(up.v)
              t1.v     = p17.v*down.v
              m0.m     = (t0.v<t1.v)
              if(all(m0.m)) then
                  w.v   = up.v/down.v
                  wsq.v = w.v*w.v
                  xnum.v= wsq.v+caljy0_qlg(0).v
                  xnum.v= caljy0_plg(0).v
                  xden.v= wsq.v+caljy0_qlg(0).v
                  xnum.v= xnum.v*wsq.v+caljy0_plg(1).v
                  xden.v= xden.v*wsq.v+caljy0_qlg(1).v
                  xnum.v= xnum.v*wsq.v+caljy0_plg(2).v
                  xden.v= xden.v*wsq.v+caljy0_qlg(2).v
                  xnum.v= xnum.v*wsq.v+caljy0_plg(3).v
                  xden.v= xden.v*wsq.v+caljy0_qlg(3).v
                  t0.v  = xnum.v/xden.v
                  t1.v  = pi.v*val.v
                  resj.v= t1.v*w.v*t0.v
              else
                  t0.v  = xnum.v/xden.v
                  t1.v  = pi.v*val.v
                  resj.v= t1.v*log(t0.v)
              end if
!               /*
!                           Now calculate Y0 for appropriate interval, preserving
!                           !  accuracy near the zero of Y0.
!                       */
              m0.v = (ax.v<=three.v)
              m1.v = (ax.v<=five5.v)
              if(all(m0.m)) then
                  xnum.v = caljy0_py0(5).v*zsq.v+ &
                           caljy0_py0(0).v
                  xden.v = zsq.v+caljy0_qy0(0).v
                  xnum.v = xnum.v*zsq.v+caljy0_py0(1).v
                  xden.v = xden.v*zsq.v+caljy0_qy0(1).v
                  xnum.v = xnum.v*zsq.v+caljy0_py0(2).v
                  xden.v = xden.v*zsq.v+caljy0_qy0(2).v
                  xnum.v = xnum.v*zsq.v+caljy0_py0(3).v
                  xden.v = xden.v*zsq.v+caljy0_qy0(3).v
                  xnum.v = xnum.v*zsq.v+caljy0_py0(4).v
                  xden.v = xden.v*zsq.v+caljy0_qy0(4).v
              else if(all(m1.m)) then
                  xnum.v = caljy0_py1(6).v*zsq.v+ &
                           caljy0_py1(0).v
                  xden.v = zsq.v+caljy0_qy1(0).v
                  xnum.v = xnum.v*zsq.v+caljy0_py1(1).v
                  xden.v = xden.v*zsq.v+caljy0_qy1(1).v
                  xnum.v = xnum.v*zsq.v+caljy0_py1(2).v
                  xden.v = xden.v*zsq.v+caljy0_qy1(2).v
                  xnum.v = xnum.v*zsq.v+caljy0_py1(3).v
                  xden.v = xden.v*zsq.v+caljy0_qy1(3).v
                  xnum.v = xnum.v*zsq.v+caljy0_py1(4).v
                  xden.v = xden.v*zsq.v+caljy0_qy1(4).v
                  xnum.v = xnum.v*zsq.v+caljy0_py1(5).v
                  xden.v = xden.v*zsq.v+caljy0_qy1(5).v
              else
                  xnum.v = caljy0_py2(7).v*zsq.v+ &
                           caljy0_py2(0).v
                  xden.v = zsq.v+caljy0_qy2(0).v
                  xnum.v = xnum.v*zsq.v+caljy0_py2(1).v
                  xden.v = xnum.v*zsq.v+caljy0_qy2(1).v
                  xnum.v = xnum.v*zsq.v+caljy0_py2(2).v
                  xden.v = xnum.v*zsq.v+caljy0_qy2(2).v
                  xnum.v = xnum.v*zsq.v+caljy0_py2(3).v
                  xden.v = xnum.v*zsq.v+caljy0_qy2(3).v
                  xnum.v = xnum.v*zsq.v+caljy0_py2(4).v
                  xden.v = xnum.v*zsq.v+caljy0_qy2(4).v
                  xnum.v = xnum.v*zsq.v+caljy0_py2(5).v
                  xden.v = xnum.v*zsq.v+caljy0_qy2(5).v
                  xnum.v = xnum.v*zsq.v+caljy0_py2(6).v
                  xden.v = xnum.v*zsq.v+caljy0_qy2(6).v
              end if
              t0.v  = xnum.v/xden.v
              t1.v  = up.v*down.v
              val.v = t0.v*t1.v+resj.v
              return
800           z.v   = eight.v/ax.v
              w.v   = ax.v/twopi.v
              t1.v  = int(w.v,kind=i8)
              w.v   = real(t1.v,kind=dp)+oneov8.v
              t0.v  = w.v*twopi2.v
              t1.v  = ax.v-w.v
              w.v   = t1.v*twopi1.v-t0.v
              zsq.v = z.v*z.v
              xnum.v= caljy0_p0(4).v*zsq.v+ &
                      caljy0_p0(5).v
              xden.v= zsq.v+caljy0_q0(4).v
              up.v  = caljy0_p1(4).v*zsq.v+ &
                      caljy0_p1(5).v
              down.v= zsq.v+caljy0_q1(4).v
              xnum.v= xnum.v*zsq.v+caljy0_p0(0).v
              xden.v= xden.v*zsq.v+caljy0_q0(0).v
              xnum.v= xnum.v*zsq.v+caljy0_p0(1).v
              xden.v= xden.v*zsq.v+caljy0_q0(1).v
              xnum.v= xnum.v*zsq.v+caljy0_p0(2).v
              xden.v= xden.v*zsq.v+caljy0_q0(2).v
              xnum.v= xnum.v*zsq.v+caljy0_p0(3).v
              xden.v= xden.v*zsq.v+caljy0_q0(3).v
              r0.v  = xnum.v/xden.v
              t1.v  = cos(w.v)
              r1.v  = up.v/down.v
              t0.v  = sqrt(pi2ax.v)
              t2.v  = sin(w.v)
              t3.v  = z.v*r1.v
              if(jint==0) then
                  val.v = t0.v*r0.v*t1.v- &
                          t3.v*t2.v
              else
                  val.v = t0.v*r0.v*t2.v+ &
                          t3.v*t1.v
              end if
#endif
        end subroutine caljy0_ymm4r8
        
        
#if 0
      /*
	             !*****************************************************************************80
!
!! CALJY1 computes various J1 and Y1 Bessel functions.
!
!  Discussion:
!
!    This routine computes first-order Bessel functions of the first and
!    second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
!    for Y1, and |X| <= XMAX for J1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 0, ARG
!    must satisfy
!     -XMAX < ARG < XMAX;
!    If JINT = 1, then ARG must satisfy
!      0 < ARG < XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = J1(x);
!    1, RESULT = Y1(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    0, J1(x);
!    1, Y1(x);  
#endif


          subroutine caljy1_ymm4r8(arg,jint,val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: caljy1_ymm4r8
              !dir$ attributes forceinline :: caljy1_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: caljy1_ymm4r8  
              use mod_vectypes, only : ZMM8i8_t
              type(YMM4r8_t),   intent(in)   :: arg
              type(YMM4r8_t),   intent(out)  :: val
              integer(kind=i4), intent(in)   :: jint 
              !dir$ attributes align : 64 :: eight
              !dir$ attributes align : 64 :: four
              !dir$ attributes align : 64 :: half
              !dir$ attributes align : 64 :: throv8
              !dir$ attributes align : 64 :: pi2
              !dir$ attributes align : 64 :: p17
              !dir$ attributes align : 64 :: twopi
              !dir$ attributes align : 64 :: zero
              !dir$ attributes align : 64 :: twopi1
              !dir$ attributes align : 64 :: twopi2
              !dir$ attributes align : 64 :: two56
              !dir$ attributes align : 64 :: rtpi2
              !dir$ attributes align : 64 :: xmax
              !dir$ attributes align : 64 :: xsmall
              !dir$ attributes align : 64 :: xinf
              !dir$ attributes align : 64 :: xj0
              !dir$ attributes align : 64 :: xj1
              !dir$ attributes align : 64 :: xy0
              !dir$ attributes align : 64 :: xy1
              !dir$ attributes align : 64 :: xj01
              !dir$ attributes align : 64 :: xj02
              !dir$ attributes align : 64 :: xj11
              !dir$ attributes align : 64 :: xj12
              !dir$ attributes align : 64 :: xy01
              !dir$ attributes align : 64 :: xy02
              !dir$ attributes align : 64 :: xy11
              !dir$ attributes align : 64 :: xy12
              !dir$ attributes align : 64 :: ax
              !dir$ attributes align : 64 :: down
              !dir$ attributes align : 64 :: prod
              !dir$ attributes align : 64 :: resj
              !dir$ attributes align : 64 :: r0
              !dir$ attributes align : 64 :: r1
              !dir$ attributes align : 64 :: up
              !dir$ attributes align : 64 :: w
              !dir$ attributes align : 64 :: wsq
              !dir$ attributes align : 64 :: xden
              !dir$ attributes align : 64 :: xnum
              !dir$ attributes align : 64 :: t0
              !dir$ attributes align : 64 :: t1
              !dir$ attributes align : 64 :: z
              !dir$ attributes align : 64 :: zsq
              !dir$ attributes align : 64 :: t2
              !dir$ attributes align : 64 :: t3
              type(YMM4r8_t),   parameter    :: eight = YMM4r8_t(8.0e+0_dp);
              type(YMM4r8_t),   parameter    :: four  = YMM4r8_t(4.0e+0_dp);
              type(YMM4r8_t),   parameter    :: half  = YMM4r8_t(0.5_dp);
              type(YMM4r8_t),   parameter    :: throv8= YMM4r8_t(0.375_dp);
              type(YMM4r8_t),   parameter    :: pi2   = YMM4r8_t(6.3661977236758134308e-1_dp);
              type(YMM4r8_t),   parameter    :: p17   = YMM4r8_t(1.716e-1_dp);
              type(YMM4r8_t),   parameter    :: twopi = YMM4r8_t(6.2831853071795864769e+0_dp);
              type(YMM4r8_t),   parameter    :: zero  = YMM4r8_t(0.0_dp);
              type(YMM4r8_t),   parameter    :: twopi1= YMM4r8_t(6.28125e+0_dp);
              type(YMM4r8_t),   parameter    :: twopi2= YMM4r8_t(1.9353071795864769253e-3_dp);
              type(YMM4r8_t),   parameter    :: two56 = YMM4r8_t(256.0e+0_dp);
              type(YMM4r8_t),   parameter    :: rtpi2 = YMM4r8_t(7.9788456080286535588e-1_dp);
              type(YMM4r8_t),   parameter    :: xmax  = YMM4r8_t(1.07e+9_dp);
              type(YMM4r8_t),   parameter    :: xsmall= YMM4r8_t(9.31e-10_dp);
              type(YMM4r8_t),   parameter    :: xinf  = YMM4r8_t(1.7e+38_dp);  
              type(YMM4r8_t),   parameter    :: xj0   = YMM4r8_t(3.8317059702075123156e+0_dp);
              type(YMM4r8_t),   parameter    :: xj1   = YMM4r8_t(7.0155866698156187535e+0_dp);
              type(YMM4r8_t),   parameter    :: xy0   = YMM4r8_t(2.1971413260310170351e+0_dp);
              type(YMM4r8_t),   parameter    :: xy1   = YMM4r8_t(5.4296810407941351328e+0_dp);
              type(YMM4r8_t),   parameter    :: xj01  = YMM4r8_t(981.0e+0_dp);
              type(YMM4r8_t),   parameter    :: xj02  = YMM4r8_t(-3.2527979248768438556e-4_dp);
              type(YMM4r8_t),   parameter    :: xj11  = YMM4r8_t(1796.0e+0_dp);
              type(YMM4r8_t),   parameter    :: xj12  = YMM4r8_t(-3.8330184381246462950e-5_dp);
              type(YMM4r8_t),   parameter    :: xy01  = YMM4r8_t(562.0e+0_dp);
              type(YMM4r8_t),   parameter    :: xy02  = YMM4r8_t(1.8288260310170351490e-3_dp);
              type(YMM4r8_t),   parameter    :: xy11  = YMM4r8_t(1390.0e+0_dp);
              type(YMM4r8_t),   parameter    :: xy12  = YMM4r8_t(-6.4592058648672279948e-6_dp);
              type(YMM4r8_t),   automatic    :: ax,down
              type(YMM4r8_t),   automatic    :: prod,resj
              type(YMM4r8_t),   automatic    :: r0,r1
              type(YMM4r8_t),   automatic    :: up,w
              type(YMM4r8_t),   automatic    :: wsq,xden
              type(YMM4r8_t),   automatic    :: xnum,t0
              type(YMM4r8_t),   automatic    :: t1,z
              type(YMM4r8_t),   automatic    :: zsq,t2
              type(YMM4r8_t),   automatic    :: t3
              type(Mask4_t),    automatic    :: m0,m1,m2,m3
#if (GMS_EXPLICIT_VECTORIZE) == 1
               integer(kind=i4) :: j
#endif   
              ax.v  = arg.v
              m3.m  = (xmax.v<ax.v)
              t0.v  = ax.v*xinf.v
              m0.m  = (arg.v<=zero.v)
              m1.m  = (arg.v<half.v)
              m2.m  = (t0.v<pi2.v)
              if(all(m0.m).or.                 &
                 (all(m1.m).and.all(m2.m))) then
                    val.v = -xinf.v
                    return
              else if(all(m3.m)) then
                    val.v = zero.v
                    return
              end if
              m0.m  = (eight.v<ax.v)
              m1.m  = (ax.v<=xsmall.v)
              if(all(m0.m)) then
                  goto 800
              else if(all(m1.m)) then
                   if(jint==0) then
                       val.v = arg.v*half.v
                       return
                   else
                       val.v = -pi2.v/ax.v 
                       return
                   end if
              end if 
              ! /*
              !                Calculate J1 for appropriate interval, preserving
              !                !  accuracy near the zero of J1.
              !           */
#if (GMS_EXPLICIT_VECTORIZE) == 1
           
               !dir$ loop_count(4)
               !dir$ vector aligned
               !dir$ vector vectorlength(8)
               !dir$ vector always
               do j=0,3
               
                   zsq.v(j) = ax.v(j)*ax.v(j)
                   m0.m(j)  = (ax.v(j)<=four.v(j))
                   if(all(m0.m(j))) then
                         xnum.v(j) = caljy1_pj0(6).v(j)*zsq.v(j)+ &
                                  caljy1_pj0(5).v(j)*zsq.v(j)+ &
                                  caljy1_pj0(4).v(j)
                         xden.v(j) = zsq.v(j)+caljy1_qj0(4).v(j)
                         xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_pj0(0).v(j)
                         xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qj0(0).v(j)
                         xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_pj0(1).v(j)
                         xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qj0(1).v(j)
                         xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_pj0(2).v(j)
                         xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qj0(2).v(j)
                         xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_pj0(3).v(j)
                         xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qj0(3).v(j)
                         t0.v(j)   = ax.v(j)-(xj01.v(j)/two56.v(j))
                         t1.v(j)   = ax.v(j)+xj0.v(j)
                         prod.v(j) = (t0.v(j)-xj02.v(j))*t1.v(j)
                   else
                         xnum.v(j) = caljy1_pj1(0).v(j)
                         xden.v(j) = (zsq.v(j)+caljy1_qj1(6).v(j))* &
                                     (zsq.v(j)+caljy1_qj1(0).v(j))
                         xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_pj1(1).v(j)
                         xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qj1(1).v(j)
                         xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_pj1(2).v(j)
                         xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qj1(2).v(j)
                         xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_pj1(3).v(j)
                         xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qj1(3).v(j)
                         xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_pj1(4).v(j)
                         xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qj1(4).v(j)
                         xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_pj1(5).v(j)
                         xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qj1(5).v(j)
                         t0.v(j)   = xnum.v(j)*(ax.v(j)-eight.v(j))
                         t1.v(j)   = (ax.v(j)+eight.v(j))+caljy1_pj1(6).v(j)
                         xnum.v(j) = t0.v(j)*t1.v(j)
                         t0.v(j)   = xnum.v(j)*(ax.v(j)-four.v(j))
                         t1.v(j)   = (ax.v(j)+four.v(j))+caljy1_pj1(7).v(j)
                         xnum.v(j) = t0.v(j)*t1.v(j)
                         t0.v(j)   = ax.v(j)-(xj11.v(j)/two56.v(j))-xj12.v(j)
                         t1.v(j)   = ax.v(j)+xj1.v(j)
                         prod.v(j) = arg.v(j)*t0.v(j)*t1.v(j)
                  endif
                  val.v(j) = prod.v(j)*(xnum.v(j)/xden.v(j))
                  if(jint==0) return
              ! /*
              !               Calculate Y1.  First find RESJ = pi/2 ln(x/xn) J1(x),
              !               !  where xn is a zero of Y1.
              !           */
                  m0.m(j) = (ax.v(j)<=four.v(j))
                  if(all(m0.m(j))) then
                       t0.v(j) = ax.v(j)-(xy01.v(j)/two56.v(j))
                       up.v(j) = t0.v(j)-xy02.v(j)
                       xy.v(j) = xy0.v(j)
                  else
                       t0.v(j) = ax.v(j)-(xy11.v(j)/two56.v(j))
                       up.v(j) = t0.v(j)-xy12.v(j)
                       xy.v(j) = xy01.v(j)
                  endif
                  down.v(j)   = ax.v(j)+xy.v(j)
                  t1.v(j)     = p17.v(j)*down.v(j))
                  m0.m(j)     = (abs(up.v(j))<t1.v(j))
                  if(all(m0.m(j))) then
                      w.v(j)   = up.v(j)/down.v(j)
                      wsq.v(j) = w.v*w.v(j)
                      xnum.v(j)= caljy1_plg(0).v(j)
                      xden.v(j)= wsq.v(j)+caljy1_qlg(0).v(j)
                      xnum.v(j)= xnum.v(j)*wsq.v(j)+caljy1_plg(1).v(j)
                      xden.v(j)= xden.v(j)*wsq.v(j)+caljy1_qlg(1).v(j)
                      xnum.v(j)= xnum.v(j)*wsq.v(j)+caljy1_plg(2).v(j)
                      xden.v(j)= xden.v(j)*wsq.v(j)+caljy1_qlg(2).v(j)
                      xnum.v(j)= xnum.v(j)*wsq.v(j)+caljy1_plg(3).v(j)
                      xden.v(j)= xden.v(j)*wsq.v(j)+caljy1_qlg(3).v(j)
                      t0.v(j)  = w.v(j)*(xnum.v(j)/xden.v(j))
                      t1.v(j)  = pi2.v(j)*val.v(j)
                      resj.v(j)= t0.v(j)*t1.v(j)
                  else
                      t0.v(j)  = log(ax.v(j)/xy.v(j))
                      resj.v(j)= pi12.v(j)*val.v(j)*t0.v(j)
                  end if
              !  /*
              !               Now calculate Y1 for appropriate interval, preserving
              !               !  accuracy near the zero of Y1.
              !           */   
                  m0.m(j) = (ax.v(j)<=four.v(j))
                  if(all(m0.m(j))) then
                      xnum.v(j) = caljy1_py0(6).v(j)*zsq.v(j)+ &
                                  caljy1_py0(0).v(j)         
                      xden.v(j) = zsq.v(j)+caljy1_qy0(0).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py0(1).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy0(1).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py0(2).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy0(2).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py0(3).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy0(3).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py0(4).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy0(4).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py0(5).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy0(5).v(j)  
                  else
                      xnum.v(j) = caljy1_py1(8).v(j)*zsq.v(j)+ &
                                  caljy1_py1(0).v(j)
                      xden.v(j) = zsq.v(j)+caljy1_qy1(0).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py1(1).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy1(1).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py1(2).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy1(2).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py1(3).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy1(3).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py1(4).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy1(4).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py1(5).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy1(5).v(j)
                      xnum.v(j) = xnum.v(j)*zsq.v(j)+caljy1_py1(6).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy1(6).v(j)
                      xnum.v(j) = xnum.v*zsq.v+caljy1_py1(7).v(j)
                      xden.v(j) = xden.v(j)*zsq.v(j)+caljy1_qy1(7).v(j)
              end if
              t0.v(j) = xnum.v(j)/xden.v(j)
              t1.v(j) = up.v(j)*(down.v(j)/ax.v(j))
              val.v(j)= t0.v(j)*t1.v(j)+resj.v(j)
              return
800           z.v(j)  = eight.v(j)/ax.v(j)
              t0.v(j) = ax.v(j)/twopi.v(j)
              ti.v(j) = int(t0.v(j),kind=i8)
              w.v(j)  = real(ti.v(j),kind=dp)*throv8.v(j)
              t0.v(j) = ax.v(j)-w.v(j)
              w.v(j)  = t0.v(j)*twopi1.v(j)-w.v(j)*twopi2.v(j)
              zsq.v(j)= z.v(j)*z.v(j)
              xnum.v(j)=caljy1_p0(5).v(j)
              xden.v(j)=zsq.v(j)+caljy1_q0(5).v
              up.v(j)  = caljy1_p1(5).v(j)
              xnum.v(j)= xnum.v(j)*zsq.v(j)+caljy1_p0(0).v(j)
              t0.v(j)  = rtpi2.v(j)/sqrt(ax.v(j))
              xden.v(j)= xden.v(j)*zsq.v(j)+caljy1_q0(0).v(j)
              xnum.v(j)= xnum.v(j)*zsq.v(j)+caljy1_p0(1).v(j)
              xden.v(j)= xden.v(j)*zsq.v(j)+caljy1_q0(1).v(j)
              xnum.v(j)= xnum.v(j)*zsq.v(j)+caljy1_p0(2).v(j)
              xden.v(j)= xden.v(j)*zsq.v(j)+caljy1_q0(2).v(j)
              xnum.v(j)= xnum.v(j)*zsq.v(j)+caljy1_p0(3).v(j)
              xden.v(j)= xden.v(j)*zsq.v(j)+caljy1_q0(3).v(j)
              xnum.v(j)= xnum.v(j)*zsq.v(j)+caljy1_p0(4).v(j)
              t1.v(j)  = sin(w.v(j))
              xden.v(j)= xden.v(j)*zsq.v(j)+caljy1_q0(4).v(j)
              r0.v(j)  = xnum.v(j)/xden.v(j)
              r1.v(j)  = up.v(j)/down.v(j)
              t2.v(j)  = cos(w.v(j))
              t3.v(j)  = z.v(j)*r1.v(j)
              if(jint==1) then
                   val.v = t0.v(j)*r0.v(j)*t2.v(j)- &
                           t3.v(j)*t1.v(j)
              else
                   val.v(j) = t0.v(j)*r0.v(j)*t1.v(j)- &
                           t3.v(j)*t2.v(j)
              end if
              m0.m(j)  = (arg.v(j)<zero.v(j))
              if(jint==0.and.all(m0.m(j))) val.v(j) = -val.v(j)
           end do
#else                 
              zsq.v = ax.v*ax.v
              m0.m  = (ax.v<=four.v)
              if(all(m0.m)) then
                   xnum.v = caljy1_pj0(6).v*zsq.v+ &
                                  caljy1_pj0(5).v*zsq.v+ &
                                  caljy1_pj0(4).v
                   xden.v = zsq.v+caljy1_qj0(4).v
                   xnum.v = xnum.v*zsq.v+caljy1_pj0(0).v
                   xden.v = xden.v*zsq.v+caljy1_qj0(0).v
                   xnum.v = xnum.v*zsq.v+caljy1_pj0(1).v
                   xden.v = xden.v*zsq.v+caljy1_qj0(1).v
                   xnum.v = xnum.v*zsq.v+caljy1_pj0(2).v
                   xden.v = xden.v*zsq.v+caljy1_qj0(2).v
                   xnum.v = xnum.v*zsq.v+caljy1_pj0(3).v
                   xden.v = xden.v*zsq.v+caljy1_qj0(3).v
                   t0.v   = ax.v-(xj01.v/two56.v)
                   t1.v   = ax.v+xj0.v
                   prod.v = (t0.v-xj02.v)*t1.v
              else
                   xnum.v = caljy1_pj1(0).v
                   xden.v = (zsq.v+caljy1_qj1(6).v)* &
                            (zsq.v+caljy1_qj1(0).v)
                   xnum.v = xnum.v*zsq.v+caljy1_pj1(1).v
                   xden.v = xden.v*zsq.v+caljy1_qj1(1).v
                   xnum.v = xnum.v*zsq.v+caljy1_pj1(2).v
                   xden.v = xden.v*zsq.v+caljy1_qj1(2).v
                   xnum.v = xnum.v*zsq.v+caljy1_pj1(3).v
                   xden.v = xden.v*zsq.v+caljy1_qj1(3).v
                   xnum.v = xnum.v*zsq.v+caljy1_pj1(4).v
                   xden.v = xden.v*zsq.v+caljy1_qj1(4).v
                   xnum.v = xnum.v*zsq.v+caljy1_pj1(5).v
                   xden.v = xden.v*zsq.v+caljy1_qj1(5).v
                   t0.v   = xnum.v*(ax.v-eight.v)
                   t1.v   = (ax.v+eight.v)+caljy1_pj1(6).v
                   xnum.v = t0.v*t1.v
                   t0.v   = xnum.v*(ax.v-four.v)
                   t1.v   = (ax.v+four.v)+caljy1_pj1(7).v
                   xnum.v = t0.v*t1.v
                   t0.v   = ax.v-(xj11.v/two56.v)-xj12.v
                   t1.v   = ax.v+xj1.v
                   prod.v = arg.v*t0.v*t1.v
              endif
              val.v = prod.v*(xnum.v/xden.v)
              if(jint==0) return
              ! /*
              !               Calculate Y1.  First find RESJ = pi/2 ln(x/xn) J1(x),
              !               !  where xn is a zero of Y1.
              !           */
              m0.m = (ax.v<=four.v)
              if(all(m0.m)) then
                  t0.v = ax.v-(xy01.v/two56.v)
                  up.v = t0.v-xy02.v
                  xy.v = xy0.v
              else
                  t0.v = ax.v-(xy11.v/two56.v)
                  up.v = t0.v-xy12.v
                  xy.v = xy01.v
              endif
              down.v   = ax.v+xy.v
              t1.v     = p17.v*down.v)
              m0.m     = (abs(up.v)<t1.v)
              if(all(m0.m)) then
                  w.v   = up.v/down.v
                  wsq.v = w.v*w.v
                  xnum.v= caljy1_plg(0).v
                  xden.v= wsq.v+caljy1_qlg(0).v
                  xnum.v= xnum.v*wsq.v+caljy1_plg(1).v
                  xden.v= xden.v*wsq.v+caljy1_qlg(1).v
                  xnum.v= xnum.v*wsq.v+caljy1_plg(2).v
                  xden.v= xden.v*wsq.v+caljy1_qlg(2).v
                  xnum.v= xnum.v*wsq.v+caljy1_plg(3).v
                  xden.v= xden.v*wsq.v+caljy1_qlg(3).v
                  t0.v  = w.v*(xnum.v/xden.v)
                  t1.v  = pi2.v*val.v
                  resj.v= t0.v*t1.v
              else
                  t0.v  = log(ax.v/xy.v)
                  resj.v= pi12.v*val.v*t0.v
              end if
              !  /*
              !               Now calculate Y1 for appropriate interval, preserving
              !               !  accuracy near the zero of Y1.
              !           */   
              m0.m = (ax.v<=four.v)
              if(all(m0.m)) then
                   xnum.v = caljy1_py0(6).v*zsq.v+ &
                            caljy1_py0(0).v         
                   xden.v = zsq.v+caljy1_qy0(0).v
                   xnum.v = xnum.v*zsq.v+caljy1_py0(1).v
                   xden.v = xden.v*zsq.v+caljy1_qy0(1).v
                   xnum.v = xnum.v*zsq.v+caljy1_py0(2).v
                   xden.v = xden.v*zsq.v+caljy1_qy0(2).v
                   xnum.v = xnum.v*zsq.v+caljy1_py0(3).v
                   xden.v = xden.v*zsq.v+caljy1_qy0(3).v
                   xnum.v = xnum.v*zsq.v+caljy1_py0(4).v
                   xden.v = xden.v*zsq.v+caljy1_qy0(4).v
                   xnum.v = xnum.v*zsq.v+caljy1_py0(5).v
                   xden.v = xden.v*zsq.v+caljy1_qy0(5).v  
              else
                   xnum.v = caljy1_py1(8).v*zsq.v+ &
                            caljy1_py1(0).v
                   xden.v = zsq.v+caljy1_qy1(0).v
                   xnum.v = xnum.v*zsq.v+caljy1_py1(1).v
                   xden.v = xden.v*zsq.v+caljy1_qy1(1).v
                   xnum.v = xnum.v*zsq.v+caljy1_py1(2).v
                   xden.v = xden.v*zsq.v+caljy1_qy1(2).v
                   xnum.v = xnum.v*zsq.v+caljy1_py1(3).v
                   xden.v = xden.v*zsq.v+caljy1_qy1(3).v
                   xnum.v = xnum.v*zsq.v+caljy1_py1(4).v
                   xden.v = xden.v*zsq.v+caljy1_qy1(4).v
                   xnum.v = xnum.v*zsq.v+caljy1_py1(5).v
                   xden.v = xden.v*zsq.v+caljy1_qy1(5).v
                   xnum.v = xnum.v*zsq.v+caljy1_py1(6).v
                   xden.v = xden.v*zsq.v+caljy1_qy1(6).v
                   xnum.v = xnum.v*zsq.v+caljy1_py1(7).v
                   xden.v = xden.v*zsq.v+caljy1_qy1(7).v
              end if
              t0.v = xnum.v/xden.v
              t1.v = up.v*(down.v/ax.v)
              val.v= t0.v*t1.v+resj.v
              return
800           z.v  = eight.v/ax.v
              t0.v = ax.v/twopi.v
              ti.v = int(t0.v,kind=i8)
              w.v  = real(ti.v,kind=dp)*throv8.v
              t0.v = ax.v-w.v
              w.v  = t0.v*twopi1.v-w.v*twopi2.v
              zsq.v= z.v*z.v
              xnum.v=caljy1_p0(5).v
              xden.v=zsq.v+caljy1_q0(5).v
              up.v  = caljy1_p1(5).v
              xnum.v= xnum.v*zsq.v+caljy1_p0(0).v
              t0.v  = rtpi2.v/sqrt(ax.v)
              xden.v= xden.v*zsq.v+caljy1_q0(0).v
              xnum.v= xnum.v*zsq.v+caljy1_p0(1).v
              xden.v= xden.v*zsq.v+caljy1_q0(1).v
              xnum.v= xnum.v*zsq.v+caljy1_p0(2).v
              xden.v= xden.v*zsq.v+caljy1_q0(2).v
              xnum.v= xnum.v*zsq.v+caljy1_p0(3).v
              xden.v= xden.v*zsq.v+caljy1_q0(3).v
              xnum.v= xnum.v*zsq.v+caljy1_p0(4).v
              t1.v  = sin(w.v)
              xden.v= xden.v*zsq.v+caljy1_q0(4).v
              r0.v  = xnum.v/xden.v
              r1.v  = up.v/down.v
              t2.v  = cos(w.v)
              t3.v  = z.v*r1.v
              if(jint==1) then
                   val.v = t0.v*r0.v*t2.v- &
                           t3.v*t1.v
              else
                   val.v = t0.v*r0.v*t1.v- &
                           t3.v*t2.v
              end if
              m0.m  = (arg.v<zero.v)
              if(jint==0.and.all(m0.m)) val.v = -val.v
#endif                  
       end subroutine caljy1_ymm4r8
       
#if 0
   *****************************************************************************80
!
!! CALCK0 computes various K0 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order zero, K0(X) and EXP(X)*K0(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  0 < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K0(x);
!    2, RESULT = exp(x) * K0(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K0(x);
!    2, exp(x) * K0(x);
#endif


          subroutine calck0_ymm4r8(arg,jint,val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: calck0_ymm4r8
              !dir$ attributes forceinline :: calck0_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: calck0_ymm4r8  
              type(YMM4r8_t),   intent(in)   :: arg
              type(YMM4r8_t),   intent(out)  :: val
              integer(kind=i4), intent(in)   :: jint 
              !dir$ attributes align : 64 :: one
              !dir$ attributes align : 64 :: zero
              !dir$ attributes align : 64 :: xsmall
              !dir$ attributes align : 64 :: xinf
              !dir$ attributes align : 64 :: xmax
              !dir$ attributes align : 64 :: sumf
              !dir$ attributes align : 64 :: sumg
              !dir$ attributes align : 64 :: sump
              !dir$ attributes align : 64 :: sumq
              !dir$ attributes align : 64 :: x
              !dir$ attributes align : 64 :: xx
              !dir$ attributes align : 64 :: temp
              !dir$ attributes align : 64 :: t0
              !dir$ attributes align : 64 :: t1
              !dir$ attributes align : 64 :: t2
              type(YMM4r8_t),   parameter    :: one    = YMM4r8_t(1.0e+0_dp);
              type(YMM4r8_t),   parameter    :: zero   = YMM4r8_t(0.0e+0_dp);
              type(YMM4r8_t),   parameter    :: xsmall = YMM4r8_t(1.11e-16_dp);
              type(YMM4r8_t),   parameter    :: xinf   = YMM4r8_t(1.79e+308_dp);
              type(YMM4r8_t),   parameter    :: xmax   = YMM4r8_t(705.342e+0_dp);
              type(YMM4r8_t),   automatic    :: sumf,sumg
              type(YMM4r8_t),   automatic    :: sump,sumq
              type(YMM4r8_t),   automatic    :: x,xx
              type(YMM4r8_t),   automatic    :: temp,t0
              type(YMM4r8_t),   automatic    :: t1,t2
              type(Mask4_t),    automatic    :: m0,m1,m2
#if (GMS_EXPLICIT_VECTORIZE) == 1
               integer(kind=i4) :: j
#endif   
#if (GMS_EXPLICIT_VECTORIZE) == 1
           
              !dir$ loop_count(84
              !dir$ vector aligned
              !dir$ vector vectorlength(8)
              !dir$ vector always
              do j=0,3
                   x.v(j)  = arg.v(j)
                   m0.m(j) = (zero.v(j)<x.v(j))
                   if(all(m0.m(j))) then
                      m1.m(j) = (x.v(j)<one.v(j))
                      m0.m(j) = (xmax.v(j)<x.v(j))
                      if(all(m1.m(j))) then
                         temp.v(j) = log(x.v(j))
                         m2.m(j)   = (x.v(j)<small.v(j))
                         if(all(m2.m(j))) then
                            val.v(j) = (calck0_p(5).v(j)/calck0_p(1).v(j))- &
                                       temp.v(j)
                         else
                            xx.v(j)  = x.v(j)*x.v(j)
                            sump.v(j)= ((((                 &
                                   calck0_p(0).v(j)         &
                                * xx.v(j)+calck0_p(1).v(j)) &
                                * xx.v(j)+calck0_p(2).v(j)) &
                                * xx.v(j)+calck0_p(3).v(j)) &
                                * xx.v(j)+calck0_p(4).v(j)) &
                                * xx.v(j)+calck0_p(5).v(j)
                            sumq.v(j)= (xx.v(j)+calck0_q(0).v(j))* &
                                  xx.v(j)+calck0_q(1).v(j)
                            sumf.v(j)= ((    &
                                    calck0_f(0).v(j)      &
                                *  xx.v(j)+calck0_f(1).v(j)) &
                                *  xx.v(j)+calck0_f(2).v(j)) &
                                *  xx.v(j)+calck0_f(3).v(j)
                            sumg.v(j)= ((xx.v(j)+calck0_g(0).v(j)) &
                                *   xx.v(j)+calck0_g(1).v(j)) &
                                *   xx.v(j)+calck0_g(2).v(j)
                            t0.v(j)  = sump.v(j)/sumq.v(j)
                            t1.v(j)  = xx.v(j)*sumf.v(j)
                            t2.v(j)  = temp.v(j)/sumg.v(j)-temp.v(j)
                            val.v(j) = t0.v(j)-t1.v(j)*t2.v(j)
                            if(jint==2) val.v(j) = exp(val.v(j))
                       end if
                  else if(jint==1.and.all(m0.m(j))) then
                        val.v(j) = zero.v(j)
                  else
                        xx.v(j)   = one.v(j)/x.v(j)
                        sumq.v(j) = xx.v(j)
                        sump.v(j) = calck0_pp(0).v(j)
                        sump.v(j) = sump.v(j)*xx.v(j)+calck0_pp(1).v(j)
                        sumq.v(j) = (sumq.v(j)+calck0_qq(0).v(j))*xx.v(j)
                        sump.v(j) = sump.v(j)*xx.v(j)+calck0_pp(2).v(j)
                        sumq.v(j) = (sumq.v(j)+calck0_qq(1).v(j))*xx.v(j)
                        sump.v(j) = sump.v(j)*xx.v(j)+calck0_pp(3).v(j)
                        sumq.v(j) = (sumq.v(j)+calck0_qq(2).v(j))*xx.v(j)
                        sump.v(j) = sump.v(j)*xx.v(j)+calck0_pp(4).v(j)
                        sumq.v(j) = (sumq.v(j)+calck0_qq(3).v(j))*xx.v(j)
                        sump.v(j) = sump.v(j)*xx.v(j)+calck0_pp(5).v(j)
                        sumq.v(j) = (sumq.v(j)+calck0_qq(4).v(j))*xx.v(j)
                        sump.v(j) = sump.v(j)*xx.v(j)+calck0_pp(6).v(j)
                        sumq.v(j) = (sumq.v(j)+calck0_qq(5).v(j))*xx.v(j)
                        sump.v(j) = sump.v(j)*xx.v(j)+calck0_pp(7).v(j)
                        sumq.v(j) = (sumq.v(j)+calck0_qq(6).v(j))*xx.v(j)
                        sump.v(j) = sump.v(j)*xx.v(j)+calck0_pp(8).v(j)
                        sumq.v(j) = (sumq.v(j)+calck0_qq(7).v(j))*xx.v(j)
                        sump.v(j) = sump.v(j)*xx.v(j)+calck0_pp(9).v(j)
                        sumq.v(j) = (sumq.v(j)+calck0_qq(8).v(j))*xx.v(j)
                        sumq.v(j) = sumq.v(j)+calck0_qq(9).v(j)
                        t0.v(j)   = sqrt(x.v(j))
                        t1.v(j)   = sump.v(j)/sumq.v(j)
                        val.v(j)  = t1.v(j)/t0.v(j)
                        if(jint==1) val.v(j) = val.v(j)*exp(-x.v(j))
                 end if
             else
                 val.v(j) = xinf.v(j)
              end if
         end do      
#else         
              x.v  = arg.v
              m0.m = (zero.v<x.v)
              if(all(m0.m)) then
                  m1.m = (x.v<one.v)
                  m0.m = (xmax.v<x.v)
                  if(all(m1.m)) then
                      temp.v = log(x.v)
                      m2.m   = (x.v<small.v)
                      if(all(m2.m)) then
                          val.v = (calck0_p(5).v/calck0_p(1).v)- &
                                   temp.v
                      else
                          xx.v  = x.v*x.v
                          sump.v= ((((                &
                                   calck0_p(0).v      &
                                * xx.v+calck0_p(1).v) &
                                * xx.v+calck0_p(2).v) &
                                * xx.v+calck0_p(3).v) &
                                * xx.v+calck0_p(4).v) &
                                * xx.v+calck0_p(5).v
                          sumq.v= (xx.v+calck0_q(0).v)* &
                                   xx.v+calck0_q(1).v
                          sumf.v= ((    &
                                    calck0_f(0).v      &
                                *  xx.v+calck0_f(1).v) &
                                *  xx.v+calck0_f(2).v) &
                                *  xx.v+calck0_f(3).v
                          sumg.v= ((xx.v+calck0_g(0).v) &
                                *   xx.v+calck0_g(1).v) &
                                *   xx.v+calck0_g(2).v
                          t0.v  = sump.v/sumq.v
                          t1.v  = xx.v*sumf.v
                          t2.v  = temp.v/sumg.v-temp.v
                          val.v = t0.v-t1.v*t2.v
                          if(jint==2) val.v = exp(val.v)
                     end if
                 else if(jint==1.and.all(m0.m)) then
                     val.v = zero.v
                 else
                     xx.v   = one.v/x.v
                     sumq.v = xx.v
                     sump.v = calck0_pp(0).v
                     sump.v = sump.v*xx.v+calck0_pp(1).v
                     sumq.v = (sumq.v+calck0_qq(0).v)*xx.v
                     sump.v = sump.v*xx.v+calck0_pp(2).v
                     sumq.v = (sumq.v+calck0_qq(1).v)*xx.v
                     sump.v = sump.v*xx.v+calck0_pp(3).v
                     sumq.v = (sumq.v+calck0_qq(2).v)*xx.v
                     sump.v = sump.v*xx.v+calck0_pp(4).v
                     sumq.v = (sumq.v+calck0_qq(3).v)*xx.v
                     sump.v = sump.v*xx.v+calck0_pp(5).v
                     sumq.v = (sumq.v+calck0_qq(4).v)*xx.v
                     sump.v = sump.v*xx.v+calck0_pp(6).v
                     sumq.v = (sumq.v+calck0_qq(5).v)*xx.v
                     sump.v = sump.v*xx.v+calck0_pp(7).v
                     sumq.v = (sumq.v+calck0_qq(6).v)*xx.v
                     sump.v = sump.v*xx.v+calck0_pp(8).v
                     sumq.v = (sumq.v+calck0_qq(7).v)*xx.v
                     sump.v = sump.v*xx.v+calck0_pp(9).v
                     sumq.v = (sumq.v+calck0_qq(8).v)*xx.v
                     sumq.v = sumq.v+calck0_qq(9).v
                     t0.v   = sqrt(x.v)
                     t1.v   = sump.v/sumq.v
                     val.v  = t1.v/t0.v
                     if(jint==1) val.v = val.v*exp(-x.v)
                 end if
             else
                 val.v = xinf.v
             end if
#endif
         end subroutine calck0_ymm4r8
         
         
#if 0
/*
!*****************************************************************************80
!
!! CALCK1 computes various K1 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order one, K1(X) and EXP(X)*K1(X), for real arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  XLEAST < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K1(x);
!    2, RESULT = exp(x) * K1(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K1(x);
!    2, exp(x) * K1(x);
*/      
#endif    


          subroutine calck1_ymm4r8(arg,jint,val)
               
              !dir$ optimize:3
              !dir$ attributes code_align : 32 :: calck1_ymm4r8
              !dir$ attributes forceinline :: calck1_ymm4r8
              !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: calck1_ymm4r8  
              type(YMM4r8_t),   intent(in)   :: arg
              type(YMM4r8_t),   intent(out)  :: val
              integer(kind=i4), intent(in)   :: jint  
              !dir$ attributes align : 64 :: one
              !dir$ attributes align : 64 :: zero
              !dir$ attributes align : 64 :: xleast
              !dir$ attributes align : 64 :: xsmall
              !dir$ attributes align : 64 :: xinf
              !dir$ attributes align : 64 :: xmax
              !dir$ attributes align : 64 :: sumf
              !dir$ attributes align : 64 :: sumg
              !dir$ attributes align : 64 :: sump
              !dir$ attributes align : 64 :: sumq
              !dir$ attributes align : 64 :: x
              !dir$ attributes align : 64 :: xx
              !dir$ attributes align : 64 :: t0
              !dir$ attributes align : 64 :: t1
              !dir$ attributes align : 64 :: t2
              type(YMM4r8_t),   parameter :: one   = YMM4r8_t(1.0e+0_dp)
              type(YMM4r8_t),   parameter :: zero  = YMM4r8_t(0.0e+0_dp)
              type(YMM4r8_t),   parameter :: xleast= YMM4r8_t(2.23e-308_dp)
              type(YMM4r8_t),   parameter :: xsmall= YMM4r8_t(1.11e-16_dp)
              type(YMM4r8_t),   parameter :: xinf  = YMM4r8_t(1.79e+308_dp)
              type(YMM4r8_t),   parameter :: xmax  = YMM4r8_t(705.343e+0_dp)
              type(YMM4r8_t),   automatic :: sumf,sumg
              type(YMM4r8_t),   automatic :: sump,sumq
              type(YMM4r8_t),   automatic :: x,xx
              type(YMM4r8_t),   automatic :: t0,t1
              type(YMM4r8_t),   automatic :: t2 
              type(Mask4_t),    automatic :: m0,m1,m2 
#if (GMS_EXPLICIT_VECTORIZE) == 1
               integer(kind=i4) :: j
#endif   
#if (GMS_EXPLICIT_VECTORIZE) == 1
           
              !dir$ loop_count(4)
              !dir$ vector aligned
              !dir$ vector vectorlength(8)
              !dir$ vector always
              do j=0,3
                    x.v(j)  = arg.v(j) 
                    m0.m(j) = (x.v(j)<xleast.v(j))
                    m1.m(j) = (x.v(j)<=one.v(j))
                    if(all(m0.m(j))) then
                        val.v(j) = xinf.v(j)
                        m0.m(j) = (xmax.v(j)<x.v(j))
                    else if(m1.m)) then
                        m2.m(j) = (x.v(j)<small.v(j))
                        if(all(m2.m(j))) then
                            val.v(j) = one.v(j)/x.v(j)
                        else
                            xx.v(j)   = x.v(j)*x.v(j)
                            sump.v(j) = ((((                &
                                    calck1_p(0).v(j)  &
                             * xx.v(j)+calck1_p(1).v(j)) &
                             * xx.v(j)+calck1_p(2).v(j)) &
                             * xx.v(j)+calck1_p(3).v(j)) &
                             * xx.v(j)+calck1_p(4).v(j)) &
                             * xx.v(j)+calck1_q(2).v(j)
                             sumq.v(j) = ((                  &
                               xx.v(j)+calck1_q(0).v(j)) &
                             * xx.v(j)+calck1_q(1).v(j)) &
                             * xx.v(j)+calck1_q(2).v(j))
                             t2.v(j)   = sump.v(j)/sumq.v(j)
                             sumf.v(j) = (((                 &
                                    calck1_f(0).v(j)  &
                             * xx.v(j)+calck1_f(1).v(j)) &
                             * xx.v(j)+calck1_f(2).v(j)) &
                             * xx.v(j)+calck1_f(3).v(j))
                             * xx.v(j)+calck1_f(4).v(j)
                             sumg.v(j) = ((                  &
                                    xx.v(j)+calck1_g(1).v(j)) &
                                  * xx.v(j)+calck1_g(2).v(j)) &
                                  * xx.v(j)+calck1_g(3).v(j)
                             t0.v(j)   = xx.v(j)+log(x.v(j))
                             t1.v(j)   = sumf.v(j)/sumg.v(j)
                             val.v(j)  = (t0.v(j)*t1.v(j)+t2.v(j))/x.v(j)
                             if(jint==2) val.v(j) = val.v(j)*exp(x.v(j))
                       end if
              else if(jint==1.and.all(m0.m(j))) then
                      val.v(j) = zero.v(j)
              else
                      xx.v(j)  = one.v(j)/x.v(j)
                      sumq.v(j)= xx.v(j)
                      sump.v(j)= calck1_pp(0).v(j)
                      sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(1).v(j)
                      sumq.v(j)= (sumq.v(j)+calck1_qq(0).v(j))*xx.v(j)
                      sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(2).v(j)
                      sumq.v(j)= (sumq.v(j)+calck1_qq(1).v(j))*xx.v(j)
                      sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(3).v(j)
                      sumq.v(j)= (sumq.v(j)+calck1_qq(2).v(j))*xx.v(j)
                      sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(4).v(j)
                      sumq.v(j)= (sumq.v(j)+calck1_qq(3).v(j))*xx.v(j)
                      sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(5).v(j)
                      sumq.v(j)= (sumq.v(j)+calck1_qq(4).v(j))*xx.v(j)
                      sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(6).v(j)
                      sumq.v(j)= (sumq.v(j)+calck1_qq(5).v(j))*xx.v(j)
                      sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(7).v(j)
                      sumq.v(j)= (sumq.v(j)+calck1_qq(6).v(j))*xx.v(j)
                      sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(8).v(j)
                      sumq.v(j)= (sumq.v(j)+calck1_qq(7).v(j))*xx.v(j)
                      sump.v(j)= sump.v(j)*xx.vv+calck1_pp(9).v(j)
                      sump.v(j)= sump.v(j)*xx.v(j)+calck1_pp(10).v(j)
                      sumq.v(j)= sumq.v(j)+calck1_qq(8).v(j)
                      t0.v(j)  = sqrt(x.v(j))
                      t1.v(j)  = sumq.v(j)/sump.v(j)
                      if(jint==1) val.v(j) = val.v(j)*exp(-x.v(j))
               end if
          end do
#else
              x.v  = arg.v
              m0.m = (x.v<xleast.v)
              m1.m = (x.v<=one.v)
              if(all(m0.m)) then
                  val.v = xinf.v
                  m0.m = (xmax.v<x.v)
              else if(m1.m)) then
                  m2.m = (x.v<small.v)
                  if(all(m2.m)) then
                      val.v = one.v/x.v
                  else
                      xx.v   = x.v*x.v
                      sump.v = ((((                &
                                    calck1_p(0).v  &
                             * xx.v+calck1_p(1).v) &
                             * xx.v+calck1_p(2).v) &
                             * xx.v+calck1_p(3).v) &
                             * xx.v+calck1_p(4).v) &
                             * xx.v+calck1_q(2).v
                      sumq.v = ((                  &
                               xx.v+calck1_q(0).v) &
                             * xx.v+calck1_q(1).v) &
                             * xx.v+calck1_q(2).v)
                      t2.v   = sump.v/sumq.v
                      sumf.v = (((                 &
                                    calck1_f(0).v  &
                             * xx.v+calck1_f(1).v) &
                             * xx.v+calck1_f(2).v) &
                             * xx.v+calck1_f(3).v)
                             * xx.v+calck1_f(4).v
                      sumg.v = ((                  &
                                    xx.v+calck1_g(1).v) &
                                  * xx.v+calck1_g(2).v) &
                                  * xx.v+calck1_g(3).v
                      t0.v   = xx.v+log(x.v)
                      t1.v   = sumf.v/sumg.v
                      val.v  = (t0.v*t1.v+t2.v)/x.v
                      if(jint==2) val.v = val.v*exp(x.v)
                  end if
              else if(jint==1.and.all(m0.m)) then
                  val.v = zero.v
              else
                  xx.v  = one.v/x.v
                  sumq.v= xx.v
                  sump.v= calck1_pp(0).v
                  sump.v= sump.v*xx.v+calck1_pp(1).v
                  sumq.v= (sumq.v+calck1_qq(0).v)*xx.v
                  sump.v= sump.v*xx.v+calck1_pp(2).v
                  sumq.v= (sumq.v+calck1_qq(1).v)*xx.v
                  sump.v= sump.v*xx.v+calck1_pp(3).v
                  sumq.v= (sumq.v+calck1_qq(2).v)*xx.v
                  sump.v= sump.v*xx.v+calck1_pp(4).v
                  sumq.v= (sumq.v+calck1_qq(3).v)*xx.v
                  sump.v= sump.v*xx.v+calck1_pp(5).v
                  sumq.v= (sumq.v+calck1_qq(4).v)*xx.v
                  sump.v= sump.v*xx.v+calck1_pp(6).v
                  sumq.v= (sumq.v+calck1_qq(5).v)*xx.v
                  sump.v= sump.v*xx.v+calck1_pp(7).v
                  sumq.v= (sumq.v+calck1_qq(6).v)*xx.v
                  sump.v= sump.v*xx.v+calck1_pp(8).v
                  sumq.v= (sumq.v+calck1_qq(7).v)*xx.v
                  sump.v= sump.v*xx.v+calck1_pp(9).v
                  sump.v= sump.v*xx.v+calck1_pp(10).v
                  sumq.v= sumq.v+calck1_qq(8).v
                  t0.v  = sqrt(x.v)
                  t1.v  = sumq.v/sump.v
                  if(jint==1) val.v = val.v*exp(-x.v)
               end if
#endif
          end subroutine calck1_ymm4r8 
          
          
          

end module spec_funcs_ymm4r8
