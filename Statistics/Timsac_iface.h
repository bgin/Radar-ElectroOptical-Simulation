
#ifndef __TIMSAC_IFACE_H__
#define __TIMSAC_IFACE_H__







#if defined(__cplusplus)

extern "C"  {


	/*
			
			    THIS PROGRAM COMPUTES POWER SPECTRUM ESTIMATES FOR TWO
			    TRIGONOMETRIC WINDOWS OF BLACKMAN-TUKEY TYPE BY GOERTZEL METHOD.
			    ONLY ONE CARD OF LAGH(MAXIMUM LAG OF COVARIANCES TO BE USED FOR
			    POWER SPECTRUM COMPUTATION) SHOULD BE ADDED ON TOP OF THE OUTPUT
			    OF PROGRAM 5.1.1 AUTCOR TO FORM INPUT TO THIS PROGRAM.
			    OUTPUTS ARE ESTIMATES P1(I),P2(I) FOR FREQUENCIES I/(2LAGH*DELTAT)
			    AND THE TEST STATISTICS Q(I) FOR THE DIFFERENCES BETWEEN P1(I) AND
			    P2(I).   Q(I) GREATER THAN 1 MEANS SIGNIFICANT DIFFERENCE.
	*/
	 void auspecf_(int *, 
			      int *, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict );

	 /*

			     THIS PROGRAM PROVIDES AN AUTOMATIC AR-MA MODEL FITTING PROCEDURE.
				 MODELS WITH VARIOUS ORDERS ARE FITTED AND THE BEST CHOICE IS DETER
		         WITH THE AID OF THE STATISTICS AIC.
				 THE MAXIMUM LIKELIHOOD ESTIMATES OF THE COEFFICIENTS OF A SCALAR
			     AUTOREGRESSIVE MOVING AVERAGE MODEL Y(I)+B(1)Y(I-1)+...+B(IQ)Y(I-I
			     =X(I)+A(1)X(I-1)+...+A(IP)X(I-IP) OF A TIME SERIES Y(I)
				 ARE OBTAINED BY USING DAVIDON'S VARIANCE ALGORITHM.
				 PURE AUTOREGRESSION IS NOT ALLOWED.
				 FOR AR-MODELS USE THE INTERMEDIATE OUTPUTS OF CANARM.
	 */
	 void autarmf_(
				     int *,
					 int *,
					 double * __restrict,
					 int *,
					 int * __restrict,
					 double * __restrict,
					 int * __restrict,
					 double * __restrict,
					 int *,
					 int * __restrict,
					 double * __restrict,
					 int * __restrict,
					 double * __restrict,
					 double * __restrict,
					 double * __restrict,
					 double * __restrict,
					 double * __restrict,
					 double *,
					 int *,
					 int *,
					 int *,
					 int *,
					 int * );


	 /*
				     THIS PROGRAM REQUIRES FOLLOWING INPUTS:
				     N: LENGTH OF DATA
					 LAGH: MAXIMUM LAG
					 DFORM: INPUT FORMAT SPECIFICATION STATEMENT IN ONE CARD,
					 FOR EXAMPLE
					 (8F10.4)
					 (X(I),I=1,N): ORIGINAL DATA.
					 THE OUTPUTS ARE AUTOCOVARIANCES (CXX(I); I=0,LAGH) AND
					 AUTO CORRELATIONS (NORMALIZED COVARIANCES).
	 */
	 void autcorf_(double * __restrict, 
				  int *, 
				  double * __restrict, 
				  double * __restrict, 
				  int *, 
				  double *);

	 /*
				     THIS PROGRAM REALIZES A DECOMPOSITION OF TIME SERIES Y            
					 INTO THE FORM                                                     
					 Y(I) = T(I) +S(I)+I(I)+TDC(I)+OCF(I)                              
					 WHERE  T(I)=TREND  S(I)=SEASONAL  I(I)=IRREGULAR                  
				     TDC(I)=TRADING DAY COMPONENT     AND                       
					 OCF(I)=OUTLIER CORRECTION FACTOR                           
                                                                       
				     THE PROCEDURE IS BASED ON A BAYESIAN MODEL AND ITS                
					 PERFORMANCE IS CONTROLLED BY THE SELECTION OF THE PARAMETERS OF   
					 THE PRIOR DISTRIBUTION.  
	 */
	 void bayseaf_(double * __restrict, 
				  int *, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict, 
				  double * ,
				  int * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  int *, 
				  int *, 
				  int * );

	 /*
				     THIS PROGRAM COMPUTES BISPECTRUM USING THE DIRECT FOURIER TRANSFOR
			         OF SAMPLE THIRD ORDER MOMENTS.
					 THIS PROGRAM REQUIRES THE FOLLOWING INPUTS;
					 OUTPUTS OF THE PROGRAM THIRMO:
					 N; DATA LENGTH,
					 MH; MAXIMUM LAG,
					 CC(I); AUTOCOVARIANCES,
					 C(I,J); THIRD ORDER MOMENTS.   
	 */
	 void bispecf_(int *, 
			      int *, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double *);

	 /*
					   BAYESIAN METHOD OF LOCALLY STATIONARY AR MODEL FITTING; SCALAR CAS
                                                                       
                       THIS PROGRAM LOCALLY FITS AUTOREGRESSIVE MODELS TO NON-STATIONARY 
                       SERIES BY A BAYESIAN PROCEDURE.  POWER SPECTRA FOR STATIONARY SPAN
                       ARE GRAPHICALLY PRINTED OUT.  (THIS PROGRAM IS TENTATIVE.)        
                                                                       
                       INPUTS REQUIRED:                                                  
                       MT:       INPUT DEVICE FOR ORIGINAL DATA (MT=5 : CARD READ
                       LAG:      UPPER LIMIT OF THE ORDER OF AR MODEL, MUST BE LE
                                 OR EQUAL TO 50.                                 
                       NS:       LENGTH OF BASIC LOCAL SPAN                      
                       KSW:      =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESS
                                 =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REG
                                                                       
                       -- THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDA
                       TITLE:    SPECIFICATION OF DATA                           
                       N:        DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000
                       DFORM:    INPUT DATA SPECIFICATION STATEMENT.             
                       -- EXAMPLE  --     (8F10.5)                     
                       (Z(I),I=1,N):  ORIGINAL DATA         
	 */
	 void blocarf_(double * __restrict, 
				  int * , 
				  int *, 
				  int *, 
				  int *, 
				  double *, 
				  double *,
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  int * __restrict, 
				  int * __restrict,
				  double * __restrict);

	 /*
					    BAYESIAN METHOD OF LOCALLY STATIONARY MULTIVARIATE AR MODEL FITTIN
                                                                       
					    THIS PROGRAM LOCALLY FITS MULTI-VARIATE AUTOREGRESSIVE MODELS TO  
						NON-STATIONARY TIME SERIES BY A BAYESIAN PROCEDURE.               
                                                                       
     
					    THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
                        MRDATA                                                    
                        MNONSB                                                    
                        --------------------------------------------------------------- 
                        INPUTS REQUIRED;                                                
                        MT:    INPUT DEVICE FOR ORIGINAL DATA (MT=5: CARD READER).   
                        LAG:   UPPER LIMIT OF THE ORDER OF AR-MODEL, MUST BE LESS THA
                               OR EQUAL TO 50.                                       
                        NS:    LENGTH OF BASIC LOCAL SPAN.                           
                        KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR    
                               =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR
                                                                        
                               -- THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE MRDATA 
                        TITLE: SPECIFICATION OF DATA                                 
                        N:     DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 1000.      
                        ID:    DIMENSION OF DATA,  MUST BE LESS THAN 6               
                               < ID*(LAG+1)+KSW MUST BE LESS THAN 101 >        
                        IFM:   INPUT FORMAT                                          
                        FORM:  INPUT DATA FORMAT SPECIFICATION STATEMENT.            
                              -- EXAMPLE --     (8F10.5)                            
                        C(J):  CALIBRATION CONSTANT FOR CHANNEL J (J=1,ID)           
                        Z(I,J): ORIGINAL DATA            
	 */
	 void blomarf_(double * __restrict, 
				  int *, 
				  int *, 
				  double *, 
				  int *, 
				  int *, 
				  int *, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  int * __restrict, 
				  int * __restrict, 
				  int *);

	 /*
						THIS PROGRAM PRODUCES BAYESIAN ESTIMATES OF TIME SERIES MODELS SUC
					    PURE AR MODELS, AR-MODELS WITH NON-LINEAR TERMS, AR-MODELS WITH PO
						TYPE MEAN VALUE FUNCTIONS, ETC.  THE GOODNESS OF FIT OF A MODEL IS
					    CHECKED BY THE ANALYSIS OF SEVERAL STEPS AHEAD PREDICTION ERRORS. 
						BY PREPARING AN EXTERNAL SUBROUTINE SETX PROPERLY, ANY TIME SERIES
					    WHICH IS LINEAR IN PARAMETERS CAN BE TREATED.                     
                        ----------------------------------------------------------------
                        THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
                        REDATA                                                    
                        REDLAG                                                    
                        SETLAG                                                    
                        REDREG                                                    
                        REDUCT                                                    
                        ARMFIT                                                    
                        SBBAYS                                                    
                        CHECK                                                     
                        ----------------------------------------------------------------
                                                                       
						INPUTS REQUIRED:                                                
                        MT:    ORIGINAL DATA INPUT DEVICE SPECIFICATION              
                        IMODEL:=1  AUTOREGRESSIVE MODEL                              
                               =2  POLYNOMIAL TYPE NON-LINEAR MODEL (LAG'S READ IN ) 
                               =3  POLYNOMIAL TYPE NON-LINEAR MODEL (LAG'S AUTOMATICA
                               =4  AR-MODEL WITH POLYNOMIAL MEAN VALUE FUNCTION      
                               =5  ANY NON-LINEAR MODEL                              
                               =6  POLYNOMIAL TYPE EXPONENTIALLY DAMPED NON-LINEAR MO
                               =7  THIS MODEL IS RESERVED FOR THE USER'S OPTIONAL USE
                       LAG:   MAXIMUM TIME LAG USED IN THE MODEL                    
                       K:     NUMBER OF REGRESSORS                                  
                       IL:    PREDICTION ERRORS CHECKING (UP TO IL-STEPS AHEAD) IS R
                              N*IL SHOULD BE LESS THAN OR EQUAL TO 20000            
                                                                       
                              THE FOLLOWING INPUTS ARE REQUIRED AT SUBROUTINE REDATA   --
                                                                       
                      TITLE:   ORIGINAL DATA SPECIFICATION                         
                      N:       DATA LENGTH                                         
                      DFORM:   INPUT DATA FORMAT SPECIFICATION STATEMENT           
                               -- EXAMPLE --  (8F10.5 )                            
                               X(I) (I=1,N):   ORIGINAL DATA            
	 */
	 void bsubstf_(
								 double * __restrict,
								 int *,
								 int *,
								 int *,
								 int *,
								 int *,
								 int *__restrict,
								 int *__restrict,
								 double *,
								 double *,
								 int *,
								 double *,
								 double *,
								 double * __restrict,
								 double * __restrict,
								 double * __restrict,
								 double * __restrict,
								 double *,
								 double *,
								 double *,
								 double * __restrict,
								 int * __restrict,
								 double * __restrict,
								 double * __restrict,
								 double * __restrict,
								 double * __restrict,
								 double *,
								 double * __restrict,
								 double *,
								 double *,
								 double * __restrict,
								 double * __restrict,
								 double * __restrict,
								 double * __restrict,
								 double * __restrict,
								 double * __restrict,
								 double * __restrict );

	 /*
				     THIS PROGRAM FITS AN AR-MA MODEL TO STATIONARY SCALAR TIME SERIES
					 THROUGH THE ANALYSIS OF CANONICAL CORRELATIONS
					 BETWEEN THE FUTURE AND PAST SETS OF OBSERVATIONS.
					 THE OUTPUTS OF THIS PROGRAM SHOULD BE ADDED TO THE INPUTS
					 TO THIS PROGRAM TO FORM AN INPUT TO THE PROGRAM AUTARM.
    
        		     INPUTS REQUIRED:
					 (N,LAGH0): N, LENGTH OF ORIGINAL DATA Y(I) (I=1,N)
							    LAGH0, MAXIMUM LAG OF COVARIANCE
					 CYY(I),I=0,LAGH0: AUTOCOVARIANCE SEQUENCE OF Y(I)
    
				     OUTPUTS:
                     NEWL: NEWL=1, FOR DIRECT INPUT TO PROGRAM AUTARM
                     M1M: ORDER OF AR
                     BETA(I)(I=1,M1M): AR-COEFFICIENTS
                     M1N: ORDER OF MA (=M1M-1)
                     ALPHA(I)(I=1,M1N): MA-COEFFICIENTS
    
                    THE AR-MA MODEL IS GIVEN BY
                    Y(N)+BETA(1)Y(N-1)+...+BETA(M1M)Y(N-M1M) = X(N)+ALPHA(1)X(N-1)+...
                                                 ...+ALPHA(M1N)X(N-M1N)    
	 */
	 void canarmf_(int *, 
				  int *, 
				  double * __restrict, 
				  double * __restrict, 
				  int *, 
				  double * __restrict,
			      double * __restrict, 
				  double *, 
				  int *, 
				  double * __restrict, 
				  int *, 
				  int * __restrict,
				  int * __restrict, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  int * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  int * __restrict, 
				  int *, 
				  double * __restrict, 
				  int *, 
				  double * __restrict, 
				  int *,
				  int *);

	 /*
					THIS PROGRAM DOES CANONICAL CORRELATION ANALYSIS OF AN IR-DIMENSIO
				    MULTIVARIATE TIME SERIES Y(I) (I=1,N).
    
				    FIRST AR-MODEL IS FITTED BY THE MINIMUM  A I C  PROCEDURE.
					THE RESULTS ARE USED TO ORTHO-NORMALIZE THE PRESENT AND PAST VARIA
					THE PRESENT AND FUTURE VARIABLES ARE TESTED SUCCESSIVELY TO DECIDE
					ON THE DEPENDENCE OF THEIR PREDICTORS. WHEN THE LAST DIC (AN INFOR
					CRITERION) IS NEGATIVE THE PREDICTOR OF THE VARIABLE IS DECIDED
					TO BE LINEARLY DEPENDENT ON THE ANTECEDENTS. 
				    THE STRUCTURAL CHARACTERISTIC VECTOR H OF THE CANONICAL MARKOVIAN
					REPRESENTATION AND THE ESTIMATE OF THE TRANSITION MATRIX F, IN
					VECTOR FORM, ARE PUNCHED OUT. THE ESTIMATE OF THE INPUT MATRIX G A
					THE COVARIANCE MATRIX C OF THE INNOVATION, OBTAINED BY USING
					THE F-MATRIX AND THE AR-MODEL, ARE ALSO PUNCHED OUT.
	 */
	 void canocaf_(int * , 
				  int * __restrict, 
				  int *, 
				  int *, 
				  int *, 
				  double * __restrict, 
				  int *,
			      double * __restrict, 
				  double * , 
				  int *, 
				  double * __restrict, 
				  double * __restrict,
				  int *, 
				  int * __restrict, 
				  int * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  int * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  int *, 
				  int * __restrict, 
				  double * __restrict,
				  int *, 
				  double * __restrict, 
				  int *, 
				  int *, 
				  int *);

	 /*
				     THIS PROGRAM PRODUCES THE FOURIER TRANSFORM OF A POWER
			         GAIN FUNCTION IN THE FORM OF AN AUTOCOVARIANCE SEQUENCE.
				     THE GAIN FUNCTION IS DEFINED AS A RECTILINEAR FUNCTION WITH
			         THE VALUES G(I) SPECIFIED AT THE FREQUENCIES F(I),I=1,K.
				     THE OUTPUTS OF THIS PROGRAM ARE USED AS THE INPUTS TO THE CANONICA
				     CORRELATION ANALYSIS PROGRAM CANARM, TO REALIZE A FILTER WITH
			         THE DESIRED GAIN FUNCTION.
    
			         THE FOLLOWING INPUTS ARE REQUIRED:
			         (L,K): L, DESIRED MAXIMUM LAG OF COVARIANCE (AT MOST 1024)
		              K, NUMBER OF DATA POINTS (LESS THAN OR EQUAL TO 500)
			         (F(I),G(I))(I=1,K): F(I), FREQUENCY. BY DEFINITION F(1)=0.0 AND F(
                      F(I)'S ARE ARRANGED IN INCREASING ORDER.
                     G(I), POWER GAIN OF THE FILTER AT THE FREQUEN
    
				     OUTPUTS:
                     (N,LAGH): N=2048
			         C(I)(I=0,LAGH): 
	 */
	 void covgenf_(int *, 
				  int *, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict );

	 /*
					   ...  TIME SERIES DECOMPOSITION (SEASONAL ADJUSTMENT) ...             
                                                                       
                       THE BASIC MODEL:                                                  
                                                                       
                              y(n) = T(n) + AR(n) + S(n) + TD(n) + R(n) + W(n)             
                                                                       
                              where                                                           
                              T(n):       trend component                                
                              AR(n):      AR process                                     
                              S(n):       seasonal component                             
                              TD(n):      trading day factor                             
                              R(n):       any other explanetory variables                
                              W(n):       observational noise                            
                                                                       
                     COMPONENT MODELS:                                                 
                                                                       
                     Trend component                                                 
                     T(N) =  T(N-1) + V1(N)                             :M1 = 1 
                     T(N) = 2T(N-1) - T(N-2) + V1(N)                    :M1 = 2 
                     T(N) = 3T(N-1) -3T(N-2) + T(N-2) + V1(N)           :M1 = 3 
                                                                       
                     AR componet:                                                    
                     AR(n) = a(1)AR(n-1) + ... + a(m2)AR(n-m2) + V2(n)          
                                                                       
                     Seasonal component:                                             
                     S(N) =  -S(N-1) - ... - S(N-PERIOD+1) + V3(N)    :SORDER=1 
                     S(N) = -2S(N-1) - ... -PERIOB*S(N-PERIOD+1)      :SORDER=2 
                            - ... - S(n-2PERIOD+2) + V3(n)             
                     Trading day effect:                                             
                     TD(n) = b(1)*TRADE(n,1) + ... + b(7)*TRADE(n,7)            
                                                                      
                     TRADE(n,i):  number of i-th days of the week in n-th data  
                     b(1) + ... + b(7) = 0           
	 */
	 void decompf_(double * __restrict, 
				  int *, 
				  int * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  int *, 
				  double *, 
				  int * );

	 /*
					  EXACT MAXIMUM LIKELIHOOD METHOD OF SCALAR AR-MODEL FITTING        
                                                                       
					  THIS PROGRAM PRODUCES EXACT MAXIMUM LIKELIHOOD ESTIMATES OF THE   
                      PARAMETERS OF A SCALAR AR-MODEL.                                  
                                                                       
                      THE AR-MODEL IS GIVEN BY                                          
                                                                       
                      Z(I) = A(1)*Z(I-1) + ... + A(K)*Z(I-K) + E(I)           
                                                                       
                      WHERE E(I) IS A ZERO MEAN WHITE NOISE.                            
                                                                       
					 --------------------------------------------------------------    
					 THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:    
                     REDATA                                                    
                     REDUCT                                                    
                     ARMFIT                                                    
                     RECOEF                                                    
                     ARMLE                                                     
                     PRINTA                                                    
					 --------------------------------------------------------------    
					 INPUTS REQUIRED:                                                  
                     MT:      INPUT DEVICE SPECIFICATION (MT=5: CARD READER)      
                     LAG:     UPPER LIMIT OF AR-ORDER, MUST BE LESS THAN 51       
                                                                       
                 --  THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDATA  --   
                     TITLE:  TITLE OF DATA                                        
                     N:      DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000     
                     DFORM:  INPUT DATA FORMAT SPECIFICATION STATEMENT            
                 -- EXAMPLE --     (8F10.5)                           
                    (Z(I),I=1,N):  ORIGINAL DATA           
	 */
	 void exsarf_(double * __restrict, 
				 int *, 
				 int *, 
				 double *, 
				 double *, 
				 double * __restrict,
				 double * __restrict, 
				 double * __restrict, 
				 int *, 
				 double *, 
				 double *,
				 double * __restrict, 
				 double *, 
				 double * __restrict, 
				 int *);

	 /*
				    THIS PROGRAM COMPUTES AUTO AND/OR CROSS
					COVARIANCES AND CORRELATIONS VIA FFT.
					IT REQUIRES FOLLOWING INPUTS:
					ISW: ISW=1...AUTO CORRELATION OF X (ONE-CHANNEL)
					     ISW=2...AUTO CORRELATIONS OF X AND Y (TWO-CHANNEL)
					     ISW=4...AUTO,CROSS CORRELATIONS OF X AND Y (TWO-CHANNEL)
			        LD: LENGTH OF DATA
                    LAGH: MAXIMUM LAG
                    DFORM: INPUT FORMAT SPECIFICATION STATEMENT IN ONE CARD,
                    FOR EXAMPLE
                   (8F10.4)
                   (X(I); I=1,LD): DATA OF CHANNEL X
                   (Y(I); I=1,LD): DATA OF CHANNEL Y (FOR ISW=2 OR 4 ONLY)
	 */
	 void fftcorf_(int *, 
				  int *, 
				  int *, 
				  int *, 
				  int *, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict);

	 /*
				    THIS PROGRAM PERFORMS FPE(FINAL PREDICTION ERROR) COMPUTATION FOR
				    ONE-DIMENSIONAL AR-MODEL. A CARD CONTAINING THE FOLLOWING
				     !INFORMATION OF L, UPPER LIMIT OF MODEL ORDER, SHOULD BE ADDED ON
                    TOP OF THE OUTPUT OF PROGRAM 5.1.1 AUTCOR TO FORM THE INPUT TO
                    THIS PROGRAM.
                    CXX(0) IS READ AS INITIAL SD.
                    THE OUTPUTS ARE THE COEFFICIENTS A(I) OF AR-PROCESS
                    X(N)=A(1)X(N-1)+...+A(M)X(N-M)+E(N)
                    AND THE VARIANCE SIGMA**2 OF E(N).
                    CHI**2 SHOWS THE SIGNIFICANCE OF PARCOR=A(M) AS A CHI-SQUARED
                    VARIABLE WITH D.F.=1.
	 */
	 void fpeautf_(int *, 
				  int *, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict, 
				  double * __restrict,
				  double *, 
				  double *, 
				  double *, 
				  int *, 
				  double *, 
				  double * __restrict, 
				  double * __restrict );

	 /*
				     THIS PROGRAM PERFORMS FPEC(AR-MODEL FITTING FOR CONTROL)
					 COMPUTATION.
					 BESIDES THE OUTPUTS OF PROGRAM 5.1.2   MULCOR, THE FOLLOWING
					 INPUTS ARE REQUIRED:
				     L: UPPER LIMIT OF MODEL ORDER M (LESS THAN 30)
					 IR: NUMBER OF CONTROLLED VARIABLES
					 IL: NUMBER OF MANINPULATED VARIABLES, IL=0 FOR MFPE COMPUTATION
					 INW(I): INDICATOR; FIRST IR INDICATE THE CONTROLLED VARIABLES
							 AND THE REST THE MANIPULATE VARIABLES WITHIN THE IP0 VARIABLES
					 IN THE OUTPUT OF PROGRAM 5.1.2   MULCOR.
					 THE OUTPUTS ARE THE PREDICTION ERROR COVARIANCE MATRIX OSD AND
					 THE SET OF COEFFICIENT MATRICES A AND B TO BE USED IN
					 PROGRAM 5.5.1   OPTIMAL CONTROLLER DESIGN.
	 */
	 void fpec7f_(int *, 
				 int *, 
				 int *, 
				 int *, 
				 int *, 
				 int * __restrict,
			     double * __restrict, 
				 double * __restrict, 
				 double * __restrict,
				 double * __restrict, 
				 double * __restrict, 
				 int *, 
				 double *,
				 double *, 
				 double *, 
				 double * __restrict, 
				 double * __restrict);

	 /*
					     IN THIS PROGRAM THE MATRICES A,B,C, STAND FOR TRANSITION MATRIX (F
						 ! INPUT MATRIX (G), OUTPUT MATRIX (H), RESPECTIVELY.
						 !
					     ! THE INPUTS REQUIRED ARE AS FOLLOWS:
					     ! (N,LAGH0,ID0):
						 !     N, LENGTH OF ORIGINAL DATA
						 !     LAGH0, MAXIMUM LAG OF COVARIANCE
					     !     ID0, DIMENSION OF Y(I)
					     ! (CYY(I)(I=0,LAGH0): COVARIANCE MATRIX SEQUENCE OF Y(I). CYY(I) ARE
						 ! !                      THE OUTPUTS OF THE PROGRAM MULCOR OF TIMSAC AND
					     !                       ARE USED AS THE INPUT TO THE PROGRAM CANOCA OF
						 !                       TIMSAC-74. FOR THE CONSTRUCTION OF CYY(I),
					     !                       SEE PROGRAM CANOCA.
					     ! THE OUTPUTS OF PROGRAM CANACA:
						 !      (ID,K):
					     !          ID, DIMENSION OF THE TIME SERIES Y(I) (NOT GREATER THAN 5)
						 !          K, DIMENSION OF THE STATE VECTOR (NOT GREATER THAN 10)
					     !      (NH(I))(I=1,K): STRUCTURAL CHARACTERISTIC VECTOR
						 !      (AW(I))(I=1,IAW): INITIAL ESTIMATE OF THE VECTOR OF FREE
						 !                           PARAMETERS IN F (=A)
						 !      B(I,J)(I=1,ID+1;J=1,ID): INITIAL ESTIMATES OF THE FREE
						 !                                   PARAMETERS IN G (=B)
						 ! ICONT: OUTPUT CONTROL
						 !      = 0, FOR AR-MA COEFFICIENTS
						 !      = 1, FOR SIMCON INPUT
						 !      = 2, FOR BOTH
						 !
						 ! THE OUTPUTS OF THIS PROGRAM ARE;
						 ! THE MAXIMUM LIKELIHOOD ESTIMATES
						 !      (ID,K):
						 !      (NH(I))(I=1,K):
						 !      (AW(I))(I=1,IAW):
						 !      (B(I,J))(I=1,ID+1;J=1,ID):
						 ! AND THE OUTPUTS TO BE USED BY THE PROGRAMS PRDCTR AND SIMCON
					     !      (ID,Q,Q-1):
						 !      (B(I,J,L))(I,J=1,ID,L=1,Q): AR-COEFFICIENT MATRICES
						 !      AND, WHEN ICONT=1 OR 2,
						 !     (W(I,J,L))(I,J=1,ID,L=1,Q-1): IMPULSE RESPONSE MATRICES
						 !      AND, WHEN ICONT=0 OR 1,
						 !      (A(I,J,L))(I,J=1,ID;L=1,Q-1): MA-COEFFICIENT MATRICES
					     ! AND
						 !      (C0(I,J))(I,J=1,ID): INNOVATION COVARIANCE.
	 */
	 void markovf_(int *, 
				  int *, 
				  int *,
				  double * __restrict,
				  int *,
				  int * __restrict,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  int *,
				  int * __restrict,
				  int * __restrict,
				  int * __restrict,
				  int * __restrict,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * ,
				  int *,
				  int *,
				  int *,
				  int * );

	 /*
				'MLOCARF'  
	 */
	 void mlocarf_(double * __restrict,
				  int *,
				  int *,
				  int *,
				  int *,
				  int *,
				  double *,
				  double *,
				  double * __restrict,
				  int * __restrict,
				  double * __restrict,
				  int * __restrict,
				  int * __restrict,
				  double * __restrict,
				  int * __restrict,
				  int * __restrict,
				  int * __restrict,
				  double * __restrict,
				  double * __restrict,
				  int * __restrict,
				  double * __restrict,
				  double * __restrict);

	 /*
			  MINIMUM AIC METHOD OF LOCALLY STATIONARY MULTIVARIATE AR MODEL FIT
			  !                                                                   
			  !     THIS PROGRAM LOCALLY FITS MULTI-VARIATE AUTOREGRESSIVE MODELS TO  
			  !     NON-STATIONARY TIME SERIES BY THE MINIMUM AIC PROCEDURE USING THE 
			  !     HOUSEHOLDER TRANSFORMATION.                                       
			  !                                                                   
			  !     BY THIS PROCEDURE, THE DATA OF LENGTH N ARE DIVIDED INTO J LOCALLY
			  !     STATIONARY SPANS                                                  
			  !                                                                   
			  !            <-- N1 --> <-- N2 --> <-- N3 -->          <-- NJ -->   
			  !           !----------!----------!----------!--------!----------!  
			  !            <-----------------------  N  ---------------------->   
		      !                                                                   
			  !     WHERE NI (I=1,...,J) DENOTES THE NUMBER OF BASIC SPANS, EACH OF   
			  !     LENGTH NS, WHICH CONSTITUTE THE I-TH LOCALLY STATIONARY SPAN.     
			  !     AT EACH LOCAL SPAN, THE PROCESS IS REPRESENTED BY A STATIONARY    
			  !     AUTOREGRESSIVE MODEL. 
	 */
	 void mlomarf_(double * __restrict,
				  int *,
				  int *,
				  double * __restrict,
				  int *,
				  int *,
				  int *,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  int * __restrict,
				  int * __restrict,
				  int * __restrict,
				  double * __restrict,
				  int * __restrict,
				  double * __restrict,
				  int * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  int * __restrict,
				  int * __restrict,
				  int * );

	 /*
				  MULTIVARIATE BAYESIAN METHOD OF AR MODEL FITTING                  
				  !                                                                   
				  ! THIS PROGRAM DETERMINES MULTI-VARIATE AUTOREGRESSIVE MODELS BY A  
				  ! BAYESIAN PROCEDURE.  THE BASIC LEAST SQUARES ESTIMATES OF THE PARA
				  ! ARE OBTAINED BY THE HOUSEHOLDER TRANSFORMATION.                   
				  !                                                                   
				  ! THE STATISTIC AIC IS DEFINED BY                                   
				  !                                                                  
				  !        AIC  =  N * LOG( DET(SD) ) + 2 * (NUMBER OF PARAMETERS)    
				  !                                                                   
				  !   WHERE                                                           
				  !       N:    NUMBER OF DATA,                                       
				  !       SD:   ESTIMATE OF INNOVATION VARIANCE MATRIX                
		          !       DET:  DETERMINANT,                                          
                  !       K:    NUMBER OF FREE PARAMETERS.                            
                  !                                                                   
                  ! BAYESIAN WEIGHT OF THE M-TH ORDER MODEL IS DEFINED BY             
                  !     W(M)  = CONST * C(M) / (M+1)                                  
                  ! WHERE                                                             
                  !     CONST = NORMALIZING CONSTANT                                  
                  !     C(M)  = EXP( -0.5*AIC(M) ).                                   
                  ! THE BAYESIAN ESTIMATES OF PARTIAL AUTOREGRESSION COEFFICIENT MATRI
                  ! OF FORWARD AND BACKWARD MODELS ARE OBTAINED BY (M=1,...,LAG)      
                  !     G(M)  = G(M)*D(M)                                             
                  !     H(M)  = H(M)*D(M),                                            
                  ! WHERE THE ORIGINAL G(M) AND H(M) ARE THE (CONDITIONAL) MAXIMUM    
                  ! LIKELIHOOD ESTIMATES OF THE HIGHEST ORDER COEFFICIENT MATRICES OF 
                  ! FORWARD AND BACKWARD AR MODELS OF ORDER M AND D(M) IS DEFINED BY  
                  !     D(M)  = W(M) + ... + W(LAG).                                  
                  !                                                                   
                  ! THE EQUIVALENT NUMBER OF PARAMETERS FOR THE BAYESIAN MODEL IS     
                  ! !DEFINED BY                                                        
                  !     EK = (D(1)**2 + ... + D(LAG)**2)*ID + ID*(ID+1)/2             
                  ! WHERE ID DENOTES DIMENSION OF THE PROCESS.     
	 */
	 void mulbarf_(double * __restrict,
				  int *,
				  int *,
				  double * __restrict,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  int *,
				  double *,
				  double *,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * );

	 /*
				  PROGRAM 5.1.2   MULTIPLE CORRELATION
                  THIS PROGRAM REQUIRES FOLLOWING INPUTS:
				  ! N: LENGTH OF DATA
				  ! K: DIMENSION OF THE OBSERVATION VECTOR
				  ! LAGH: MAXIMUM LAG
				  ! ISW: ISW=1...ROWWISE DATA INPUT
				  !       ISW=2...COLUMNWISE DATA INPUT
				  ! DFORM: INPUT FORMAT SPECIFICATION STATEMENT IN ONE CARD,
				  ! FOR EXAMPLE
				  ! (8F10.4)
				  ! (X1(S,I); S=1,...,N, I=1,...,K): ORIGINAL DATA MATRIX.
				  ! THE OUTPUTS ARE (CIJ(L): L=0,1,...,LAGH) (I=1,...,K; J=1,...,K),
				  ! WHERE CIJ(L)=COVARIANCE(XI(S+L),XJ(S)),
				  ! !AND THEIR NORMALIZED (CORRELATION) VALUES.   
	 */
	 void mulcorf_(double * __restrict,
				  int *,
				  int *,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict);

	 /*
				   MULTIVARIATE CASE OF MINIMUM AIC METHOD OF AR MODEL FITTING.      
    !                                                                   
    ! THIS PROGRAM FITS A MULTI-VARIATE AUTOREGRESSIVE MODEL BY THE MINI
    ! AIC PROCEDURE.  ONLY THE POSSIBILITIES OF ZERO COEFFICIENTS AT THE
    ! BEGINNING AND END OF THE MODEL ARE CONSIDERED. THE LEAST SQUARES E
    !OF THE PARAMETERS ARE OBTAINED BY THE HOUSEHOLDER TRANSFORMATION. 
    !AIC IS DEFINED BY                                                 
    !                                                                   
    !        AIC  =  N * LOG( DET(SD) ) + 2 * (NUMBER OF PARAMETERS)    
    !                                                                   
    !   WHERE                                                           
    !       N:    NUMBER OF DATA,                                       
    !       SD:   ESTIMATE OF INNOVATION VARIANCE MATRIX                
    !      DET:  DETERMINANT,                                          
    !       K:    NUMBER OF FREE PARAMETERS.                            
    !                                                                  
    !                                                                   
    !   --------------------------------------------------------------- 
    !   THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM.  
    !       MRDATA                                                      
    !       MREDCT                                                      
    !       MARFIT                                                      
    !   --------------------------------------------------------------- 
	 */
	 void mulmarf_(double * __restrict,
				  int * ,
				  int *,
				  double * __restrict,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  int * __restrict,
				  double * __restrict,
				  double * __restrict,
				  int * __restrict,
				  int * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  int *,
				  double * );

	 /*
					 THIS PROGRAM COMPUTES RELATIVE POWER CONTRIBUTIONS IN DIFFERENTIAL     AND INTEGRATED FORM, ASSUMING THE ORTHOGONALITY BETWEEN NOISE
    !     SOURCES.
    !     THE PROGRAM OPERATES ON THE OUTPUT OF PROGRAM 5.3.2 FPEC WITH
    !     IL=0.
    !     THE RESULTS ARE GIVEN AT FREQUIENCIES I/(2*H).
	 */
	 void mulnosf_(int *,
				  int *,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict);

	 /*
			  THIS PROGRAM COMPUTES MULTIPLE SPECTRUM ESTIMATES FROM THE OUTPUT
    ! OF PROGRAM 5.1.2 MULCOR, USING WINDOWS W1 AND W2.
    ! ONLY ONE CARD OF LAGH(MAXIMUM LAG OF COVARIANCES TO BE USED FOR
    ! SPECTRUM COMPUTATION) SHOULD BE ADDED ON TOP OF THE OUTPUT OF
    ! PROGRAM 5.1.2 MULCOR TO FORM THE INPUT TO THIS PROGRAM.
    ! IN THE CARD OUTPUT OF SPECTRUM MATRIX ON AND LOWER DIAGONAL ARE
    ! REAL PARTS AND UPPER DIAGONAL ARE IMAGINARY PARTS OF ON AND LOWER
    ! DIAGONAL SPECTRAL ELEMENTS.
	 */
	 void mulspef_(int *,
				  int *,
				  int *,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict );

	 /*
			THIS PROGRAM LOCALLY FITS AUTOREGRESSIVE MODELS TO NON-STATIONARY
    ! TIME SERIES BY AIC CRITERION.
    ! POWER SPECTRA FOR STATIONARY SPANS ARE GRAPHICALLY PRINTED OUT.
    ! THE FOLLOWING INPUTS ARE REQUIRED;
    !     N: LENGTH OF DATA
    !     ISTP : LENGTH OF THE BASIC LOCAL SPAN
    !     DFORM : INPUT FORMAT SPECIFICATION IN ONE CARD, FOR EXAMPLE,'(8
    !     (X(I),I=1,N) : ORIGINAL DATA.
	 */
	 void nonstf_(int *,
				 int *,
				 double * __restrict,
				 int *,
				 int *,
				 int * __restrict,
				 double * __restrict,
				 double * __restrict,
				 double * __restrict,
				 double * __restrict,
				 double * __restrict,
				 int * __restrict,
				 int * __restrict,
				 double * __restrict );

	 /*
			  THIS PROGRAM OPERATES ON A REAL RECORD OF A VECTOR PROCESS
    !     Y(I) (I=1,N) AND COMPUTES PREDICTED VALUES. ONE STEP AHEAD
    !     PREDICTION STARTS AT TIME P AND ENDS AT TIME Q. PREDICTION IS
    !     CONTINUED WITHOUT NEW OBSERVATIONS UNTIL TIME Q+H.
    !     BASIC MODEL IS THE AUTOREGRESSIVE MOVING AVERAGE
    !     MODEL OF Y(I) WHICH IS GIVEN BY
    !     Y(I)+B(1)Y(I-1)+...+B(K)Y(I-K) = X(I)+A(1)X(I-1)+...+A(L)X(I-L).
    !
    !     THE FOLLOWING INPUTS ARE REQUIRED:
    !     (N,P,Q,H):
    !              N, LENGTH OF DATA
    !              P, ONE STEP AHEAD PREDICTION STARTING POSITION
    !              Q, LONG RANGE FORECAST STARTING POSITION
    !              H, MAXIMUM SPAN OF FORECAST (LESS THAN OR EQUAL TO 100)
    !              (Q+H MUST BE LESS THAN 1001)
    !     JSW: JSW=0 FOR DIRECT LOADING OF AR-MA COEFFICIENTS,
    !              THE OUTPUTS OF PROGRAM MARKOV WITH ICONT=0.
    !        JSW=1 FOR LOADING OF THE OUTPUTS OF PROGRAM MARKOV,
    !              THE OUTPUTS OF PROGRAM MARKOV WITH ICONT=1.
    !    (D,K,L):
    !          D, DIMENSION OF THE VECTOR Y(I)
    !          K, AR-ORDER (LESS THAN OR EQUAL TO 10)
    !          L, MA-ORDER (LESS THAN OR EQUAL TO 10)
    !          N,L,K,H,P,Q,D,JSW,ARE ALL INTEGERS
    !   (DFORM(I),I=1,20): INPUT FORMAT STATEMENT IN ONE CARD,
    !                         FOR EXAMPLE, (8F10.4)
    !    (NAME(I,J),I=1,20,J=1,D): NAME OF THE I-TH COMPONENT
    !    (Y(I,J),I=1,N;J=1,D): ORIGINAL DATA
    !    (B(I1,I2,J),I1=1,D,I2=1,D,J=1,K): AR-COEFFICIENT MATRICES.
    !    FOR JSW=0,
    !         (A(I1,I2,J),I1=1,D,I2=1,D,J=1,L): MA-COEFFICIENT MATRICES.
    !    FOR JSW=1,
    !     (W(I1,I2,J),I1=1,D,I2=1,D,J=1,L): IMPULSE RESPONSE MATRICES.
    !     (S(I,J),I=1,D,J=1,D): INNOVATION VARIANCE MATRIX
    !
    !      THE OUTPUTS OF THIS PROGRAM ARE THE REAL
    !!     AND PREDICTED VALUES OF Y(I).
	 */
	 void prdctrf_(int * ,
				  int *,
				  int *,
				  int *,
				  int *,
				  int *,
				  int *,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict, 
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict );

	 /*
			  EXACT MAXIMUM LIKELIHOOD METHOD OF SCALAR AR-MA MODEL FITTING     
    !                                                                   
    !-----------------------------------------------------------------------
    ! THIS PROGRAM PRODUCES EXACT MAXIMUM LIKELIHOOD ESTIMATES OF THE   
    ! PARAMETERS OF A SCALAR AR-MA MODEL.                               
    !                                                                   
    ! THE AR-MA MODEL IS GIVEN BY                                       
    !                                                                   
    ! Y(I)+B(1)Y(I-1)+...+B(IQ)Y(I-IQ)=X(I)+A(1)X(I-1)+...+A(IP)X(I-IP),
    !                                                                   
    ! WHERE X(I) IS A ZERO MEAN WHITE NOISE.
	 */
	 void xsarmaf_(double * __restrict,
				  int *,
				  int *,
				  int *,
				  double * __restrict,
				  double * __restrict,
				  double *,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double * __restrict,
				  double *,
				  double * );
}


#endif


#endif /*__TIMSAC_IFACE_H__*/
