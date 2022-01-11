    MODULE quadpack2003
      use mod_kinds, only : i4,dp
      IMPLICIT NONE

      PRIVATE ! Key exported symbols are declared PUBLIC

! The components of this derived type define a number of parameters that
! comntrol the approximation process. Extending this defined type allows
! for example, parameters to be passed to the function evaluation routine, 
! the use of recursion to evaluate multiple integrals, etc.

! Input: INTEGER(kind=i4) :: key
!        defines the order of the local Gauss-Kronrod integration rule to be
!        used on each subinterval.
!         <=1: 7 Gauss points, 15 Gauss-Kronrod points,
!           2: 10 Gauss points, 21 Gauss-Kronrod points,
!           3: 15 Gauss points, 31 Gauss-Kronrod points,
!           4: 20 Gauss points, 41 Gauss-Kronrod points,
!           5: 25 Gauss points, 51 Gauss-Kronrod points,
!         >=6: 30 Gauss points, 61 Gauss-Kronrod points.
!        Default value: 6

! Input: INTEGER(kind=i4) :: limit
!        defines the maximum number of subdivisions allowed in the subdivision
!        process. If this value is exceeded the routine will return with an
!        error code (quadpackBase%ier) of 1.
!        Default value: 50

! Input: REAL(dp) :: epsabs, epsrel
!        define the absolute and relative error tolerances required. The routine
!        will attempt to return an approximation satisfying
!             ||I - intApprox|| <= max(epsabs, epsrel*||I||)
!        where I is the exact value of the integral.
!        Default values: epsabs = 0.0_dp, epsrel = sqrt(machine precision)

! Output: INTEGER(kind=i4) :: nfunction_evaluations
!        provides the total number of function evaluations performed.

! Output: INTEGER(kind=i4) :: nvector_evaluations
!        provides the total number of subintervals used in the evaluation.

! Output: REAL(dp) :: abserr
!        provides an estimate of the modulus of the absolute error which should
!        equal or exceed abs(I-intApprox)

! Output: INTEGER(kind=i4) :: ier
!        returns the error status. 0 signifies a successful calculation. The
!        following non-zero values all signify a problem with the computation.
!           1: maximum number of subdivisions allowed has been achieved.  One can
!              allow more subdivisions by increasing the value of LIMIT in QAG.
!              However, if this yields no improvement it is advised to analyze
!              the integrand to determine the integration difficulties.
!              If the position of a local difficulty can be determined, such
!              as a singularity or discontinuity within the interval) one will
!              probably gain from splitting up the interval at this point
!              and calling the integrator on the subranges.  If possible,
!              an appropriate special-purpose integrator should be used
!              which is designed for handling the type of difficulty involved.
!           2: the occurrence of roundoff error is detected, which prevents 
!              the requested tolerance from being achieved.
!           3: extremely bad integrand behavior occurs at some points of the
!              integration interval.
!           4: roundoff error on extrapolation.
!           5: divergent integral (or slowly convergent integral).
!           6: the input is invalid, because EPSABS < 0 and EPSREL < 0.
!           7: limiting number of cycles attained (QAWF only).

! Input/Output: INTEGER(kind=i4) ::ncomponents
!         signifies the number of components in the vector integral. This is
!         used to define the abstract interface for vector integration.
!         It is unused for scalar function integration.

! PROCEDURE(ERRNORM), POINTER :: enorm
!         points to a user-supplied function of the form

!            FUNCTION errnorm(e,extend) RESULT(norm)
!              IMPORT dp, quadpackbase
!              REAL(dp), INTENT(IN) :: e(:)
!              CLASS(quadpackbase), INTENT(INOUT) :: extend
!              REAL(dp) :: norm
!            END FUNCTION

!         which is then used to perform norm calculations on the vector of
!         integrands. It is unused for scalar function integration.
!         If the pointer is not associated then the packaged routine MAXNORM is
!         called. This uses MAXVAL(ABS(e)).
!         Default: enorm => NULL()
! A fuller description of all the components is included in the
! comments immediately preceding this type definition.
      TYPE, PUBLIC :: quadpackbase
  ! Define Gauss rule to be used on each subinterval
        INTEGER(kind=i4) :: key = 6
  ! Maximum number of subdivisions allowed
        INTEGER(kind=i4) :: limit = 50
  ! Absolute and relative error tolerances required
        REAL (dp) :: epsabs = 0.0E0_dp
        REAL (dp) :: epsrel = sqrt(epsilon(1.0E0_dp))
  ! Total number of function evaluations and subintervals used
        INTEGER(kind=i4) :: nfunction_evaluations
        INTEGER(kind=i4) :: nvector_evaluations
  ! Estimate of absolute error in result
        REAL (dp) :: abserr
  ! Error status flag
        INTEGER(kind=i4) :: ier
  ! Number of components in vector integration 
        INTEGER(kind=i4) :: ncomponents
  ! Pointer to user-supplied function for norm calculations
        PROCEDURE(errnorm), POINTER, NOPASS :: enorm=>NULL()
      END TYPE quadpackbase

      ABSTRACT INTERFACE
        RECURSIVE FUNCTION fs(xval,extend) RESULT (yval)
          IMPORT dp, quadpackbase
          REAL (dp), INTENT (IN) :: xval(:)
          CLASS (quadpackbase), INTENT (INOUT) :: extend
          REAL (dp) :: yval(size(xval))
        END FUNCTION fs

        RECURSIVE FUNCTION fv(xval,extend) RESULT (yval)
          IMPORT dp, quadpackbase
          REAL (dp), INTENT (IN) :: xval(:)
          CLASS (quadpackbase), INTENT (INOUT) :: extend
          REAL (dp) :: yval(size(xval),extend%ncomponents)
        END FUNCTION fv

        FUNCTION errnorm(e,extend) RESULT (norm)
          IMPORT dp, quadpackbase
          REAL (dp), INTENT (IN) :: e(:)
          CLASS (quadpackbase), INTENT (INOUT) :: extend
          REAL (dp) :: norm
        END FUNCTION errnorm
    END INTERFACE

! This is a generic interface definition.  Either scalar or vector
! integration can use the same routine name, QAG2003.  The routine called
! is either QAG2003S or QAG2003V.  This choice depends on the type and 
! rank of intApprox. If intApprox is scalar, QAG2003S is called.
! If intApprox(:) is rank 1, QAG2003V is called.

      INTERFACE qag2003
        MODULE PROCEDURE qag2003s, qag2003v
      END INTERFACE

      PUBLIC qag2003, qag2003s, qag2003v, fs, fv, errnorm, dp

! Named constants
      REAL (dp), PARAMETER, PUBLIC :: half = 0.5E0_dp, zero = 0.0E0_dp, &
        one = 1.0E0_dp, two = 2.0E0_dp, onept5 = 1.5_dp, twohun = 200.0_dp, &
        fifty = 50.0E0_dp, cc = 200.0E0_dp, p99 = 0.99_dp, tenm5 = 0.00001_dp, &
        ten4 = 10000.0E0_dp, ten3 = 1000.0E0_dp

! Named constants for linking key values with rules
      INTEGER(kind=i4), PARAMETER :: gk15rule = 1, gk21rule = 2, gk31rule = 3, &
        gk41rule = 4, gk51rule = 5, gk61rule = 6, maxrules = 6

! Array defining number of points associated with each rule
      INTEGER(kind=i4), PARAMETER :: xrule(1:maxrules) = (/ 15, 21, 31, 41, 51, 61 /)

! Indexes for pointers into the one dimension arrays
! to extract each rule
      INTEGER(kind=i4), PARAMETER :: gkindex(1:maxrules+1) = (/ 1, 9, 20, 36, 57, 83, &
        114 /)
      INTEGER(kind=i4), PARAMETER :: gaussindex(1:maxrules+1) = (/ 1, 5, 10, 18, 28, &
        41, 56 /)

! Data for 15-point Gauss-Kronrod rule
      REAL (dp), PARAMETER :: xwg15(1:4) = (/ &
        0.129484966168869693270611432679082E0_dp, &
        0.279705391489276667901467771423780E0_dp, &
        0.381830050505118944950369775488975E0_dp, &
        0.417959183673469387755102040816327E0_dp/)

      REAL (dp), PARAMETER :: xxgk15(1:8) = (/ &
        0.949107912342758524526189684047851E0_dp, &
        0.741531185599394439863864773280788E0_dp, &
        0.405845151377397166906606412076961E0_dp, &
        0.991455371120812639206854697526329E0_dp, &
        0.864864423359769072789712788640926E0_dp, &
        0.586087235467691130294144838258730E0_dp, &
        0.207784955007898467600689403773245E0_dp, &
        0.000000000000000000000000000000000E0_dp/)

      REAL (dp), PARAMETER :: xwgk15(1:8) = (/ &
        0.063092092629978553290700663189204E0_dp, &
        0.140653259715525918745189590510238E0_dp, &
        0.190350578064785409913256402421014E0_dp, &
        0.022935322010529224963732008058970E0_dp, &
        0.104790010322250183839876322541518E0_dp, &
        0.169004726639267902826583426598550E0_dp, &
        0.204432940075298892414161999234649E0_dp, &
        0.209482141084727828012999174891714E0_dp/)

! Data for 21-point Gauss-Kronrod rule
      REAL (dp), PARAMETER :: xwg21(1:5) = (/ &
        0.066671344308688137593568809893332E0_dp, &
        0.149451349150580593145776339657697E0_dp, &
        0.219086362515982043995534934228163E0_dp, &
        0.269266719309996355091226921569469E0_dp, &
        0.295524224714752870173892994651338E0_dp/)

      REAL (dp), PARAMETER :: xxgk21(1:11) = (/ &
        0.973906528517171720077964012084452E0_dp, &
        0.865063366688984510732096688423493E0_dp, &
        0.679409568299024406234327365114874E0_dp, &
        0.433395394129247190799265943165784E0_dp, &
        0.148874338981631210884826001129720E0_dp, &
        0.995657163025808080735527280689003E0_dp, &
        0.930157491355708226001207180059508E0_dp, &
        0.780817726586416897063717578345042E0_dp, &
        0.562757134668604683339000099272694E0_dp, &
        0.294392862701460198131126603103866E0_dp, &
        0.000000000000000000000000000000000E0_dp/)

      REAL (dp), PARAMETER :: xwgk21(1:11) = (/ &
        0.032558162307964727478818972459390E0_dp, &
        0.075039674810919952767043140916190E0_dp, &
        0.109387158802297641899210590325805E0_dp, &
        0.134709217311473325928054001771707E0_dp, &
        0.147739104901338491374841515972068E0_dp, &
        0.011694638867371874278064396062192E0_dp, &
        0.054755896574351996031381300244580E0_dp, &
        0.093125454583697605535065465083366E0_dp, &
        0.123491976262065851077958109831074E0_dp, &
        0.142775938577060080797094273138717E0_dp, &
        0.149445554002916905664936468389821E0_dp/)

! Data for 31-point Gauss-Kronrod rule
      REAL (dp), PARAMETER :: xwg31(1:8) = (/ &
        0.030753241996117268354628393577204E0_dp, &
        0.070366047488108124709267416450667E0_dp, &
        0.107159220467171935011869546685869E0_dp, &
        0.139570677926154314447804794511028E0_dp, &
        0.166269205816993933553200860481209E0_dp, &
        0.186161000015562211026800561866423E0_dp, &
        0.198431485327111576456118326443839E0_dp, &
        0.202578241925561272880620199967519E0_dp/)

      REAL (dp), PARAMETER :: xxgk31(1:16) = (/ &
        0.987992518020485428489565718586613E0_dp, &
        0.937273392400705904307758947710209E0_dp, &
        0.848206583410427216200648320774217E0_dp, &
        0.724417731360170047416186054613938E0_dp, &
        0.570972172608538847537226737253911E0_dp, &
        0.394151347077563369897207370981045E0_dp, &
        0.201194093997434522300628303394596E0_dp, &
        0.998002298693397060285172840152271E0_dp, &
        0.967739075679139134257347978784337E0_dp, &
        0.897264532344081900882509656454496E0_dp, &
        0.790418501442465932967649294817947E0_dp, &
        0.650996741297416970533735895313275E0_dp, &
        0.485081863640239680693655740232351E0_dp, &
        0.299180007153168812166780024266389E0_dp, &
        0.101142066918717499027074231447392E0_dp, &
        0.000000000000000000000000000000000E0_dp/)

      REAL (dp), PARAMETER :: xwgk31(1:16) = (/ &
        0.015007947329316122538374763075807E0_dp, &
        0.035346360791375846222037948478360E0_dp, &
        0.053481524690928087265343147239430E0_dp, &
        0.069854121318728258709520077099147E0_dp, &
        0.083080502823133021038289247286104E0_dp, &
        0.093126598170825321225486872747346E0_dp, &
        0.099173598721791959332393173484603E0_dp, &
        0.005377479872923348987792051430128E0_dp, &
        0.025460847326715320186874001019653E0_dp, &
        0.044589751324764876608227299373280E0_dp, &
        0.062009567800670640285139230960803E0_dp, &
        0.076849680757720378894432777482659E0_dp, &
        0.088564443056211770647275443693774E0_dp, &
        0.096642726983623678505179907627589E0_dp, &
        0.100769845523875595044946662617570E0_dp, &
        0.101330007014791549017374792767493E0_dp/)

! Data for 41-point Gauss-Kronrod rule
      REAL (dp), PARAMETER :: xwg41(1:10) = (/ &
        0.017614007139152118311861962351853E0_dp, &
        0.040601429800386941331039952274932E0_dp, &
        0.062672048334109063569506535187042E0_dp, &
        0.083276741576704748724758143222046E0_dp, &
        0.101930119817240435036750135480350E0_dp, &
        0.118194531961518417312377377711382E0_dp, &
        0.131688638449176626898494499748163E0_dp, &
        0.142096109318382051329298325067165E0_dp, &
        0.149172986472603746787828737001969E0_dp, &
        0.152753387130725850698084331955098E0_dp/)

      REAL (dp), PARAMETER :: xxgk41(1:21) = (/ &
        0.993128599185094924786122388471320E0_dp, &
        0.963971927277913791267666131197277E0_dp, &
        0.912234428251325905867752441203298E0_dp, &
        0.839116971822218823394529061701521E0_dp, &
        0.746331906460150792614305070355642E0_dp, &
        0.636053680726515025452836696226286E0_dp, &
        0.510867001950827098004364050955251E0_dp, &
        0.373706088715419560672548177024927E0_dp, &
        0.227785851141645078080496195368575E0_dp, &
        0.076526521133497333754640409398838E0_dp, &
        0.998859031588277663838315576545863E0_dp, &
        0.981507877450250259193342994720217E0_dp, &
        0.940822633831754753519982722212443E0_dp, &
        0.878276811252281976077442995113078E0_dp, &
        0.795041428837551198350638833272788E0_dp, &
        0.693237656334751384805490711845932E0_dp, &
        0.575140446819710315342946036586425E0_dp, &
        0.443593175238725103199992213492640E0_dp, &
        0.301627868114913004320555356858592E0_dp, &
        0.152605465240922675505220241022678E0_dp, &
        0.000000000000000000000000000000000E0_dp/)

      REAL (dp), PARAMETER :: xwgk41(1:21) = (/ &
        0.008600269855642942198661787950102E0_dp, &
        0.020388373461266523598010231432755E0_dp, &
        0.031287306777032798958543119323801E0_dp, &
        0.041668873327973686263788305936895E0_dp, &
        0.050944573923728691932707670050345E0_dp, &
        0.059111400880639572374967220648594E0_dp, &
        0.065834597133618422111563556969398E0_dp, &
        0.071054423553444068305790361723210E0_dp, &
        0.074582875400499188986581418362488E0_dp, &
        0.076377867672080736705502835038061E0_dp, &
        0.003073583718520531501218293246031E0_dp, &
        0.014626169256971252983787960308868E0_dp, &
        0.025882133604951158834505067096153E0_dp, &
        0.036600169758200798030557240707211E0_dp, &
        0.046434821867497674720231880926108E0_dp, &
        0.055195105348285994744832372419777E0_dp, &
        0.062653237554781168025870122174255E0_dp, &
        0.068648672928521619345623411885368E0_dp, &
        0.073030690332786667495189417658913E0_dp, &
        0.075704497684556674659542775376617E0_dp, &
        0.076600711917999656445049901530102E0_dp/)

! Data for 51-point Gauss-Kronrod rule
      REAL (dp), PARAMETER :: xwg51(1:13) = (/ &
        0.011393798501026287947902964113235E0_dp, &
        0.026354986615032137261901815295299E0_dp, &
        0.040939156701306312655623487711646E0_dp, &
        0.054904695975835191925936891540473E0_dp, &
        0.068038333812356917207187185656708E0_dp, &
        0.080140700335001018013234959669111E0_dp, &
        0.091028261982963649811497220702892E0_dp, &
        0.100535949067050644202206890392686E0_dp, &
        0.108519624474263653116093957050117E0_dp, &
        0.114858259145711648339325545869556E0_dp, &
        0.119455763535784772228178126512901E0_dp, &
        0.122242442990310041688959518945852E0_dp, &
        0.123176053726715451203902873079050E0_dp/)

      REAL (dp), PARAMETER :: xxgk51(1:26) = (/ &
        0.995556969790498097908784946893902E0_dp, &
        0.976663921459517511498315386479594E0_dp, &
        0.942974571228974339414011169658471E0_dp, &
        0.894991997878275368851042006782805E0_dp, &
        0.833442628760834001421021108693570E0_dp, &
        0.759259263037357630577282865204361E0_dp, &
        0.673566368473468364485120633247622E0_dp, &
        0.577662930241222967723689841612654E0_dp, &
        0.473002731445714960522182115009192E0_dp, &
        0.361172305809387837735821730127641E0_dp, &
        0.243866883720988432045190362797452E0_dp, &
        0.122864692610710396387359818808037E0_dp, &
        0.999262104992609834193457486540341E0_dp, &
        0.988035794534077247637331014577406E0_dp, &
        0.961614986425842512418130033660167E0_dp, &
        0.920747115281701561746346084546331E0_dp, &
        0.865847065293275595448996969588340E0_dp, &
        0.797873797998500059410410904994307E0_dp, &
        0.717766406813084388186654079773298E0_dp, &
        0.626810099010317412788122681624518E0_dp, &
        0.526325284334719182599623778158010E0_dp, &
        0.417885382193037748851814394594572E0_dp, &
        0.303089538931107830167478909980339E0_dp, &
        0.183718939421048892015969888759528E0_dp, &
        0.061544483005685078886546392366797E0_dp, &
        0.000000000000000000000000000000000E0_dp/)

      REAL (dp), PARAMETER :: xwgk51(1:26) = (/ &
        0.005561932135356713758040236901066E0_dp, &
        0.013236229195571674813656405846976E0_dp, &
        0.020435371145882835456568292235939E0_dp, &
        0.027475317587851737802948455517811E0_dp, &
        0.034002130274329337836748795229551E0_dp, &
        0.040083825504032382074839284467076E0_dp, &
        0.045502913049921788909870584752660E0_dp, &
        0.050277679080715671963325259433440E0_dp, &
        0.054251129888545490144543370459876E0_dp, &
        0.057437116361567832853582693939506E0_dp, &
        0.059720340324174059979099291932562E0_dp, &
        0.061128509717053048305859030416293E0_dp, &
        0.001987383892330315926507851882843E0_dp, &
        0.009473973386174151607207710523655E0_dp, &
        0.016847817709128298231516667536336E0_dp, &
        0.024009945606953216220092489164881E0_dp, &
        0.030792300167387488891109020215229E0_dp, &
        0.037116271483415543560330625367620E0_dp, &
        0.042872845020170049476895792439495E0_dp, &
        0.047982537138836713906392255756915E0_dp, &
        0.052362885806407475864366712137873E0_dp, &
        0.055950811220412317308240686382747E0_dp, &
        0.058689680022394207961974175856788E0_dp, &
        0.060539455376045862945360267517565E0_dp, &
        0.061471189871425316661544131965264E0_dp, &
        0.061580818067832935078759824240066E0_dp/)

! Data for 61-point Gauss-Kronrod rule
      REAL (dp), PARAMETER :: xwg61(1:15) = (/ &
        0.007968192496166605615465883474674E0_dp, &
        0.018466468311090959142302131912047E0_dp, &
        0.028784707883323369349719179611292E0_dp, &
        0.038799192569627049596801936446348E0_dp, &
        0.048402672830594052902938140422808E0_dp, &
        0.057493156217619066481721689402056E0_dp, &
        0.065974229882180495128128515115962E0_dp, &
        0.073755974737705206268243850022191E0_dp, &
        0.080755895229420215354694938460530E0_dp, &
        0.086899787201082979802387530715126E0_dp, &
        0.092122522237786128717632707087619E0_dp, &
        0.096368737174644259639468626351810E0_dp, &
        0.099593420586795267062780282103569E0_dp, &
        0.101762389748405504596428952168554E0_dp, &
        0.102852652893558840341285636705415E0_dp/)

      REAL (dp), PARAMETER :: xxgk61(1:31) = (/ &
        0.996893484074649540271630050918695E0_dp, &
        0.983668123279747209970032581605663E0_dp, &
        0.960021864968307512216871025581798E0_dp, &
        0.926200047429274325879324277080474E0_dp, &
        0.882560535792052681543116462530226E0_dp, &
        0.829565762382768397442898119732502E0_dp, &
        0.767777432104826194917977340974503E0_dp, &
        0.697850494793315796932292388026640E0_dp, &
        0.620526182989242861140477556431189E0_dp, &
        0.536624148142019899264169793311073E0_dp, &
        0.447033769538089176780609900322854E0_dp, &
        0.352704725530878113471037207089374E0_dp, &
        0.254636926167889846439805129817805E0_dp, &
        0.153869913608583546963794672743256E0_dp, &
        0.051471842555317695833025213166723E0_dp, &
        0.999484410050490637571325895705811E0_dp, &
        0.991630996870404594858628366109486E0_dp, &
        0.973116322501126268374693868423707E0_dp, &
        0.944374444748559979415831324037439E0_dp, &
        0.905573307699907798546522558925958E0_dp, &
        0.857205233546061098958658510658944E0_dp, &
        0.799727835821839083013668942322683E0_dp, &
        0.733790062453226804726171131369528E0_dp, &
        0.660061064126626961370053668149271E0_dp, &
        0.579345235826361691756024932172540E0_dp, &
        0.492480467861778574993693061207709E0_dp, &
        0.400401254830394392535476211542661E0_dp, &
        0.304073202273625077372677107199257E0_dp, &
        0.204525116682309891438957671002025E0_dp, &
        0.102806937966737030147096751318001E0_dp, &
        0.000000000000000000000000000000000E0_dp/)

      REAL (dp), PARAMETER :: xwgk61(1:31) = (/ &
        0.003890461127099884051267201844516E0_dp, &
        0.009273279659517763428441146892024E0_dp, &
        0.014369729507045804812451432443580E0_dp, &
        0.019414141193942381173408951050128E0_dp, &
        0.024191162078080601365686370725232E0_dp, &
        0.028754048765041292843978785354334E0_dp, &
        0.032981447057483726031814191016854E0_dp, &
        0.036882364651821229223911065617136E0_dp, &
        0.040374538951535959111995279752468E0_dp, &
        0.043452539701356069316831728117073E0_dp, &
        0.046059238271006988116271735559374E0_dp, &
        0.048185861757087129140779492298305E0_dp, &
        0.049795683427074206357811569379942E0_dp, &
        0.050881795898749606492297473049805E0_dp, &
        0.051426128537459025933862879215781E0_dp, &
        0.001389013698677007624551591226760E0_dp, &
        0.006630703915931292173319826369750E0_dp, &
        0.011823015253496341742232898853251E0_dp, &
        0.016920889189053272627572289420322E0_dp, &
        0.021828035821609192297167485738339E0_dp, &
        0.026509954882333101610601709335075E0_dp, &
        0.030907257562387762472884252943092E0_dp, &
        0.034979338028060024137499670731468E0_dp, &
        0.038678945624727592950348651532281E0_dp, &
        0.041969810215164246147147541285970E0_dp, &
        0.044814800133162663192355551616723E0_dp, &
        0.047185546569299153945261478181099E0_dp, &
        0.049055434555029778887528165367238E0_dp, &
        0.050405921402782346840893085653585E0_dp, &
        0.051221547849258772170656282604944E0_dp, &
        0.051494729429451567558340433647099E0_dp/)

      REAL (dp), TARGET :: xwg(1:gaussindex(maxrules+1)-1) = (/ xwg15, xwg21, &
        xwg31, xwg41, xwg51, xwg61/)
      REAL (dp), TARGET :: xwgk(1:gkindex(maxrules+1)-1) = (/ xwgk15, xwgk21, &
        xwgk31, xwgk41, xwgk51, xwgk61/)
      REAL (dp), TARGET :: xxgk(1:gkindex(maxrules+1)-1) = (/ xxgk15, xxgk21, &
        xxgk31, xxgk41, xxgk51, xxgk61/)

    CONTAINS

      RECURSIVE SUBROUTINE qag2003v(f,extype,a,b,intapprox)

!*****************************************************************************80

!! QAG2003V approximates a vector integral over a finite interval.

!  Discussion:

!    The routine calculates an approximation intApprox to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - intApprox || <= max ( EPSABS, EPSREL * ||I|| ).

!    QAG2003 is a simple globally adaptive integrator using the strategy of 
!    Aind (Piessens, 1973).  It is possible to choose between 6 pairs of
!    Gauss-Kronrod quadrature formulae for the rule evaluation component. 
!    The pairs of high degree of precision are suitable for handling
!    integration difficulties due to a strongly oscillating integrand.

!  Authors:

!    Fortran 77 version:
!      Robert Piessens, Elise de Doncker-Kapenger, 
!      Christian Ueberhuber, David Kahaner
!    Fortran 2003 version:
!      Richard Hanson, Tim Hopkins

!  Reference:

!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983

!  Parameters:

!   External REAL(dp), FUNCTION f 
!     User provided function, of the form

!          [RECURSIVE] FUNCTION f ( x, extype) RESULT(funVals)
!          REAL, INTENT(IN) :: x(:), f(size(x))
!          CLASS(quadpackBase) :: extype
!          REAL(dp) :: funVals(SIZE(x))

!    which evaluates the integrand function at the vector of points given by x.

!    Input/Output: CLASS(quadpackBase) :: extype
!      Used to set the required accuracy tolerances, order of rule to be used,
!      limit on the number of subdivisions, etc.
!      See description of components in definition of type quadpackBase at the
!      head of this file.

!    Input: REAL(dp) :: a, b,
!      define the lower and upper limits of the integration respectively.

!    Output: REAL(dp) :: intApprox
!      returns the computed approximations to the vector of integrals.

! Type and define dummy function values -
  PROCEDURE(fv) :: f

        CLASS (quadpackbase), INTENT (INOUT) :: extype
        REAL (dp), INTENT (IN) :: a, b
        REAL (dp), INTENT (OUT) :: intapprox(:)
        INTEGER(kind=i4) last

  PROCEDURE(errnorm), POINTER  :: norm

! Used only with vector integrals.  If the pointer is not
! associated - i.e. still points to NULL() - then a packaged
! routine MAXNORM(E(:),EXTEND) is called.  Otherwise
! a user can write another norm with these same arguments.
! They would then point to this new norm in this component.    
        IF ( .NOT. associated(extype%enorm)) THEN
! For vector integrals use the max abs norm.  
          norm => maxnorm
        ELSE
! Otherwise the calling program has defined a norm.  
          norm => extype%enorm
        END IF

! Needed to define assumed size dimension of user-written
! evaluation code.  
        extype%ncomponents = size(intapprox)
        CALL qage2003v(f,extype,a,b,intapprox,last,norm)
      END SUBROUTINE qag2003v

      RECURSIVE SUBROUTINE qage2003v(f,extype,a,b,intapprox,last,norm)
        IMPLICIT NONE
!*****************************************************************************80

!! QAGE2003V estimates a vector-valued definite integral.

!  Discussion:

!    The routine calculates an approximation intApprox to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - intApprox || <= max ( EPSABS, EPSREL * ||I|| ).

!  Author:

!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner

!  Reference:

!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983

!  Parameters:

!    Input, external real F, the name of the function routine, of the form
!      [recursive] function f ( x, EXTYPE)
!      real x(:), f(size(x))
!    which evaluates the integrand function.

!    Input, real A, B, the limits of integration.

!    Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.

!    Input, integer KEY, chooses the order of the local integration rule:
!    1,  7 Gauss points, 15 Gauss-Kronrod points,
!    2, 10 Gauss points, 21 Gauss-Kronrod points,
!    3, 15 Gauss points, 31 Gauss-Kronrod points,
!    4, 20 Gauss points, 41 Gauss-Kronrod points,
!    5, 25 Gauss points, 51 Gauss-Kronrod points,
!    6, 30 Gauss points, 61 Gauss-Kronrod points.

!    Input, integer LIMIT, the maximum number of subintervals that
!    can be used.

!    Output, real intApprox, the estimated value of the integral.

!    Output, real ABSERR, an estimate of || I - intApprox ||.

!    Output, integer IER, return code.
!    0, normal and reliable termination of the routine.  It is assumed that the 
!      requested accuracy has been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can 
!      allow more subdivisions by increasing the value of LIMIT in QAG. 
!      However, if this yields no improvement it is advised to analyze the
!      integrand to determine the integration difficulties.  If the position
!      of a local difficulty can be determined, such as a singularity or
!      discontinuity within the interval) one will probably gain from 
!      splitting up the interval at this point and calling the integrator 
!      on the subranges.  If possible, an appropriate special-purpose 
!      integrator should be used which is designed for handling the type 
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.

!    Workspace, real ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
!    through LAST the left and right ends of the partition subintervals.

!    Workspace, real RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.

!    Workspace, real ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.

!    Output, integer IORD(LIMIT), the first K elements of which are pointers 
!    to the error estimates over the subintervals, such that
!    elist(iord(1)), ..., elist(iord(k)) form a decreasing sequence, with
!    k = last if last <= (limit/2+2), and k = limit+1-last otherwise.

!    Output, integer LAST, the number of subintervals actually produced 
!    in the subdivision process.

!  Local parameters:

!    alist     - list of left end points of all subintervals
!                       considered up to now
!    blist     - list of right end points of all subintervals
!                       considered up to now
!    elist(i)  - error estimate applying to rlist(i)
!    maxerr    - pointer to the interval with largest error estimate
!    errmax    - elist(maxerr)
!    area      - sum of the integrals over the subintervals
!    errsum    - sum of the errors over the subintervals
!    errbnd    - requested accuracy max(epsabs,epsrel*abs(intApprox))
!    defab1    - variable for the left subinterval
!    defab2    - variable for the right subinterval
!    last      - index for subdivision


! Type and define dummy function values -
  PROCEDURE(fv) :: f
        CLASS (quadpackbase), INTENT (INOUT) :: extype
        REAL (dp), INTENT (IN) :: a, b
        REAL (dp), INTENT (INOUT) :: intapprox(:)
        INTEGER(kind=i4), INTENT (INOUT) :: last
  PROCEDURE(errnorm), POINTER :: norm

        REAL (dp) a12, rl
        REAL (dp) a1, a2, b1, b2, c
        REAL (dp) errbnd, errmax, errsum
        REAL (dp) error1, error2, erro12
        REAL (dp), ALLOCATABLE :: alist(:), blist(:), area(:), area1(:), &
          area2(:), area12(:), resabs(:), elist(:), rlist(:,:), defabs(:), &
          defab1(:), defab2(:)
        INTEGER(kind=i4) iroff1, iroff2, keyf, maxerr
        INTEGER(kind=i4) nrmax, i
        INTEGER(kind=i4), ALLOCATABLE :: iord(:)

        ASSOCIATE (limit => extype%limit,&
                   epsabs=> extype%epsabs,&
                   epsrel=> extype%epsrel,&
                   abserr=> extype%abserr,&
                   key   => extype%key,&
                   ier   => extype%ier,&
                   nvec  => extype%nvector_evaluations, &
                   neval => extype%nfunction_evaluations, &
                     M   => size(intApprox))
! Allocate space for lists -- use limit to set sizes
        ALLOCATE (alist(limit),blist(limit),elist(limit),rlist(limit,m), &
          area(m),area12(m),area1(m),area2(m),resabs(m),defabs(m),defab1(m), &
          defab2(m),iord(limit))

!  Test on validity of parameters.

        ier = 0
        last = 0
        intapprox = zero
        abserr = zero
        alist(1) = a
        blist(1) = b
        rlist(1,:) = zero
        elist(1) = zero
        iord(1) = 0
        nvec = 0
        neval = 0

        IF (epsabs<zero .AND. epsrel<zero) THEN
          ier = 6
          RETURN
        END IF

!  First approximation to the integral.


        keyf = min(max(key,1),6)
        c = keyf

        CALL qkgenv(f,extype,a,b,intapprox,defabs,resabs,norm)
        nvec = nvec + 1
        last = 1
        rlist(1,:) = intapprox
        elist(1) = abserr
        iord(1) = 1

!  Test on accuracy.

        errbnd = max(epsabs,epsrel*norm(abs(intapprox),extype))

        IF (abserr<=half*epsilon(abserr)*norm(defabs,extype) .AND. &
            errbnd<abserr) THEN
          ier = 2
        END IF

        IF (limit==1) THEN
          ier = 1
        END IF

        IF (ier/=0 .OR. (abserr<=errbnd .AND. abserr/=norm(resabs, &
            extype)) .OR. abserr==zero) THEN

          neval = nvec*xrule(keyf)
          RETURN

        END IF

!  Initialization.

        errmax = abserr
        maxerr = 1
        area = intapprox
        errsum = abserr
        nrmax = 1
        iroff1 = 0
        iroff2 = 0

        DO last = 2, limit

!  Bisect the subinterval with the largest error estimate.

          a1 = alist(maxerr)
          b1 = half*(alist(maxerr)+blist(maxerr))
          a2 = b1
          b2 = blist(maxerr)

          CALL qkgenv(f,extype,a1,b1,area1,resabs,defab1,norm)
          error1 = abserr
          CALL qkgenv(f,extype,a2,b2,area2,resabs,defab2,norm)
          error2 = abserr
          nvec = nvec + 2

!  Improve previous approximations to integral and error and
!  test for accuracy.

          area12 = area1 + area2
          erro12 = error1 + error2
          errsum = errsum + erro12 - errmax
          area = area + area12 - rlist(maxerr,:)

          IF (norm(defab1,extype)/=error1 .AND. norm(defab2,extype)/=error2) &
              THEN
            a12 = norm(abs(area12),extype)
            rl = norm(abs(rlist(maxerr,:)-area12),extype)
            IF (rl<=tenm5*a12 .AND. p99*errmax<=erro12) THEN
              iroff1 = iroff1 + 1
            END IF

            IF (10<last .AND. errmax<erro12) THEN
              iroff2 = iroff2 + 1
            END IF

          END IF

          rlist(maxerr,:) = area1
          rlist(last,:) = area2
          errbnd = max(epsabs,epsrel*norm(abs(area),extype))

!  Test for roundoff error and eventually set error flag.

          IF (errbnd<errsum) THEN

            IF (6<=iroff1 .OR. 20<=iroff2) THEN
              ier = 2
            END IF

!  Set error flag in the case that the number of subintervals
!  equals limit.

            IF (last==limit) THEN
              ier = 1
            END IF

!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.

            IF (max(abs(a1),abs(b2))<=(one+c*ten3*epsilon(a1))*(abs( &
                a2)+ten4*tiny(a2))) THEN
              ier = 3
            END IF

          END IF

!  Append the newly-created intervals to the list.

          IF (error2<=error1) THEN
            alist(last) = a2
            blist(maxerr) = b1
            blist(last) = b2
            elist(maxerr) = error1
            elist(last) = error2
          ELSE
            alist(maxerr) = a2
            alist(last) = a1
            blist(last) = b1
            rlist(maxerr,:) = area2
            rlist(last,:) = area1
            elist(maxerr) = error2
            elist(last) = error1
          END IF

!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with the largest error estimate (to be bisected next).

          CALL qsort_2003(limit,last,maxerr,errmax,elist,iord,nrmax)
          IF (ier/=0 .OR. errsum<=errbnd) THEN
            EXIT
          END IF

        END DO

!  Compute final result.

        intapprox = zero
        DO i = 1, last
          intapprox = intapprox + rlist(i,:)
        END DO
        abserr = errsum
        neval = nvec*xrule(keyf)

END ASSOCIATE

      END SUBROUTINE qage2003v


      RECURSIVE SUBROUTINE qag2003s(f,extype,a,b,intapprox)

!*****************************************************************************80

!! QAG2003S approximates a scalar integral over a finite interval.

!  Discussion:

!    The routine calculates an approximation intApprox to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - intApprox || <= max ( EPSABS, EPSREL * ||I|| ).

!    QAG2003 is a simple globally adaptive integrator using the strategy of 
!    Aind (Piessens, 1973).  It is possible to choose between 6 pairs of
!    Gauss-Kronrod quadrature formulae for the rule evaluation component. 
!    The pairs of high degree of precision are suitable for handling
!    integration difficulties due to a strongly oscillating integrand.

!  Author:

!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner

!  Reference:

!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!    Modification using Fortran 2003 by R. J. Hanson
!  Parameters:

!    Input, external F, the name of the function routine, of the form
!      [recursive] function f ( x, EXTYPE)
!      real x(:), f(size(x))
!    which evaluates the integrand function.

!    Input, real A, B, the limits of integration.

!    Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.

!    Input, integer KEY, chooses the order of the local integration rule:
!    1,  7 Gauss points, 15 Gauss-Kronrod points,
!    2, 10 Gauss points, 21 Gauss-Kronrod points,
!    3, 15 Gauss points, 31 Gauss-Kronrod points,
!    4, 20 Gauss points, 41 Gauss-Kronrod points,
!    5, 25 Gauss points, 51 Gauss-Kronrod points,
!    6, 30 Gauss points, 61 Gauss-Kronrod points.

!    Output, real intApprox, the estimated value of the integral.

!    Output, real ABSERR, an estimate of || I - intApprox ||.

!    Output, integer IER, return code.
!    0, normal and reliable termination of the routine.  It is assumed that the 
!      requested accuracy has been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can 
!      allow more subdivisions by increasing the value of LIMIT in QAG. 
!      However, if this yields no improvement it is advised to analyze the
!      integrand to determine the integration difficulties.  If the position
!      of a local difficulty can be determined, such as a singularity or
!      discontinuity within the interval) one will probably gain from 
!      splitting up the interval at this point and calling the integrator 
!      on the subranges.  If possible, an appropriate special-purpose 
!      integrator should be used which is designed for handling the type 
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.

!  Local parameters:

!    LIMIT is the maximum number of subintervals allowed in
!    the subdivision process of QAGE.

! Type and define dummy function values -
  PROCEDURE(fs) :: f

        CLASS (quadpackbase), INTENT (INOUT) :: extype
        REAL (dp), INTENT (IN) :: a, b
        REAL (dp), INTENT (OUT) :: intapprox

        INTEGER(kind=i4) :: last

        CALL qage2003s(f,extype,a,b,intapprox,last)

      END SUBROUTINE qag2003s

      RECURSIVE SUBROUTINE qage2003s(f,extype,a,b,intapprox,last)

!*****************************************************************************80

!! QAGE2003 estimates a definite integral.

!  Discussion:

!    The routine calculates an approximation intApprox to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - intApprox || <= max ( EPSABS, EPSREL * ||I|| ).

!  Author:

!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner

!  Reference:

!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983

!  Parameters:

!    Input, external real F, the name of the function routine, of the form
!      [recursive] function f ( x, EXTYPE)
!      real x(:), f(size(x))
!    which evaluates the integrand function.

!    Input, real A, B, the limits of integration.

!    Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.

!    Input, integer KEY, chooses the order of the local integration rule:
!    1,  7 Gauss points, 15 Gauss-Kronrod points,
!    2, 10 Gauss points, 21 Gauss-Kronrod points,
!    3, 15 Gauss points, 31 Gauss-Kronrod points,
!    4, 20 Gauss points, 41 Gauss-Kronrod points,
!    5, 25 Gauss points, 51 Gauss-Kronrod points,
!    6, 30 Gauss points, 61 Gauss-Kronrod points.

!    Input, integer LIMIT, the maximum number of subintervals that
!    can be used.

!    Output, real intApprox, the estimated value of the integral.

!    Output, real ABSERR, an estimate of || I - intApprox ||.

!    Output, integer IER, return code.
!    0, normal and reliable termination of the routine.  It is assumed that the 
!      requested accuracy has been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can 
!      allow more subdivisions by increasing the value of LIMIT in QAG. 
!      However, if this yields no improvement it is advised to analyze the
!      integrand to determine the integration difficulties.  If the position
!      of a local difficulty can be determined, such as a singularity or
!      discontinuity within the interval) one will probably gain from 
!      splitting up the interval at this point and calling the integrator 
!      on the subranges.  If possible, an appropriate special-purpose 
!      integrator should be used which is designed for handling the type 
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.

!    Workspace, real ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
!    through LAST the left and right ends of the partition subintervals.

!    Workspace, real RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.

!    Workspace, real ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.

!    Output, integer IORD(LIMIT), the first K elements of which are pointers 
!    to the error estimates over the subintervals, such that
!    elist(iord(1)), ..., elist(iord(k)) form a decreasing sequence, with
!    k = last if last <= (limit/2+2), and k = limit+1-last otherwise.

!    Output, integer LAST, the number of subintervals actually produced 
!    in the subdivision process.

!  Local parameters:

!    alist     - list of left end points of all subintervals
!                       considered up to now
!    blist     - list of right end points of all subintervals
!                       considered up to now
!    elist(i)  - error estimate applying to rlist(i)
!    maxerr    - pointer to the interval with largest error estimate
!    errmax    - elist(maxerr)
!    area      - sum of the integrals over the subintervals
!    errsum    - sum of the errors over the subintervals
!    errbnd    - requested accuracy max(epsabs,epsrel*abs(intApprox))
!    defab1    - variable for the left subinterval
!    defab2    - variable for the right subinterval
!    last      - index for subdivision


! Type and define dummy function values -
  PROCEDURE(FS) :: F
        CLASS (quadpackbase), INTENT (INOUT) :: extype
        REAL (dp), INTENT (IN) :: a, b
        REAL (dp), INTENT (OUT) :: intapprox
        INTEGER(kind=i4), INTENT (OUT) :: last

        REAL (dp) area, area1, area12, area2
        REAL (dp) a1, a2, b1, b2, c, defabs, defab1, defab2
        REAL (dp) errbnd, errmax, error1
        REAL (dp) error2, erro12, errsum, resabs
        REAL (dp), ALLOCATABLE :: alist(:), blist(:), elist(:), rlist(:)
        INTEGER(kind=i4) iroff1, iroff2, keyf, maxerr
        INTEGER(kind=i4) nrmax
        INTEGER(kind=i4), ALLOCATABLE :: iord(:)

        ASSOCIATE (limit => extype%limit,&
                   epsabs=> extype%epsabs,&
                   epsrel=> extype%epsrel,&
                   abserr=> extype%abserr,&
                   key   => extype%key,&
                   nvec  => extype%nvector_evaluations, &
                   neval => extype%nfunction_evaluations, &
                   ier   => extype%ier)

! Allocate space for lists -- use limit to set sizes
        ALLOCATE (alist(limit),blist(limit),elist(limit),rlist(limit), &
          iord(limit))

!  Test on validity of parameters.

        ier = 0
        last = 0
        intapprox = zero
        abserr = zero
        alist(1) = a
        blist(1) = b
        rlist(1) = zero
        elist(1) = zero
        iord(1) = 0
        nvec = 0
        neval = 0

        IF (epsabs<zero .AND. epsrel<zero) THEN
          ier = 6
          RETURN
        END IF

!  First approximation to the integral.


!!!!!! NEED TO DO SOMETHING WITH KEYF HERE !!!!!
        keyf = min(max(key,1),6)
        c = keyf

        CALL qkgens(f,extype,a,b,intapprox,defabs,resabs)
        nvec = 1
        last = 1
        rlist(1) = intapprox
        elist(1) = abserr
        iord(1) = 1

!  Test on accuracy.

        errbnd = max(epsabs,epsrel*abs(intapprox))

        IF (abserr<=half*epsilon(defabs)*defabs .AND. errbnd<abserr) THEN
          ier = 2
        END IF

        IF (limit==1) THEN
          ier = 1
        END IF

        IF (ier/=0 .OR. (abserr<=errbnd .AND. abserr/=resabs) .OR. &
            abserr==zero) THEN

          neval = nvec*xrule(keyf)
          RETURN

        END IF

!  Initialization.

        errmax = abserr
        maxerr = 1
        area = intapprox
        errsum = abserr
        nrmax = 1
        iroff1 = 0
        iroff2 = 0

        DO last = 2, limit

!  Bisect the subinterval with the largest error estimate.

          a1 = alist(maxerr)
          b1 = half*(alist(maxerr)+blist(maxerr))
          a2 = b1
          b2 = blist(maxerr)


          CALL qkgens(f,extype,a1,b1,area1,resabs,defab1)
          error1 = abserr
          CALL qkgens(f,extype,a2,b2,area2,resabs,defab2)
          error2 = abserr
          nvec = nvec + 2

!  Improve previous approximations to integral and error and
!  test for accuracy.

          area12 = area1 + area2
          erro12 = error1 + error2
          errsum = errsum + erro12 - errmax
          area = area + area12 - rlist(maxerr)

          IF (defab1/=error1 .AND. defab2/=error2) THEN

            IF (abs(rlist(maxerr)-area12)<=tenm5*abs(area12) .AND. p99*errmax &
                <=erro12) THEN
              iroff1 = iroff1 + 1
            END IF

            IF (10<last .AND. errmax<erro12) THEN
              iroff2 = iroff2 + 1
            END IF

          END IF

          rlist(maxerr) = area1
          rlist(last) = area2
          errbnd = max(epsabs,epsrel*abs(area))

!  Test for roundoff error and eventually set error flag.

          IF (errbnd<errsum) THEN

            IF (6<=iroff1 .OR. 20<=iroff2) THEN
              ier = 2
            END IF

!  Set error flag in the case that the number of subintervals
!  equals limit.

            IF (last==limit) THEN
              ier = 1
            END IF

!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.

            IF (max(abs(a1),abs(b2))<=(one+c*ten3*epsilon(a1))*(abs( &
                a2)+ten4*tiny(a2))) THEN
              ier = 3
            END IF

          END IF

!  Append the newly-created intervals to the list.

          IF (error2<=error1) THEN
            alist(last) = a2
            blist(maxerr) = b1
            blist(last) = b2
            elist(maxerr) = error1
            elist(last) = error2
          ELSE
            alist(maxerr) = a2
            alist(last) = a1
            blist(last) = b1
            rlist(maxerr) = area2
            rlist(last) = area1
            elist(maxerr) = error2
            elist(last) = error1
          END IF

!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with the largest error estimate (to be bisected next).

          CALL qsort_2003(limit,last,maxerr,errmax,elist,iord,nrmax)
          IF (ier/=0 .OR. errsum<=errbnd) THEN
            EXIT
          END IF

        END DO

!  Compute final result.

        intapprox = sum(rlist(1:min(last,limit)))

        abserr = errsum
        neval = nvec*xrule(keyf)

END ASSOCIATE
      END SUBROUTINE qage2003s

      RECURSIVE SUBROUTINE qkgens(f,extype,a,b,intapprox,resabs,resasc)

        REAL (dp), INTENT (IN) :: a, b
        REAL (dp), INTENT (INOUT) :: intapprox, resabs, resasc
      PROCEDURE(fs) :: f
        CLASS (quadpackbase), INTENT (INOUT) :: extype

!            To compute I = Integral of F over (A,B), with error
!                           estimate
!                       J = Integral of ABS(F) over (A,B)
!           Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven

!           PARAMETERS
!            ON ENTRY
!              F      - Real
!                       Function subprogram defining the scalar integrand
!                       FUNCTION F(X,EXTYPE). The actual name for F needs to be
!                       Declared E X T E R N A L in the driver program or an
!                       equivalent via use association.  Note that the argument
!                       X(:) is an array so the result is an array of the same
!                       size.

!              EXTYPE   The extended or base data type of class
!                       QUADPACKBASE.  This allows a user program
!                       to pass data to the evaluation function, F.            

!              A      - Real
!                       Lower limit of integration

!              B      - Real
!                       Upper limit of integration

!            ON RETURN
!              intApprox - Real
!                       Approximation to the integral I
!                       intApprox is computed by applying a
!                       KRONROD RULE (RESK) obtained by optimal addition
!                       of abscissae to a GAUSS RULE (RESG).

!              ABSERR - Real
!                       Estimate of the modulus of the absolute error,
!                       which should not exceed ABS(I-intApprox)

!              RESABS - Real
!                       Approximation to the integral J

!              RESASC - Real
!                       Approximation to the integral of ABS(F-I/(B-A))
!                       over (A,B)





!           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
!           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
!           CORRESPONDING WEIGHTS ARE GIVEN.

!           XGK    - ABSCISSAE OF THE KRONROD RULE
!                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 
!                    GAUSS RULE
!                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!                    ADDED TO THE  GAUSS RULE

!           WGK    - WEIGHTS OF THE KRONROD RULE

!           WG     - WEIGHTS OF THE GAUSS RULE



!           LIST OF MAJOR VARIABLES
!           -----------------------

!           CENTR  - MID POINT OF THE INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTERVAL
!           ABSC   - ABSCISSA
!           FVALUE - FUNCTION VALUE
!           RESG   - RESULT OF THE GAUSS FORMULA
!           RESK   - RESULT OF THE KRONROD FORMULA
!           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
!                    I.E. TO I/(B-A)
!     Note: xrule(:) is an array of integers in the module
!     set_gk_data.  Its value, via the call, is used to 
!     set the size of the automatic arrays XVALUE, FVALUE.

        REAL (dp) centr, hlgth, resg, resk, reskh
        INTEGER(kind=i4) ng, nk, nt
        REAL (dp), POINTER :: fk1(:), fk2(:), fg1(:), fg2(:)
        REAL (dp), TARGET, DIMENSION (xrule(extype%key)) :: xvalue, fvalue

      ASSOCIATE(key => extype%key, abserr => extype%abserr)

      ASSOCIATE(WG =>  XWG(gaussIndex(key):gaussIndex(key+1)-1),&
                WGK => XWGK(gkIndex(key):gkIndex(key+1)-1),&
                XGK => XXGK(gkIndex(key):gkIndex(key+1)-1),&
                RULE => XRULE(key))

        nt = rule
        nk = nt/2
        ng = nk/2

        centr = half*(a+b)
        hlgth = half*(b-a)

! Compute all the x-values and then call the user-supplied function
! to compute all the f values in a single call
        xvalue = (/ centr, centr - hlgth*xgk(1:nk), centr + hlgth*xgk(1:nk) /)

        fvalue = f(xvalue,extype)

        fk1 => fvalue(2:nk+1)
        fk2 => fvalue(nk+2:nt)
        fg1 => fvalue(2:ng+1)
        fg2 => fvalue(nk+2:nk+ng+1)

!           COMPUTE THE KRONROD APPROXIMATION TO
!           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.


        IF (mod(nt/2,2)==1) THEN
          resg = fvalue(1)*wg(ng+1)
        ELSE
          resg = zero
        END IF

        resg = resg + dot_product(wg(1:ng),fg1+fg2)

        resk = dot_product(wgk(1:nk),fk1+fk2) + fvalue(1)*wgk(nk+1)

        resabs = dot_product(wgk(1:nk),abs(fk1)+abs(fk2)) + &
          abs(fvalue(1))*wgk(nk+1)
        reskh = resk*half
        resasc = dot_product(wgk(1:nk),abs(fk1-reskh)+abs(fk2-reskh)) + &
          abs(fvalue(1)-reskh)*wgk(nk+1)
      END ASSOCIATE
        intapprox = resk*hlgth
        resabs = resabs*abs(hlgth)
        resasc = resasc*abs(hlgth)
        abserr = abs((resk-resg)*hlgth)
        IF (resasc/=zero .AND. abserr/=zero) THEN
          abserr = resasc*min(one,(twohun*abserr/resasc)**onept5)
        END IF
        IF (resabs*fifty*epsilon(one)>tiny(one)) THEN
          abserr = max((epsilon(one)*fifty)*resabs,abserr)
        END IF
      END ASSOCIATE
      END SUBROUTINE qkgens

      RECURSIVE SUBROUTINE qkgenv(f,extype,a,b,intapprox,resabs,resasc,norm)

        REAL (dp), INTENT (IN) :: a, b
        REAL (dp), INTENT (INOUT), DIMENSION (:) :: intapprox, resabs, resasc
      PROCEDURE(fv) :: f
        CLASS (quadpackbase), INTENT (INOUT) :: extype
      PROCEDURE(errnorm),POINTER :: norm

!            To compute I = Integral of F over (A,B), with error
!                           estimate
!                       J = Integral of ABS(F) over (A,B)
!           Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven

!           PARAMETERS
!            ON ENTRY
!              F      - Real
!                       Function subprogram defining the scalar integrand
!                       FUNCTION F(X,EXTYPE). The actual name for F needs to be
!                       Declared E X T E R N A L in the driver program or an
!                       equivalent via use association.  Note that the argument
!                       X(:) is an array so the result is an array of 
!                       size (SIZE(X), EXTYPE%NCOMPONENTS).

!              EXTYPE   The extended or base data type of class
!                       QUADPACKBASE.  This allows a user program
!                       to pass data to the evaluation function, F.            

!              A      - Real
!                       Lower limit of integration

!              B      - Real
!                       Upper limit of integration

!            ON RETURN
!              intApprox - Real
!                       Approximation to the integral I
!                       intApprox is computed by applying a
!                       KRONROD RULE (RESK) obtained by optimal addition
!                       of abscissae to a GAUSS RULE (RESG).

!              ABSERR - Real
!                       Estimate of the modulus of the absolute error,
!                       which should not exceed ABS(I-intApprox)

!              RESABS - Real
!                       Approximation to the integral J

!              RESASC - Real
!                       Approximation to the integral of ABS(F-I/(B-A))
!                       over (A,B)


        REAL (dp) centr, hlgth, resg(size(intapprox)), resk(size(intapprox)), &
          reskh(size(intapprox))
        INTEGER(kind=i4) ng, nk, nt, i
        REAL (dp), POINTER :: fk1(:,:), fk2(:,:), fg1(:,:), fg2(:,:)


!           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
!           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
!           CORRESPONDING WEIGHTS ARE GIVEN.

!           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
!                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
!                    GAUSS RULE
!                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!                    ADDED TO THE 10-POINT GAUSS RULE

!           WGK    - WEIGHTS OF THE KRONROD RULE

!           WG     - WEIGHTS OF THE GAUSS RULE



!           LIST OF MAJOR VARIABLES
!           -----------------------

!           CENTR  - MID POINT OF THE INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTERVAL
!           ABSC   - ABSCISSA
!           FVALUE - FUNCTION VALUE
!           RESG   - RESULT OF THE GAUSS FORMULA
!           RESK   - RESULT OF THE KRONROD FORMULA
!           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
!                    I.E. TO I/(B-A)
!     Note: xrule(:) is an array of integers in the module
!     set_gk_data.  Its value, via the call, is used to 
!     set the size of the automatic arrays XVALUE, FVALUE.

        REAL (dp), TARGET, DIMENSION (xrule(extype%key),size(intapprox)) :: &
          xvalue(xrule(extype%key)), fvalue

! It would be cleaner to amalgamate the two ASSOCIATE statements but the NAG
! compiler does not like this
      ASSOCIATE(key => extype%key, abserr => extype%abserr)

      ASSOCIATE(WG =>  XWG(gaussIndex(key):gaussIndex(key+1)-1),&
                WGK => XWGK(gkIndex(key):gkIndex(key+1)-1),&
                XGK => XXGK(gkIndex(key):gkIndex(key+1)-1),&
                RULE => XRULE(key))

        nt = rule
        nk = nt/2
        ng = nk/2

        centr = half*(a+b)
        hlgth = half*(b-a)

! Compute all the x-values and then call the user supplied function
! to compute all the f values in a single call
        xvalue = (/ centr, centr - hlgth*xgk(1:nk), centr + hlgth*xgk(1:nk) /)

        fvalue = f(xvalue,extype)

        fk1 => fvalue(2:nk+1,:)
        fk2 => fvalue(nk+2:nt,:)
        fg1 => fvalue(2:ng+1,:)
        fg2 => fvalue(nk+2:nk+ng+1,:)

! Compute the Kronrod approximation to the integral, and 
! estimate the absolute error.


        IF (mod(nt/2,2)==1) THEN
          resg = fvalue(1,:)*wg(ng+1)
        ELSE
          resg = zero
        END IF

        resg = resg + matmul(wg(1:ng),(fg1+fg2))
        resk = matmul(wgk(1:nk),(fk1+fk2)) + fvalue(1,:)*wgk(nk+1)
        resabs = matmul(wgk(1:nk),(abs(fk1)+abs(fk2))) + &
          abs(fvalue(1,:))*wgk(nk+1)
        reskh = resk*half

        DO i = 1, nk
          fk1(i,:) = fk1(i,:) - reskh
          fk2(i,:) = fk2(i,:) - reskh
        END DO

        resasc = matmul(wgk(1:nk),abs(fk1)+abs(fk2)) + &
          abs(fvalue(1,:)-reskh)*wgk(nk+1)
      END ASSOCIATE
        intapprox = resk*hlgth
        resabs = resabs*abs(hlgth)
        resasc = resasc*abs(hlgth)
        abserr = norm(abs((resk-resg)*hlgth),extype)

        IF (all(resasc/=zero) .AND. abserr/=zero) THEN
          abserr = norm(resasc*min(one,(twohun*abserr/ &
            maxval(resasc))**onept5),extype)
        END IF
        IF (norm(resabs,extype)*fifty*epsilon(one)>tiny(one)) THEN
          abserr = max((epsilon(one)*fifty)*norm(resabs,extype),abserr)
        END IF
      END ASSOCIATE
      END SUBROUTINE qkgenv

      SUBROUTINE qsort_2003(limit,last,maxerr,ermax,elist,iord,nrmax)

!*****************************************************************************80

!! QSORT_2003 maintains the order of a list of local error estimates.

!  Discussion:

!    This routine maintains the descending ordering in the list of the 
!    local error estimates resulting from the interval subdivision process. 
!    At each call two error estimates are inserted using the sequential 
!    search top-down for the largest error estimate and bottom-up for the
!    smallest error estimate.

!  Author:

!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner

!  Reference:

!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983

!  Parameters:

!    Input, integer LIMIT, the maximum number of error estimates the list can
!    contain.

!    Input, integer LAST, the current number of error estimates.

!    Input/output, integer MAXERR, the index in the list of the NRMAX-th 
!    largest error.

!    Output, real ERMAX, the NRMAX-th largest error = ELIST(MAXERR).

!    Input, real ELIST(LIMIT), contains the error estimates.

!    Input/output, integer IORD(LAST).  The first K elements contain 
!    pointers to the error estimates such that ELIST(IORD(1)) through
!    ELIST(IORD(K)) form a decreasing sequence, with
!      K = LAST 
!    if 
!      LAST <= (LIMIT/2+2), 
!    and otherwise
!      K = LIMIT+1-LAST.

!    Input/output, integer NRMAX.

        INTEGER(kind=i4), INTENT (IN) :: limit, last
        INTEGER(kind=i4), INTENT (INOUT) :: iord(last), maxerr, nrmax
        REAL (dp), INTENT (INOUT) :: elist(last)
        REAL (dp), INTENT (OUT) :: ermax

        INTEGER(kind=i4) i, ibeg, isucc, j, jbnd
        INTEGER(kind=i4) jupbn, k, istate
        REAL (dp) errmax, errmin


!  Check whether the list contains more than two error estimates.

        istate = 1
        DO
          SELECT CASE (istate)
          CASE (1)
            IF (last<=2) THEN
              iord(1) = 1
              iord(2) = 2
              istate = 90
              CYCLE
            END IF

!  This part of the routine is only executed if, due to a
!  difficult integrand, subdivision increased the error
!  estimate. In the normal case the insert procedure should
!  start after the nrmax-th largest error estimate.

            errmax = elist(maxerr)

            DO i = 1, nrmax - 1

              isucc = iord(nrmax-1)

              IF (errmax<=elist(isucc)) THEN
                EXIT
              END IF

              iord(nrmax) = isucc
              nrmax = nrmax - 1

            END DO

!  Compute the number of elements in the list to be maintained
!  in descending order.  This number depends on the number of
!  subdivisions still allowed.

            jupbn = last

            IF ((limit/2+2)<last) THEN
              jupbn = limit + 3 - last
            END IF

            errmin = elist(last)

!  Insert errmax by traversing the list top-down, starting
!  comparison from the element elist(iord(nrmax+1)).

            jbnd = jupbn - 1
            ibeg = nrmax + 1

            DO i = ibeg, jbnd
              isucc = iord(i)
              IF (elist(isucc)<=errmax) THEN
                istate = 60
                CYCLE
              END IF
              iord(i-1) = isucc
            END DO

            iord(jbnd) = maxerr
            iord(jupbn) = last
            istate = 90
            CYCLE

!  Insert errmin by traversing the list bottom-up.

          CASE (60)

            iord(i-1) = maxerr
            k = jbnd

            DO j = i, jbnd
              isucc = iord(k)
              IF (errmin<elist(isucc)) THEN
                istate = 80
              END IF
              iord(k+1) = isucc
              k = k - 1
            END DO

            iord(i) = last
            istate = 90
            CYCLE

          CASE (80)

            iord(k+1) = last

!  Set maxerr and ermax.

            istate = 90
            CYCLE

          CASE (90)

            maxerr = iord(nrmax)
            ermax = elist(maxerr)
            EXIT
          END SELECT
        END DO
      END SUBROUTINE qsort_2003


      FUNCTION maxnorm(e,extend) RESULT (norm)
! This is the default norm used for vector integration.
! It is the max abs norm.  
        REAL (dp), INTENT (IN) :: e(:)

        CLASS (quadpackbase), INTENT (INOUT) :: extend
        REAL (dp) norm

        norm = maxval(abs(e))
      END FUNCTION maxnorm

    END MODULE quadpack2003
