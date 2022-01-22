!** DDEABM
SUBROUTINE DDEABM(DF,Neq,T,Y,Tout,Info,Rtol,Atol,Idid,Rwork,Lrw,Iwork,Liw)
    use mod_kinds, only : i4,dp
  !> Solve an initial value problem in ordinary differential
  !  equations using an Adams-Bashforth method.
  !***
  ! **Library:**   SLATEC (DEPAC)
  !***
  ! **Category:**  I1A1B
  !***
  ! **Type:**      DOUBLE PRECISION (DEABM-S, DDEABM-D)
  !***
  ! **Keywords:**  ADAMS-BASHFORTH METHOD, DEPAC, INITIAL VALUE PROBLEMS,
  !             ODE, ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR
  !***
  ! **Author:**  Shampine, L. F., (SNLA)
  !           Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !   This is the Adams code in the package of differential equation
  !   solvers DEPAC, consisting of the codes DDERKF, DDEABM, and DDEBDF.
  !   Design of the package was by L. F. Shampine and H. A. Watts.
  !   It is documented in
  !        SAND79-2374, DEPAC - Design of a User Oriented Package of ODE
  !                              Solvers.
  !   DDEABM is a driver for a modification of the code ODE written by
  !             L. F. Shampine and M. K. Gordon
  !             Sandia Laboratories
  !             Albuquerque, New Mexico 87185
  !
  !- *********************************************************************
  !-  ABSTRACT *
  !- ***********
  !
  !   Subroutine DDEABM uses the Adams-Bashforth-Moulton
  !   Predictor-Corrector formulas of orders one through twelve to
  !   integrate a system of NEQ first order ordinary differential
  !   equations of the form
  !                         DU/DX = DF(X,U)
  !   when the vector Y(*) of initial values for U(*) at X=T is given.
  !   The subroutine integrates from T to TOUT. It is easy to continue the
  !   integration to get results at additional TOUT.  This is the interval
  !   mode of operation.  It is also easy for the routine to return with
  !   the solution at each intermediate step on the way to TOUT.  This is
  !   the intermediate-output mode of operation.
  !
  !   DDEABM uses subprograms DDES, DSTEPS, DINTP, DHSTRT, DHVNRM,
  !   D1MACH, and the error handling routine XERMSG.  The only machine
  !   dependent parameters to be assigned appear in D1MACH.
  !
  !- *********************************************************************
  !-  Description of The Arguments To DDEABM (An Overview) *
  !- *********************************************************************
  !
  !   The Parameters are
  !
  !      DF -- This is the name of a subroutine which you provide to
  !             define the differential equations.
  !
  !      NEQ -- This is the number of (first order) differential
  !             equations to be integrated.
  !
  !      T -- This is a DOUBLE PRECISION value of the independent
  !           variable.
  !
  !      Y(*) -- This DOUBLE PRECISION array contains the solution
  !              components at T.
  !
  !      TOUT -- This is a DOUBLE PRECISION point at which a solution is
  !              desired.
  !
  !      INFO(*) -- The basic task of the code is to integrate the
  !             differential equations from T to TOUT and return an
  !             answer at TOUT.  INFO(*) is an INTEGER(i4) array which is used
  !             to communicate exactly how you want this task to be
  !             carried out.
  !
  !      RTOL, ATOL -- These DOUBLE PRECISION quantities represent
  !                    relative and absolute error tolerances which you
  !                    provide to indicate how accurately you wish the
  !                    solution to be computed.  You may choose them to be
  !                    both scalars or else both vectors.
  !
  !      IDID -- This scalar quantity is an indicator reporting what
  !             the code did.  You must monitor this INTEGER(i4) variable to
  !             decide what action to take next.
  !
  !      RWORK(*), LRW -- RWORK(*) is a DOUBLE PRECISION work array of
  !             length LRW which provides the code with needed storage
  !             space.
  !
  !      IWORK(*), LIW -- IWORK(*) is an INTEGER(i4) work array of length LIW
  !             which provides the code with needed storage space and an
  !             across call flag.
  !
  !      RPAR, IPAR -- These are DOUBLE PRECISION and INTEGER(i4) parameter
  !             arrays which you can use for communication between your
  !             calling program and the DF subroutine.
  !
  !  Quantities which are used as input items are
  !             NEQ, T, Y(*), TOUT, INFO(*),
  !             RTOL, ATOL, RWORK(1), LRW and LIW.
  !
  !  Quantities which may be altered by the code are
  !             T, Y(*), INFO(1), RTOL, ATOL,
  !             IDID, RWORK(*) and IWORK(*).
  !
  !- *********************************************************************
  !-  INPUT -- What To Do On The First Call To DDEABM *
  !- *********************************************************************
  !
  !   The first call of the code is defined to be the start of each new
  !   problem.  Read through the descriptions of all the following items,
  !   provide sufficient storage space for designated arrays, set
  !   appropriate variables for the initialization of the problem, and
  !   give information about how you want the problem to be solved.
  !
  !
  !      DF -- Provide a subroutine of the form
  !                               DF(X,U,UPRIME,RPAR,IPAR)
  !             to define the system of first order differential equations
  !             which is to be solved.  For the given values of X and the
  !             vector  U(*)=(U(1),U(2),...,U(NEQ)), the subroutine must
  !             evaluate the NEQ components of the system of differential
  !             equations  DU/DX=DF(X,U)  and store the derivatives in the
  !             array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for
  !             equations I=1,...,NEQ.
  !
  !             Subroutine DF must NOT alter X or U(*).  You must declare
  !             the name df in an external statement in your program that
  !             calls DDEABM.  You must dimension U and UPRIME in DF.
  !
  !             RPAR and IPAR are DOUBLE PRECISION and INTEGER(i4) parameter
  !             arrays which you can use for communication between your
  !             calling program and subroutine DF. They are not used or
  !             altered by DDEABM.  If you do not need RPAR or IPAR,
  !             ignore these parameters by treating them as dummy
  !             arguments. If you do choose to use them, dimension them in
  !             your calling program and in DF as arrays of appropriate
  !             length.
  !
  !      NEQ -- Set it to the number of differential equations.
  !             (NEQ >= 1)
  !
  !      T -- Set it to the initial point of the integration.
  !             You must use a program variable for T because the code
  !             changes its value.
  !
  !      Y(*) -- Set this vector to the initial values of the NEQ solution
  !             components at the initial point.  You must dimension Y at
  !             least NEQ in your calling program.
  !
  !      TOUT -- Set it to the first point at which a solution
  !             is desired.  You can take TOUT = T, in which case the code
  !             will evaluate the derivative of the solution at T and
  !             return. Integration either forward in T  (TOUT > T)  or
  !             backward in T  (TOUT < T)  is permitted.
  !
  !             The code advances the solution from T to TOUT using
  !             step sizes which are automatically selected so as to
  !             achieve the desired accuracy.  If you wish, the code will
  !             return with the solution and its derivative following
  !             each intermediate step (intermediate-output mode) so that
  !             you can monitor them, but you still must provide TOUT in
  !             accord with the basic aim of the code.
  !
  !             The first step taken by the code is a critical one
  !             because it must reflect how fast the solution changes near
  !             the initial point.  The code automatically selects an
  !             initial step size which is practically always suitable for
  !             the problem. By using the fact that the code will not step
  !             past TOUT in the first step, you could, if necessary,
  !             restrict the length of the initial step size.
  !
  !             For some problems it may not be permissible to integrate
  !             past a point TSTOP because a discontinuity occurs there
  !             or the solution or its derivative is not defined beyond
  !             TSTOP.  When you have declared a TSTOP point (see INFO(4)
  !             and RWORK(1)), you have told the code not to integrate
  !             past TSTOP.  In this case any TOUT beyond TSTOP is invalid
  !             input.
  !
  !      INFO(*) -- Use the INFO array to give the code more details about
  !             how you want your problem solved.  This array should be
  !             dimensioned of length 15 to accommodate other members of
  !             DEPAC or possible future extensions, though DDEABM uses
  !             only the first four entries.  You must respond to all of
  !             the following items which are arranged as questions.  The
  !             simplest use of the code corresponds to answering all
  !             questions as YES ,i.e. setting ALL entries of INFO to 0.
  !
  !        INFO(1) -- This parameter enables the code to initialize
  !               itself.  You must set it to indicate the start of every
  !               new problem.
  !
  !            **** Is this the first call for this problem ...
  !                  YES -- set INFO(1) = 0
  !                   NO -- not applicable here.
  !                         See below for continuation calls.  ****
  !
  !        INFO(2) -- How much accuracy you want of your solution
  !               is specified by the error tolerances RTOL and ATOL.
  !               The simplest use is to take them both to be scalars.
  !               To obtain more flexibility, they can both be vectors.
  !               The code must be told your choice.
  !
  !            **** Are both error tolerances RTOL, ATOL scalars ...
  !                  YES -- set INFO(2) = 0
  !                         and input scalars for both RTOL and ATOL
  !                   NO -- set INFO(2) = 1
  !                         and input arrays for both RTOL and ATOL ****
  !
  !        INFO(3) -- The code integrates from T in the direction
  !               of TOUT by steps.  If you wish, it will return the
  !               computed solution and derivative at the next
  !               intermediate step (the intermediate-output mode) or
  !               TOUT, whichever comes first.  This is a good way to
  !               proceed if you want to see the behavior of the solution.
  !               If you must have solutions at a great many specific
  !               TOUT points, this code will compute them efficiently.
  !
  !            **** Do you want the solution only at
  !                 TOUT (and not at the next intermediate step) ...
  !                  YES -- set INFO(3) = 0
  !                   NO -- set INFO(3) = 1 ****
  !
  !        INFO(4) -- To handle solutions at a great many specific
  !               values TOUT efficiently, this code may integrate past
  !               TOUT and interpolate to obtain the result at TOUT.
  !               Sometimes it is not possible to integrate beyond some
  !               point TSTOP because the equation changes there or it is
  !               not defined past TSTOP.  Then you must tell the code
  !               not to go past.
  !
  !            **** Can the integration be carried out without any
  !                 Restrictions on the independent variable T ...
  !                  YES -- set INFO(4)=0
  !                   NO -- set INFO(4)=1
  !                         and define the stopping point TSTOP by
  !                         setting RWORK(1)=TSTOP ****
  !
  !      RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL)
  !             error tolerances to tell the code how accurately you want
  !             the solution to be computed.  They must be defined as
  !             program variables because the code may change them.  You
  !             have two choices --
  !                  Both RTOL and ATOL are scalars. (INFO(2)=0)
  !                  Both RTOL and ATOL are vectors. (INFO(2)=1)
  !             In either case all components must be non-negative.
  !
  !             The tolerances are used by the code in a local error test
  !             at each step which requires roughly that
  !                     ABS(LOCAL ERROR) <= RTOL*ABS(Y)+ATOL
  !             for each vector component.
  !             (More specifically, a Euclidean norm is used to measure
  !             the size of vectors, and the error test uses the magnitude
  !             of the solution at the beginning of the step.)
  !
  !             The true (global) error is the difference between the true
  !             solution of the initial value problem and the computed
  !             approximation.  Practically all present day codes,
  !             including this one, control the local error at each step
  !             and do not even attempt to control the global error
  !             directly.  Roughly speaking, they produce a solution Y(T)
  !             which satisfies the differential equations with a
  !             residual R(T),    DY(T)/DT = DF(T,Y(T)) + R(T)   ,
  !             and, almost always, R(T) is bounded by the error
  !             tolerances.  Usually, but not always, the true accuracy of
  !             the computed Y is comparable to the error tolerances. This
  !             code will usually, but not always, deliver a more accurate
  !             solution if you reduce the tolerances and integrate again.
  !             By comparing two such solutions you can get a fairly
  !             reliable idea of the true error in the solution at the
  !             bigger tolerances.
  !
  !             Setting ATOL=0.D0 results in a pure relative error test on
  !             that component. Setting RTOL=0. results in a pure absolute
  !             error test on that component.  A mixed test with non-zero
  !             RTOL and ATOL corresponds roughly to a relative error
  !             test when the solution component is much bigger than ATOL
  !             and to an absolute error test when the solution component
  !             is smaller than the threshold ATOL.
  !
  !             Proper selection of the absolute error control parameters
  !             ATOL  requires you to have some idea of the scale of the
  !             solution components.  To acquire this information may mean
  !             that you will have to solve the problem more than once. In
  !             the absence of scale information, you should ask for some
  !             relative accuracy in all the components (by setting  RTOL
  !             values non-zero) and perhaps impose extremely small
  !             absolute error tolerances to protect against the danger of
  !             a solution component becoming zero.
  !
  !             The code will not attempt to compute a solution at an
  !             accuracy unreasonable for the machine being used.  It will
  !             advise you if you ask for too much accuracy and inform
  !             you as to the maximum accuracy it believes possible.
  !
  !      RWORK(*) -- Dimension this DOUBLE PRECISION work array of length
  !             LRW in your calling program.
  !
  !      RWORK(1) -- If you have set INFO(4)=0, you can ignore this
  !             optional input parameter.  Otherwise you must define a
  !             stopping point TSTOP by setting   RWORK(1) = TSTOP.
  !             (for some problems it may not be permissible to integrate
  !             past a point TSTOP because a discontinuity occurs there
  !             or the solution or its derivative is not defined beyond
  !             TSTOP.)
  !
  !      LRW -- Set it to the declared length of the RWORK array.
  !             You must have  LRW >= 130+21*NEQ
  !
  !      IWORK(*) -- Dimension this INTEGER(i4) work array of length LIW in
  !             your calling program.
  !
  !      LIW -- Set it to the declared length of the IWORK array.
  !             You must have  LIW >= 51
  !
  !      RPAR, IPAR -- These are parameter arrays, of DOUBLE PRECISION and
  !             INTEGER(i4) type, respectively.  You can use them for
  !             communication between your program that calls DDEABM and
  !             the  DF subroutine.  They are not used or altered by
  !             DDEABM.  If you do not need RPAR or IPAR, ignore these
  !             parameters by treating them as dummy arguments.  If you do
  !             choose to use them, dimension them in your calling program
  !             and in DF as arrays of appropriate length.
  !
  !- *********************************************************************
  !-  OUTPUT -- After Any Return From DDEABM *
  !- *********************************************************************
  !
  !   The principal aim of the code is to return a computed solution at
  !   TOUT, although it is also possible to obtain intermediate results
  !   along the way.  To find out whether the code achieved its goal
  !   or if the integration process was interrupted before the task was
  !   completed, you must check the IDID parameter.
  !
  !
  !      T -- The solution was successfully advanced to the
  !             output value of T.
  !
  !      Y(*) -- Contains the computed solution approximation at T.
  !             You may also be interested in the approximate derivative
  !             of the solution at T.  It is contained in
  !             RWORK(21),...,RWORK(20+NEQ).
  !
  !      IDID -- Reports what the code did
  !
  !                         *** Task Completed ***
  !                   Reported by positive values of IDID
  !
  !             IDID = 1 -- A step was successfully taken in the
  !                       intermediate-output mode.  The code has not
  !                       yet reached TOUT.
  !
  !             IDID = 2 -- The integration to TOUT was successfully
  !                       completed (T=TOUT) by stepping exactly to TOUT.
  !
  !             IDID = 3 -- The integration to TOUT was successfully
  !                       completed (T=TOUT) by stepping past TOUT.
  !                       Y(*) is obtained by interpolation.
  !
  !                         *** Task Interrupted ***
  !                   Reported by negative values of IDID
  !
  !             IDID = -1 -- A large amount of work has been expended.
  !                       (500 steps attempted)
  !
  !             IDID = -2 -- The error tolerances are too stringent.
  !
  !             IDID = -3 -- The local error test cannot be satisfied
  !                       because you specified a zero component in ATOL
  !                       and the corresponding computed solution
  !                       component is zero.  Thus, a pure relative error
  !                       test is impossible for this component.
  !
  !             IDID = -4 -- The problem appears to be stiff.
  !
  !             IDID = -5,-6,-7,..,-32  -- Not applicable for this code
  !                       but used by other members of DEPAC or possible
  !                       future extensions.
  !
  !                         *** Task Terminated ***
  !                   Reported by the value of IDID=-33
  !
  !             IDID = -33 -- The code has encountered trouble from which
  !                       it cannot recover.  A message is printed
  !                       explaining the trouble and control is returned
  !                       to the calling program. For example, this occurs
  !                       when invalid input is detected.
  !
  !      RTOL, ATOL -- These quantities remain unchanged except when
  !             IDID = -2.  In this case, the error tolerances have been
  !             increased by the code to values which are estimated to be
  !             appropriate for continuing the integration.  However, the
  !             reported solution at T was obtained using the input values
  !             of RTOL and ATOL.
  !
  !      RWORK, IWORK -- Contain information which is usually of no
  !             interest to the user but necessary for subsequent calls.
  !             However, you may find use for
  !
  !             RWORK(11)--which contains the step size H to be
  !                        attempted on the next step.
  !
  !             RWORK(12)--if the tolerances have been increased by the
  !                        code (IDID = -2), they were multiplied by the
  !                        value in RWORK(12).
  !
  !             RWORK(13)--Which contains the current value of the
  !                        independent variable, i.e. the farthest point
  !                        integration has reached. This will be different
  !                        from T only when interpolation has been
  !                        performed (IDID=3).
  !
  !             RWORK(20+I)--Which contains the approximate derivative
  !                        of the solution component Y(I).  In DDEABM, it
  !                        is obtained by calling subroutine DF to
  !                        evaluate the differential equation using T and
  !                        Y(*) when IDID=1 or 2, and by interpolation
  !                        when IDID=3.
  !
  !- *********************************************************************
  !-  INPUT -- What To Do To Continue The Integration *
  !-              (calls after the first)             *
  !- *********************************************************************
  !
  !        This code is organized so that subsequent calls to continue the
  !        integration involve little (if any) additional effort on your
  !        part. You must monitor the IDID parameter in order to determine
  !        what to do next.
  !
  !        Recalling that the principal task of the code is to integrate
  !        from T to TOUT (the interval mode), usually all you will need
  !        to do is specify a new TOUT upon reaching the current TOUT.
  !
  !        Do not alter any quantity not specifically permitted below,
  !        in particular do not alter NEQ, T, Y(*), RWORK(*), IWORK(*) or
  !        the differential equation in subroutine DF. Any such alteration
  !        constitutes a new problem and must be treated as such, i.e.
  !        you must start afresh.
  !
  !        You cannot change from vector to scalar error control or vice
  !        versa (INFO(2)) but you can change the size of the entries of
  !        RTOL, ATOL.  Increasing a tolerance makes the equation easier
  !        to integrate.  Decreasing a tolerance will make the equation
  !        harder to integrate and should generally be avoided.
  !
  !        You can switch from the intermediate-output mode to the
  !        interval mode (INFO(3)) or vice versa at any time.
  !
  !        If it has been necessary to prevent the integration from going
  !        past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
  !        code will not integrate to any TOUT beyond the currently
  !        specified TSTOP.  Once TSTOP has been reached you must change
  !        the value of TSTOP or set INFO(4)=0.  You may change INFO(4)
  !        or TSTOP at any time but you must supply the value of TSTOP in
  !        RWORK(1) whenever you set INFO(4)=1.
  !
  !        The parameter INFO(1) is used by the code to indicate the
  !        beginning of a new problem and to indicate whether integration
  !        is to be continued.  You must input the value  INFO(1) = 0
  !        when starting a new problem.  You must input the value
  !        INFO(1) = 1  if you wish to continue after an interrupted task.
  !        Do not set  INFO(1) = 0  on a continuation call unless you
  !        want the code to restart at the current T.
  !
  !                         *** Following A Completed Task ***
  !         If
  !             IDID = 1, call the code again to continue the integration
  !                     another step in the direction of TOUT.
  !
  !             IDID = 2 or 3, define a new TOUT and call the code again.
  !                     TOUT must be different from T. You cannot change
  !                     the direction of integration without restarting.
  !
  !                         *** Following An Interrupted Task ***
  !                     To show the code that you realize the task was
  !                     interrupted and that you want to continue, you
  !                     must take appropriate action and reset INFO(1) = 1
  !         If
  !             IDID = -1, the code has attempted 500 steps.
  !                     If you want to continue, set INFO(1) = 1 and
  !                     call the code again. An additional 500 steps
  !                     will be allowed.
  !
  !             IDID = -2, the error tolerances RTOL, ATOL have been
  !                     increased to values the code estimates appropriate
  !                     for continuing.  You may want to change them
  !                     yourself.  If you are sure you want to continue
  !                     with relaxed error tolerances, set INFO(1)=1 and
  !                     call the code again.
  !
  !             IDID = -3, a solution component is zero and you set the
  !                     corresponding component of ATOL to zero.  If you
  !                     are sure you want to continue, you must first
  !                     alter the error criterion to use positive values
  !                     for those components of ATOL corresponding to zero
  !                     solution components, then set INFO(1)=1 and call
  !                     the code again.
  !
  !             IDID = -4, the problem appears to be stiff.  It is very
  !                     inefficient to solve such problems with DDEABM.
  !                     The code DDEBDF in DEPAC handles this task
  !                     efficiently.  If you are absolutely sure you want
  !                     to continue with DDEABM, set INFO(1)=1 and call
  !                     the code again.
  !
  !             IDID = -5,-6,-7,..,-32  --- cannot occur with this code
  !                     but used by other members of DEPAC or possible
  !                     future extensions.
  !
  !                         *** Following A Terminated Task ***
  !         If
  !             IDID = -33, you cannot continue the solution of this
  !                     problem.  An attempt to do so will result in your
  !                     run being terminated.
  !
  !- *********************************************************************
  !- Long Description:
  !
  !- *********************************************************************
  !-              DEPAC Package Overview           *
  !- *********************************************************************
  !
  ! ....   You have a choice of three differential equation solvers from
  ! ....   DEPAC. The following brief descriptions are meant to aid you in
  ! ....   choosing the most appropriate code for your problem.
  !
  ! ....   DDERKF is a fifth order Runge-Kutta code. It is the simplest of
  ! ....   the three choices, both algorithmically and in the use of the
  ! ....   code. DDERKF is primarily designed to solve non-stiff and
  ! ....   mildly stiff differential equations when derivative evaluations
  ! ....   are not expensive. It should generally not be used to get high
  ! ....   accuracy results nor answers at a great many specific points.
  ! ....   Because DDERKF has very low overhead costs, it will usually
  ! ....   result in the least expensive integration when solving
  ! ....   problems requiring a modest amount of accuracy and having
  ! ....   equations that are not costly to evaluate. DDERKF attempts to
  ! ....   discover when it is not suitable for the task posed.
  !
  ! ....   DDEABM is a variable order (one through twelve) Adams code.
  ! ....   Its complexity lies somewhere between that of DDERKF and
  ! ....   DDEBDF.  DDEABM is primarily designed to solve non-stiff and
  ! ....   mildly stiff differential equations when derivative evaluations
  ! ....   are expensive, high accuracy results are needed or answers at
  ! ....   many specific points are required. DDEABM attempts to discover
  ! ....   when it is not suitable for the task posed.
  !
  ! ....   DDEBDF is a variable order (one through five) backward
  ! ....   differentiation formula code. it is the most complicated of
  ! ....   the three choices. DDEBDF is primarily designed to solve stiff
  ! ....   differential equations at crude to moderate tolerances.
  ! ....   If the problem is very stiff at all, DDERKF and DDEABM will be
  ! ....   quite inefficient compared to DDEBDF. However, DDEBDF will be
  ! ....   inefficient compared to DDERKF and DDEABM on non-stiff problems
  ! ....   because it uses much more storage, has a much larger overhead,
  ! ....   and the low order formulas will not give high accuracies
  ! ....   efficiently.
  !
  ! ....   The concept of stiffness cannot be described in a few words.
  ! ....   If you do not know the problem to be stiff, try either DDERKF
  ! ....   or DDEABM. Both of these codes will inform you of stiffness
  ! ....   when the cost of solving such problems becomes important.
  !
  !- ********************************************************************
  !
  !***
  ! **References:**  L. F. Shampine and H. A. Watts, DEPAC - design of a user
  !                 oriented package of ODE solvers, Report SAND79-2374,
  !                 Sandia Laboratories, 1979.
  !***
  ! **Routines called:**  DDES, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891024  Changed references from DVNORM to DHVNRM.  (WRB)
  !   891024  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTERFACE
    SUBROUTINE DF(X,U,Uprime)
      IMPORT dp
      REAL(dp), INTENT(IN) :: X
      REAL(dp), INTENT(IN) :: U(:)
      REAL(dp), INTENT(OUT) :: Uprime(:)
    END SUBROUTINE DF
  END INTERFACE
  INTEGER(i4), INTENT(IN) :: Liw, Lrw, Neq
  INTEGER(i4), INTENT(OUT) :: Idid
  INTEGER(i4), INTENT(INOUT) :: Info(15), Iwork(Liw)
  REAL(dp), INTENT(IN) :: Tout
  REAL(dp), INTENT(INOUT) :: T
  REAL(dp), INTENT(INOUT) :: Atol(:), Rtol(:), Rwork(Lrw), Y(Neq)
     !DIR$ ASSUME_ALIGNED Atol:64
     !DIR$ ASSUME_ALIGNED Rtol:64
     !DIR$ ASSUME_ALIGNED Rwork:64
     !DIR$ ASSUME_ALIGNED Y:64
  !
  INTEGER(i4) :: igi, ixold, ialpha, ibeta, idelsn, ifouru, ig, ihold, ip, iphi, &
    ipsi, isig, itold, itstar, itwou, iv, iw, iwt, iyp, iypout, iyy
  LOGICAL :: start, phase1, nornd, stiff, intout
  CHARACTER(8) :: xern1
  CHARACTER(16) :: xern3
  !
  !     CHECK FOR AN APPARENT INFINITE LOOP
  !
  !* FIRST EXECUTABLE STATEMENT  DDEABM
  IF( Info(1)==0 ) Iwork(Liw) = 0
  IF( Iwork(Liw)>=5 ) THEN
    IF( T==Rwork(21+Neq) ) THEN
      WRITE (xern3,'(1PE15.6)') T
      ERROR STOP 'DDEABM : AN APPARENT INFINITE LOOP HAS BEEN DETECTED.&
        & YOU HAVE MADE REPEATED CALLS AT T AND THE INTEGRATION HAS NOT ADVANCED.&
        & CHECK THE WAY YOU HAVE SET PARAMETERS FOR THE CALL&
        & TO THE CODE, PARTICULARLY INFO(1).'
      RETURN
    END IF
  END IF
  !
  !     CHECK LRW AND LIW FOR SUFFICIENT STORAGE ALLOCATION
  !
  Idid = 0
  IF( Lrw<130+21*Neq ) THEN
    WRITE (xern1,'(I8)') Lrw
    ERROR STOP 'DDEABM : THE LENGTH OF THE RWORK ARRAY MUST BE AT LEAST 130 + 21*NEQ.'
    Idid = -33
  END IF
  !
  IF( Liw<51 ) THEN
    WRITE (xern1,'(I8)') Liw
    ERROR STOP 'DDEABM : THE LENGTH OF THE IWORK ARRAY MUST BE AT LEAST 51.'
    Idid = -33
  END IF
  !
  !     COMPUTE THE INDICES FOR THE ARRAYS TO BE STORED IN THE WORK ARRAY
  !
  iypout = 21
  itstar = Neq + 21
  iyp = 1 + itstar
  iyy = Neq + iyp
  iwt = Neq + iyy
  ip = Neq + iwt
  iphi = Neq + ip
  ialpha = (Neq*16) + iphi
  ibeta = 12 + ialpha
  ipsi = 12 + ibeta
  iv = 12 + ipsi
  iw = 12 + iv
  isig = 12 + iw
  ig = 13 + isig
  igi = 13 + ig
  ixold = 11 + igi
  ihold = 1 + ixold
  itold = 1 + ihold
  idelsn = 1 + itold
  itwou = 1 + idelsn
  ifouru = 1 + itwou
  !
  Rwork(itstar) = T
  IF( Info(1)/=0 ) THEN
    start = Iwork(21)/=(-1)
    phase1 = Iwork(22)/=(-1)
    nornd = Iwork(23)/=(-1)
    stiff = Iwork(24)/=(-1)
    intout = Iwork(25)/=(-1)
  END IF
  !
  CALL DDES(DF,Neq,T,Y,Tout,Info,Rtol,Atol,Idid,Rwork(iypout),Rwork(iyp),&
    Rwork(iyy),Rwork(iwt),Rwork(ip),Rwork(iphi),Rwork(ialpha),&
    Rwork(ibeta),Rwork(ipsi),Rwork(iv),Rwork(iw),Rwork(isig),&
    Rwork(ig),Rwork(igi),Rwork(11),Rwork(12),Rwork(13),Rwork(ixold),&
    Rwork(ihold),Rwork(itold),Rwork(idelsn),Rwork(1),Rwork(itwou),&
    Rwork(ifouru),start,phase1,nornd,stiff,intout,Iwork(26),&
    Iwork(27),Iwork(28),Iwork(29),Iwork(30),Iwork(31),Iwork(32),&
    Iwork(33),Iwork(34),Iwork(35),Iwork(45))
  !
  Iwork(21) = -1
  IF( start ) Iwork(21) = 1
  Iwork(22) = -1
  IF( phase1 ) Iwork(22) = 1
  Iwork(23) = -1
  IF( nornd ) Iwork(23) = 1
  Iwork(24) = -1
  IF( stiff ) Iwork(24) = 1
  Iwork(25) = -1
  IF( intout ) Iwork(25) = 1
  !
  IF( Idid/=(-2) ) Iwork(Liw) = Iwork(Liw) + 1
  IF( T/=Rwork(itstar) ) Iwork(Liw) = 0
  !
END SUBROUTINE DDEABM
