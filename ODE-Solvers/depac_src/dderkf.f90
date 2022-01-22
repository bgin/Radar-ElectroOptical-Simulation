!** DDERKF
SUBROUTINE DDERKF(DF,Neq,T,Y,Tout,Info,Rtol,Atol,Idid,Rwork,Lrw,Iwork,Liw)
     use mod_kinds, only : i4, dp
  !> Solve an initial value problem in ordinary differential
  !  equations using a Runge-Kutta-Fehlberg scheme.
  !***
  ! **Library:**   SLATEC (DEPAC)
  !***
  ! **Category:**  I1A1A
  !***
  ! **Type:**      DOUBLE PRECISION (DERKF-S, DDERKF-D)
  !***
  ! **Keywords:**  DEPAC, INITIAL VALUE PROBLEMS, ODE,
  !             ORDINARY DIFFERENTIAL EQUATIONS, RKF,
  !             RUNGE-KUTTA-FEHLBERG METHODS
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !           Shampine, L. F., (SNLA)
  !***
  ! **Description:**
  !
  !   This is the Runge-Kutta code in the package of differential equation
  !   solvers DEPAC, consisting of the codes DDERKF, DDEABM, and DDEBDF.
  !   Design of the package was by L. F. Shampine and H. A. Watts.
  !   It is documented in
  !        SAND-79-2374, DEPAC - Design of a User Oriented Package of ODE
  !                              Solvers.
  !   DDERKF is a driver for a modification of the code RKF45 written by
  !             H. A. Watts and L. F. Shampine
  !             Sandia Laboratories
  !             Albuquerque, New Mexico 87185
  !
  !- *********************************************************************
  !- *            DDEPAC PACKAGE OVERVIEW           **
  !- *********************************************************************
  !
  !        You have a choice of three differential equation solvers from
  !        DDEPAC.  The following brief descriptions are meant to aid you
  !        in choosing the most appropriate code for your problem.
  !
  !        DDERKF is a fifth order Runge-Kutta code. It is the simplest of
  !        the three choices, both algorithmically and in the use of the
  !        code. DDERKF is primarily designed to solve non-stiff and mild-
  !        ly stiff differential equations when derivative evaluations are
  !        not expensive.  It should generally not be used to get high
  !        accuracy results nor answers at a great many specific points.
  !        Because DDERKF has very low overhead costs, it will usually
  !        result in the least expensive integration when solving
  !        problems requiring a modest amount of accuracy and having
  !        equations that are not costly to evaluate.  DDERKF attempts to
  !        discover when it is not suitable for the task posed.
  !
  !        DDEABM is a variable order (one through twelve) Adams code. Its
  !        complexity lies somewhere between that of DDERKF and DDEBDF.
  !        DDEABM is primarily designed to solve non-stiff and mildly
  !        stiff differential equations when derivative evaluations are
  !        expensive, high accuracy results are needed or answers at
  !        many specific points are required.  DDEABM attempts to discover
  !        when it is not suitable for the task posed.
  !
  !        DDEBDF is a variable order (one through five) backward
  !        differentiation formula code.  It is the most complicated of
  !        the three choices.  DDEBDF is primarily designed to solve stiff
  !        differential equations at crude to moderate tolerances.
  !        If the problem is very stiff at all, DDERKF and DDEABM will be
  !        quite inefficient compared to DDEBDF.  However, DDEBDF will be
  !        inefficient compared to DDERKF and DDEABM on non-stiff problems
  !        because it uses much more storage, has a much larger overhead,
  !        and the low order formulas will not give high accuracies
  !        efficiently.
  !
  !        The concept of stiffness cannot be described in a few words.
  !        If you do not know the problem to be stiff, try either DDERKF
  !        or DDEABM.  Both of these codes will inform you of stiffness
  !        when the cost of solving such problems becomes important.
  !
  !- *********************************************************************
  !- * ABSTRACT **
  !- *********************************************************************
  !
  !   Subroutine DDERKF uses a Runge-Kutta-Fehlberg (4,5) method to
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
  !   DDERKF uses subprograms DRKFS, DFEHL, DHSTRT, DHVNRM, D1MACH, and
  !   the error handling routine XERMSG. The only machine dependent
  !   parameters to be assigned appear in D1MACH.
  !
  !- *********************************************************************
  !- * DESCRIPTION OF THE ARGUMENTS TO DDERKF (AN OVERVIEW) **
  !- *********************************************************************
  !
  !   The Parameters are:
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
  !             relative and absolute error tolerances which you provide
  !             to indicate how accurately you wish the solution to be
  !             computed. You may choose them to be both scalars or else
  !             both vectors.
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
  !             RTOL, ATOL, LRW and LIW.
  !
  !  Quantities which may be altered by the code are
  !             T, Y(*), INFO(1), RTOL, ATOL,
  !             IDID, RWORK(*) and IWORK(*).
  !
  !- *********************************************************************
  !- * INPUT -- What to do On The First Call To DDERKF **
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
  !             Subroutine DF must not alter X or U(*). You must declare
  !             the name DF in an external statement in your program that
  !             calls DDERKF. You must dimension U and UPRIME in DF.
  !
  !             RPAR and IPAR are DOUBLE PRECISION and INTEGER(i4) parameter
  !             arrays which you can use for communication between your
  !             calling program and subroutine DF. They are not used or
  !             altered by DDERKF.  If you do not need RPAR or IPAR,
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
  !             return.  Integration either forward in T  (TOUT > T) or
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
  !             the problem.  By using the fact that the code will not
  !             step past TOUT in the first step, you could, if necessary,
  !             restrict the length of the initial step size.
  !
  !             For some problems it may not be permissible to integrate
  !             past a point TSTOP because a discontinuity occurs there
  !             or the solution or its derivative is not defined beyond
  !             TSTOP.  Since DDERKF will never step past a TOUT point,
  !             you need only make sure that no TOUT lies beyond TSTOP.
  !
  !      INFO(*) -- Use the INFO array to give the code more details about
  !             how you want your problem solved.  This array should be
  !             dimensioned of length 15 to accommodate other members of
  !             DEPAC or possible future extensions, though DDERKF uses
  !             only the first three entries.  You must respond to all of
  !             the following items which are arranged as questions.  The
  !             simplest use of the code corresponds to answering all
  !             questions as YES ,i.e. setting all entries of INFO to 0.
  !
  !        INFO(1) -- This parameter enables the code to initialize
  !               itself.  You must set it to indicate the start of every
  !               new problem.
  !
  !            **** Is this the first call for this problem ...
  !                  YES -- Set INFO(1) = 0
  !                   NO -- Not applicable here.
  !                         See below for continuation calls.  ****
  !
  !        INFO(2) -- How much accuracy you want of your solution
  !               is specified by the error tolerances RTOL and ATOL.
  !               The simplest use is to take them both to be scalars.
  !               To obtain more flexibility, they can both be vectors.
  !               The code must be told your choice.
  !
  !            **** Are both error tolerances RTOL, ATOL scalars ...
  !                  YES -- Set INFO(2) = 0
  !                         and input scalars for both RTOL and ATOL
  !                   NO -- Set INFO(2) = 1
  !                         and input arrays for both RTOL and ATOL ****
  !
  !        INFO(3) -- The code integrates from T in the direction
  !               of TOUT by steps.  If you wish, it will return the
  !               computed solution and derivative at the next
  !               intermediate step (the intermediate-output mode).
  !               This is a good way to proceed if you want to see the
  !               behavior of the solution.  If you must have solutions at
  !               a great many specific TOUT points, this code is
  !               INEFFICIENT.  The code DDEABM in DEPAC handles this task
  !               more efficiently.
  !
  !            **** Do you want the solution only at
  !                 TOUT (and not at the next intermediate step) ...
  !                  YES -- Set INFO(3) = 0
  !                   NO -- Set INFO(3) = 1 ****
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
  !             (More specifically, a maximum norm is used to measure
  !             the size of vectors, and the error test uses the average
  !             of the magnitude of the solution at the beginning and end
  !             of the step.)
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
  !             Setting ATOL=0. results in a pure relative error test on
  !             that component.  Setting RTOL=0. yields a pure absolute
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
  !             If you want relative accuracies smaller than about
  !             10**(-8), you should not ordinarily use DDERKF. The code
  !             DDEABM in DEPAC obtains stringent accuracies more
  !             efficiently.
  !
  !      RWORK(*) -- Dimension this DOUBLE PRECISION work array of length
  !             LRW in your calling program.
  !
  !      LRW -- Set it to the declared length of the RWORK array.
  !             You must have  LRW >= 33+7*NEQ
  !
  !      IWORK(*) -- Dimension this INTEGER(i4) work array of length LIW in
  !             your calling program.
  !
  !      LIW -- Set it to the declared length of the IWORK array.
  !             You must have  LIW >= 34
  !
  !      RPAR, IPAR -- These are parameter arrays, of DOUBLE PRECISION and
  !             INTEGER(i4) type, respectively. You can use them for
  !             communication between your program that calls DDERKF and
  !             the  DF subroutine. They are not used or altered by
  !             DDERKF. If you do not need RPAR or IPAR, ignore these
  !             parameters by treating them as dummy arguments. If you do
  !             choose to use them, dimension them in your calling program
  !             and in DF as arrays of appropriate length.
  !
  !- *********************************************************************
  !- * OUTPUT -- After any return from DDERKF **
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
  !             IDID = -5 -- DDERKF is being used very inefficiently
  !                       because the natural step size is being
  !                       restricted by too frequent output.
  !
  !             IDID = -6,-7,..,-32  -- Not applicable for this code but
  !                       used by other members of DEPAC or possible
  !                       future extensions.
  !
  !                         *** Task Terminated ***
  !                   Reported by the value of IDID=-33
  !
  !             IDID = -33 -- The code has encountered trouble from which
  !                       it cannot recover.  A message is printed
  !                       explaining the trouble and control is returned
  !                       to the calling program.  For example, this
  !                       occurs when invalid input is detected.
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
  !             RWORK(12)--If the tolerances have been increased by the
  !                        code (IDID = -2), they were multiplied by the
  !                        value in RWORK(12).
  !
  !             RWORK(20+I)--which contains the approximate derivative
  !                        of the solution component Y(I).  In DDERKF, it
  !                        is always obtained by calling subroutine DF to
  !                        evaluate the differential equation using T and
  !                        Y(*).
  !
  !- *********************************************************************
  !- * INPUT -- What To Do To Continue The Integration **
  !- *             (calls after the first)             **
  !- *********************************************************************
  !
  !        This code is organized so that subsequent calls to continue the
  !        integration involve little (if any) additional effort on your
  !        part.  You must monitor the IDID parameter to determine
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
  !        The parameter INFO(1) is used by the code to indicate the
  !        beginning of a new problem and to indicate whether integration
  !        is to be continued.  You must input the value  INFO(1) = 0
  !        when starting a new problem.  You must input the value
  !        INFO(1) = 1  if you wish to continue after an interrupted task.
  !        Do not set  INFO(1) = 0  on a continuation call unless you
  !        want the code to restart at the current T.
  !
  !                         *** Following a Completed Task ***
  !         If
  !             IDID = 1, call the code again to continue the integration
  !                     another step in the direction of TOUT.
  !
  !             IDID = 2, define a new TOUT and call the code again.
  !                     TOUT must be different from T.  You cannot change
  !                     the direction of integration without restarting.
  !
  !                         *** Following an Interrupted Task ***
  !                     To show the code that you realize the task was
  !                     interrupted and that you want to continue, you
  !                     must take appropriate action and reset INFO(1) = 1
  !         If
  !             IDID = -1, the code has attempted 500 steps.
  !                     If you want to continue, set INFO(1) = 1 and
  !                     call the code again.  An additional 500 steps
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
  !                     inefficient to solve such problems with DDERKF.
  !                     The code DDEBDF in DEPAC handles this task
  !                     efficiently.  If you are absolutely sure you want
  !                     to continue with DDERKF, set INFO(1)=1 and call
  !                     the code again.
  !
  !             IDID = -5, you are using DDERKF very inefficiently by
  !                     choosing output points TOUT so close together that
  !                     the step size is repeatedly forced to be rather
  !                     smaller than necessary.  If you are willing to
  !                     accept solutions at the steps chosen by the code,
  !                     a good way to proceed is to use the intermediate
  !                     output mode (setting INFO(3)=1).  If you must have
  !                     solutions at so many specific TOUT points, the
  !                     code DDEABM in DEPAC handles this task
  !                     efficiently.  If you want to continue with DDERKF,
  !                     set INFO(1)=1 and call the code again.
  !
  !             IDID = -6,-7,..,-32  --- cannot occur with this code but
  !                     used by other members of DEPAC or possible future
  !                     extensions.
  !
  !                         *** Following a Terminated Task ***
  !         If
  !             IDID = -33, you cannot continue the solution of this
  !                     problem.  An attempt to do so will result in your
  !                     run being terminated.
  !
  !- *********************************************************************
  !- Long Description:
  !
  !- *********************************************************************
  !- *             DEPAC Package Overview           **
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
  !               L. F. Shampine and H. A. Watts, Practical solution of
  !                 ordinary differential equations by Runge-Kutta
  !                 methods, Report SAND76-0585, Sandia Laboratories,
  !                 1976.
  !***
  ! **Routines called:**  DRKFS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891024  Changed references from DVNORM to DHVNRM.  (WRB)
  !   891024  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Convert XERRWV calls to XERMSG calls, make Prologue comments
  !           consistent with DERKF.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
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
      !DIR$ ASSUME_ALIGNED Atol:64, Rtol:64, Rwork:64, Y:64
  !
  INTEGER(i4) ::  kdi, kf1, kf2, kf3, kf4, kf5, kh, krer, ktf, kto, ktstar, ku, kyp, kys
  LOGICAL :: stiff, nonstf
  CHARACTER(8) :: xern1
  CHARACTER(16) :: xern3
  !
  !     CHECK FOR AN APPARENT INFINITE LOOP
  !
  !* FIRST EXECUTABLE STATEMENT  DDERKF
  IF( Info(1)==0 ) Iwork(Liw) = 0
  IF( Iwork(Liw)>=5 ) THEN
    IF( T==Rwork(21+Neq) ) THEN
      WRITE (xern3,'(1PE15.6)') T
      ERROR STOP 'DDERKF : AN APPARENT INFINITE LOOP HAS BEEN DETECTED.&
        & YOU HAVE MADE REPEATED CALLS AT  T  AND THE INTEGRATION&
        & HAS NOT ADVANCED. CHECK THE WAY YOU HAVE SET PARAMETERS FOR THE CALL&
        & TO THE CODE, PARTICULARLY INFO(1).'
      RETURN
    END IF
  END IF
  !
  !     CHECK LRW AND LIW FOR SUFFICIENT STORAGE ALLOCATION
  !
  Idid = 0
  IF( Lrw<30+7*Neq ) THEN
    WRITE (xern1,'(I8)') Lrw
    ERROR STOP 'DDERKF : LENGTH OF RWORK ARRAY MUST BE AT LEAST 30 + 7*NEQ.'
    Idid = -33
  END IF
  !
  IF( Liw<34 ) THEN
    WRITE (xern1,'(I8)') Liw
    ERROR STOP 'DDERKF : LENGTH OF IWORK ARRAY MUST BE AT LEAST 34.'
    Idid = -33
  END IF
  !
  !     COMPUTE INDICES FOR THE SPLITTING OF THE RWORK ARRAY
  !
  kh = 11
  ktf = 12
  kyp = 21
  ktstar = kyp + Neq
  kf1 = ktstar + 1
  kf2 = kf1 + Neq
  kf3 = kf2 + Neq
  kf4 = kf3 + Neq
  kf5 = kf4 + Neq
  kys = kf5 + Neq
  kto = kys + Neq
  kdi = kto + 1
  ku = kdi + 1
  krer = ku + 1
  !
  !- *********************************************************************
  !     THIS INTERFACING ROUTINE MERELY RELIEVES THE USER OF A LONG
  !     CALLING LIST VIA THE SPLITTING APART OF TWO WORKING STORAGE
  !     ARRAYS. IF THIS IS NOT COMPATIBLE WITH THE USERS COMPILER,
  !     S/HE MUST USE DRKFS DIRECTLY.
  !- *********************************************************************
  !
  Rwork(ktstar) = T
  IF( Info(1)/=0 ) THEN
    stiff = (Iwork(25)==0)
    nonstf = (Iwork(26)==0)
  END IF
  !
  CALL DRKFS(DF,Neq,T,Y,Tout,Info,Rtol,Atol,Idid,Rwork(kh),Rwork(ktf),&
    Rwork(kyp),Rwork(kf1),Rwork(kf2),Rwork(kf3),Rwork(kf4),&
    Rwork(kf5),Rwork(kys),Rwork(kto),Rwork(kdi),Rwork(ku),&
    Rwork(krer),Iwork(21),Iwork(22),Iwork(23),Iwork(24),stiff,&
    nonstf,Iwork(27),Iwork(28))
  !
  Iwork(25) = 1
  IF( stiff ) Iwork(25) = 0
  Iwork(26) = 1
  IF( nonstf ) Iwork(26) = 0
  !
  IF( Idid/=(-2) ) Iwork(Liw) = Iwork(Liw) + 1
  IF( T/=Rwork(ktstar) ) Iwork(Liw) = 0
  !
END SUBROUTINE DDERKF
