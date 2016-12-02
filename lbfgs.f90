#define LP 6
#define XTOL 1.0e-17
#define DEBUG .FALSE.
#define UNITDEBUG 101

!     ----------------------------------------------------------------------
!     This part of file contains the LBFGS algorithm and supporting routines
!     ----------------------------------------------------------------------
      SUBROUTINE LBFGS(N,M,X,F,G,DIAG,EPS,W,IFLAG,  &
                       GTOL,STPMIN,STPMAX,STP,ITER, &
                       INFO, NFEV,                  &
                       LINE_DGINIT,LINE_FINIT,      &
                       LINE_STX,LINE_FX,LINE_DGX,   &
                       LINE_STY,LINE_FY,LINE_DGY,   &
                       LINE_STMIN,LINE_STMAX,       &
                       LINE_BRACKT,LINE_STAGE1,LINE_INFOC) BIND(C)
      use iso_c_binding
      implicit none
      real(C_DOUBLE),intent(inout) :: GTOL
      real(C_DOUBLE),intent(in),value :: STPMIN,STPMAX
      real(C_DOUBLE),intent(inout) :: STP
      integer(C_INT),intent(inout) :: ITER,IFLAG,INFO,NFEV
      integer(C_INT),intent(in),value  :: N,M
      real(C_DOUBLE),intent(in),value  :: F,EPS
      real(C_DOUBLE),intent(inout) :: X(N),DIAG(N),W(N*(2*M+1)+2*M)
      real(C_DOUBLE),intent(in)    :: G(N)
      real(C_DOUBLE),intent(inout) :: LINE_DGINIT,LINE_FINIT
      logical(C_BOOL),intent(inout):: LINE_BRACKT,LINE_STAGE1
      integer(C_INT),intent(inout) :: LINE_INFOC
      real(C_DOUBLE),intent(inout) :: LINE_STX,LINE_FX,LINE_DGX
      real(C_DOUBLE),intent(inout) :: LINE_STY,LINE_FY,LINE_DGY
      real(C_DOUBLE),intent(inout) :: LINE_STMIN,LINE_STMAX
!
!        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!                          JORGE NOCEDAL
!                        *** July 1990 ***
!
! 
!     This subroutine solves the unconstrained minimization problem
! 
!                      min F(x),    x= (x1,x2,...,xN),
!
!      using the limited memory BFGS method. The routine is especially
!      effective on problems involving a large number of variables. In
!      a typical iteration of this method an approximation Hk to the
!      inverse of the Hessian is obtained by applying M BFGS updates to
!      a diagonal matrix Hk0, using information from the previous M steps.
!      The user specifies the number M, which determines the amount of
!      storage required by the routine. The user may also provide the
!      diagonal matrices Hk0 if not satisfied with the default choice.
!      The algorithm is described in "On the limited memory BFGS method
!      for large scale optimization", by D. Liu and J. Nocedal,
!      Mathematical Programming B 45 (1989) 503-528.
! 
!      The user is required to calculate the function value F and its
!      gradient G. In order to allow the user complete control over
!      these computations, reverse  communication is used. The routine
!      must be called repeatedly under the control of the parameter
!      IFLAG. 
!
!      The steplength is determined at each iteration by means of the
!      line search routine MCVSRCH, which is a slight modification of
!      the routine CSRCH written by More' and Thuente.
! 
!      The calling statement is 
! 
!          CALL LBFGS(N,M,X,F,G,DIAG,EPS,XTOL,W,IFLAG)
! 
!      where
! 
!     N       is an integer(C_INT) variable that must be set by the user to the
!             number of variables. It is not altered by the routine.
!             Restriction: N>0.
! 
!     M       is an integer(C_INT) variable that must be set by the user to
!             the number of corrections used in the BFGS update. It
!             is not altered by the routine. Values of M less than 3 are
!             not recommended; large values of M will result in excessive
!             computing time. 3<= M <=7 is recommended. Restriction: M>0.
! 
!     X       is a real(C_DOUBLE) array of length N. On initial entry
!             it must be set by the user to the values of the initial
!             estimate of the solution vector. On exit with IFLAG=0, it
!             contains the values of the variables at the best point
!             found (usually a solution).
! 
!     F       is a real(C_DOUBLE) variable. Before initial entry and on
!             a re-entry with IFLAG=1, it must be set by the user to
!             contain the value of the function F at the point X.
! 
!     G       is a real(C_DOUBLE) array of length N. Before initial
!             entry and on a re-entry with IFLAG=1, it must be set by
!             the user to contain the components of the gradient G at
!             the point X.
! 
!     DIAG    is a real(C_DOUBLE) array of length N. 
!             Restriction: all elements of DIAG must be positive.
! 
! 
!     EPS     is a positive real(C_DOUBLE) variable that must be set by
!             the user, and determines the accuracy with which the solution
!             is to be found. The subroutine terminates when
!
!                         ||G|| < EPS max(1,||X||),
!
!             where ||.|| denotes the Euclidean norm.
! 
!     XTOL    is a  positive real(C_DOUBLE) variable that must be set by
!             the user to an estimate of the machine precision (e.g.
!             10**(-16) on a SUN station 3/60). The line search routine will
!             terminate if the relative width of the interval of uncertainty
!             is less than XTOL.
! 
!     W       is a real(C_DOUBLE) array of length N(2M+1)+2M used as
!             workspace for LBFGS. This array must not be altered by the
!             user.
! 
!     IFLAG   is an integer(C_INT) variable that must be set to 0 on initial entry
!             to the subroutine. A return with IFLAG<0 indicates an error,
!             and IFLAG=0 indicates that the routine has terminated without
!             detecting errors. On a return with IFLAG=1, the user must
!             evaluate the function F and gradient G. On a return with
!             IFLAG=2, the user must provide the diagonal matrix Hk0.
! 
!             The following negative values of IFLAG, detecting an error,
!             are possible:
! 
!              IFLAG=-1  The line search routine MCSRCH failed. The
!                        parameter INFO provides more detailed information
!                        (see also the documentation of MCSRCH):
!
!                       INFO = 0  IMPROPER INPUT PARAMETERS.
!
!                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
!                                 UNCERTAINTY IS AT MOST XTOL.
!
!                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
!                                 REQUIRED AT THE PRESENT ITERATION.
!
!                       INFO = 4  THE STEP IS TOO SMALL.
!
!                       INFO = 5  THE STEP IS TOO LARGE.
!
!                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. 
!                                 THERE MAY NOT BE A STEP WHICH SATISFIES
!                                 THE SUFFICIENT DECREASE AND CURVATURE
!                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL.
!
! 
!              IFLAG=-2  The i-th diagonal element of the diagonal inverse
!                        Hessian approximation, given in DIAG, is not
!                        positive.
! 
!
!
! 
!    COMMON:
! 
!     The subroutine contains one common area, which the user may wish to
!    reference:
! 
! 
!    MP  is an integer(C_INT) variable with default value 6. It is used as the
!        unit number for the printing of the monitoring information
!        controlled by IPRINT.
! 
!    LP  is an integer(C_INT) variable with default value 6. It is used as the
!        unit number for the printing of error messages. This printing
!        may be suppressed by setting LP to a non-positive value.
! 
!    GTOL is a real(C_DOUBLE) variable with default value 0.9, which
!        controls the accuracy of the line search routine MCSRCH. If the
!        function and gradient evaluations are inexpensive with respect
!        to the cost of the iteration (which is sometimes the case when
!        solving very large problems) it may be advantageous to set GTOL
!        to a small value. A typical small value is 0.1.  Restriction:
!        GTOL should be greater than 1.D-04.
! 
!    STPMIN and STPMAX are non-negative real(C_DOUBLE) variables which
!        specify lower and uper bounds for the step in the line search.
!        Their default values are 1.D-20 and 1.D+20, respectively. These
!        values need not be modified unless the exponents are too large
!        for the machine being used, or unless the problem is extremely
!        badly scaled (in which case the exponents should be increased).
! 
!
!  MACHINE DEPENDENCIES
!
!        The only variables that are machine-dependent are XTOL,
!        STPMIN and STPMAX.
! 
!
!    Other routines called directly:   MCSRCH
! 
!    Input/Output  :  No input; diagnostic messages on unit MP and
!                     error messages on unit LP.
! 
! 
!
!     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
!     ---------------------------------------
!     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
!         OTHER TEMPORARY DATA
!     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
!     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
!         IN THE FORMULA THAT COMPUTES H*G.
!     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
!         STEPS.
!     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
!         GRADIENT DIFFERENCES.
!
!     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
!     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
!
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      real(C_DOUBLE),parameter :: ONE = 1.0_C_DOUBLE
      real(C_DOUBLE),parameter :: ZERO = 0.0_C_DOUBLE

      real(C_DOUBLE) :: GNORM,FTOL,YS,YY,SQ,YR,BETA,XNORM
      integer(C_INT) :: POINT,ISPT,IYPT,MAXFEV, &
              BOUND,NPT,CP,I,INMC,IYCN,ISCN
!
!
!     INITIALIZE
!     ----------
!

      if( DEBUG ) then
        write(*,*) '================================================'
        write(*,*) '== ITER,IFLAG',ITER,IFLAG
        write(*,*) 
        write(*,*) '== X'
        write(*,*) X(:)
        write(*,*) '== G'
        write(*,*) G(:)
        write(*,*) '== DIAG'
        write(*,*) DIAG(:)
        write(*,*) '== W'
        write(*,*) W(:)
        write(*,*) 
        write(*,*) 
      endif

!     PARAMETERS FOR LINE SEARCH ROUTINE
      FTOL = 1.0D-4
      MAXFEV = 20

      ISPT = N + 2 * M
      IYPT = ISPT + N * M     
      POINT = MAX(0,MOD(ITER-1,M))
      NPT = POINT * N

      GO TO (10,172) IFLAG+1

     
  10  CONTINUE

      W(ISPT+1:ISPT+N) = -G(1:N) * DIAG(1:N)
      GNORM = NORM2(G(:))
!
!
!    --------------------
!     MAIN ITERATION LOOP
!    --------------------
!
 80   ITER  = ITER + 1
      INFO  = 0
      BOUND = MIN( ITER-1 , M)

      IF( ITER /= 1 ) THEN

        YS = DOT_PRODUCT( W(IYPT+NPT+1:IYPT+NPT+N) , W(ISPT+NPT+1:ISPT+NPT+N) )
        YY = DOT_PRODUCT( W(IYPT+NPT+1:IYPT+NPT+N) , W(IYPT+NPT+1:IYPT+NPT+N) )
        DIAG(1:N)= YS / YY

!
!       COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!       "Updating quasi-Newton matrices with limited storage",
!       Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!       ---------------------------------------------------------
!
        POINT = MODULO(ITER - 1,M)
        CP = POINT
        IF (POINT == 0) CP = M
        W(N+CP)= ONE/YS
        W(1:N)= -G(1:N)

        CP = POINT
        DO I= 1,BOUND
           CP = CP - 1
           IF (CP ==  -1) CP = M - 1
           SQ= DOT_PRODUCT(W(ISPT+CP*N+1:ISPT+CP*N+N),W(1:N))
           INMC=N+M+CP+1
           IYCN=IYPT+CP*N
           W(INMC)= W(N+CP+1) * SQ
           W(1:N) = W(1:N) - W(INMC) * W(IYCN+1:IYCN+N)
        ENDDO     

        W(1:N) = DIAG(1:N) * W(1:N)

        DO I=1,BOUND
           YR = DOT_PRODUCT(W(IYPT+CP*N+1:IYPT+CP*N+N),W(1:N))
           BETA = W(N+CP+1) * YR
           INMC = N + M + CP + 1
           BETA = W(INMC) - BETA
           ISCN = ISPT + CP * N
           W(1:N) = W(1:N) + BETA * W(ISCN+1:ISCN+N)
           CP = CP + 1
           IF (CP == M) CP = 0
        ENDDO
   
!
!       STORE THE NEW SEARCH DIRECTION
        W(ISPT+POINT*N+1:ISPT+POINT*N+N) = W(1:N)

      ENDIF
!
!     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION 
!     BY USING THE LINE SEARCH ROUTINE MCSRCH
!     ----------------------------------------------------
      NFEV = 0
      STP = ONE
      W(1:N) = G(1:N)

 172  CONTINUE
      if( DEBUG ) then
        write(UNITDEBUG,*) 'INFO BEFORE MCSRCH',INFO
      endif
      CALL MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,MAXFEV,INFO,NFEV, &
                  DIAG,GTOL,STPMIN,STPMAX,LINE_DGINIT,LINE_FINIT, &
                  LINE_STX,LINE_FX,LINE_DGX, &
                  LINE_STY,LINE_FY,LINE_DGY, &
                  LINE_STMIN,LINE_STMAX, &
                  LINE_BRACKT,LINE_STAGE1,LINE_INFOC)
      if( DEBUG ) then
        write(UNITDEBUG,*) 'X:',X(:)
      endif

      IF (INFO  ==  -1) THEN
        IFLAG = 1
        RETURN
      ENDIF
      IF (INFO /= 1) THEN
        IFLAG = -1
        WRITE(LP,'(/,a,i4,/,a,/,a,i4)') ' IFLAG= ',IFLAG,  &
                                   ' LINE SEARCH FAILED ', &
                                   ' ERROR RETURNED FROM LINE SEARCH ',INFO
        RETURN
      ENDIF
!
!     COMPUTE THE NEW STEP AND GRADIENT CHANGE 
!     -----------------------------------------
!
      NPT = POINT * N
      W(ISPT+NPT+1:ISPT+NPT+N) = STP * W(ISPT+NPT+1:ISPT+NPT+N)
      W(IYPT+NPT+1:IYPT+NPT+N) = G(1:N) - W(1:N)
!
!     TERMINATION TEST
!     ----------------
!
      GNORM = NORM2(G)
      XNORM = NORM2(X)
      XNORM = MAX(1.0D0,XNORM)

      IF ( GNORM/XNORM < EPS ) THEN
        IFLAG = 0
        RETURN
      ENDIF
      GO TO 80
!
!     ------------------------------------------------------------
!     END OF MAIN ITERATION LOOP. ERROR EXITS.
!     ------------------------------------------------------------
!

      END SUBROUTINE LBFGS


!
!
!     **************************
!     LINE SEARCH ROUTINE MCSRCH
!     **************************
!
      SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,MAXFEV,INFO,NFEV,WA, &
                        GTOL,STPMIN,STPMAX,DGINIT,FINIT, &
                        STX,FX,DGX,STY,FY,DGY,STMIN,STMAX, &
                        BRACKT,STAGE1,INFOC)
      use iso_c_binding
      implicit none
      real(C_DOUBLE),intent(in)     :: GTOL,STPMIN,STPMAX
      integer(C_INT),intent(in)     :: N,MAXFEV
      integer,intent(inout)         :: INFO,NFEV
      real(C_DOUBLE),intent(in)     :: F,FTOL
      real(C_DOUBLE),intent(inout)  :: STP,DGINIT,FINIT
      real(C_DOUBLE),intent(in)     :: G(N)
      real(C_DOUBLE),intent(inout)  :: X(N),S(N),WA(N)
      logical(C_BOOL),intent(inout) :: BRACKT,STAGE1
      integer(C_INT),intent(inout)  :: INFOC
      real(C_DOUBLE),intent(inout)  :: STX,FX,DGX
      real(C_DOUBLE),intent(inout)  :: STY,FY,DGY
      real(C_DOUBLE),intent(inout)  :: STMIN,STMAX
!
!                     SUBROUTINE MCSRCH
!                
!     A slight modification of the subroutine CSRCH of More' and Thuente.
!     The changes are to allow reverse communication, and do not affect
!     the performance of the routine. 
!
!     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
!     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
!
!     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
!     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
!     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
!     MINIMIZER OF THE MODIFIED FUNCTION
!
!          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
!
!     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
!     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
!     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
!     CONTAINS A MINIMIZER OF F(X+STP*S).
!
!     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
!     THE SUFFICIENT DECREASE CONDITION
!
!           F(X+STP*S) <= F(X) + FTOL*STP*(GRADF(X)'S),
!
!     AND THE CURVATURE CONDITION
!
!           ABS(GRADF(X+STP*S)'S)) <= GTOL*ABS(GRADF(X)'S).
!
!     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
!     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
!     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
!     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
!     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
!     SATISFIES THE SUFFICIENT DECREASE CONDITION.
!
!     THE SUBROUTINE STATEMENT IS
!
!        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
!     WHERE
!
!       N IS A POSITIVE integer(C_INT) INPUT VARIABLE SET TO THE NUMBER
!         OF VARIABLES.
!
!       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
!         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
!         X + STP*S.
!
!       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
!         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
!
!       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
!         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
!         OF F AT X + STP*S.
!
!       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
!         SEARCH DIRECTION.
!
!       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
!         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
!         STP CONTAINS THE FINAL ESTIMATE.
!
!       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
!         communication implementation GTOL is defined in a COMMON
!         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
!         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
!         SATISFIED.
!
!       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
!         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
!         IS AT MOST XTOL.
!
!       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
!         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
!         communication implementatin they are defined in a COMMON
!         statement).
!
!       MAXFEV IS A POSITIVE integer(C_INT) INPUT VARIABLE. TERMINATION
!         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
!         MAXFEV BY THE END OF AN ITERATION.
!
!       INFO IS AN integer(C_INT) OUTPUT VARIABLE SET AS FOLLOWS:
!
!         INFO = 0  IMPROPER INPUT PARAMETERS.
!
!         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
!
!         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
!                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
!
!         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
!                   IS AT MOST XTOL.
!
!         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
!
!         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
!
!         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
!
!         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
!                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
!                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
!                   TOLERANCES MAY BE TOO SMALL.
!
!       NFEV IS AN integer(C_INT) OUTPUT VARIABLE SET TO THE NUMBER OF
!         CALLS TO FCN.
!
!       WA IS A WORK ARRAY OF LENGTH N.
!
!     SUBPROGRAMS CALLED
!
!       MCSTEP
!
!       FORTRAN-SUPPLIED...ABS,MAX,MIN
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!     JORGE J. MORE', DAVID J. THUENTE
!
!     **********
      integer(C_INT) :: J
      real(C_DOUBLE) :: DG,DGM,DGTEST,DGXM,DGYM, &
             FTEST1,FM,FXM,FYM,WIDTH,WIDTH1
      real(C_DOUBLE),parameter :: P5     = 0.50_C_DOUBLE
      real(C_DOUBLE),parameter :: P66    = 0.66_C_DOUBLE
      real(C_DOUBLE),parameter :: XTRAPF = 4.00_C_DOUBLE
      real(C_DOUBLE),parameter :: ZERO   = 0.00_C_DOUBLE

      if( DEBUG ) then
        write(UNITDEBUG,*) 'INFO INFOC',INFO,INFOC
        write(UNITDEBUG,*) 'STP       ',STP
        write(UNITDEBUG,*) 'F,FTOL:',F,FTOL
        write(UNITDEBUG,*) 'X:',X(:)
        write(UNITDEBUG,*) 'G:',G(:)
        write(UNITDEBUG,*) 'S:',S(:)
        write(UNITDEBUG,*) 'WA:',WA(:)
      endif

      DGTEST = FTOL * DGINIT
      WIDTH = STPMAX - STPMIN
      WIDTH1 = WIDTH / P5

      IF(INFO == -1) GO TO 45
      INFOC = 1
!
!     CHECK THE INPUT PARAMETERS FOR ERRORS.
!
      IF ( STP <= ZERO .OR. FTOL < ZERO .OR.  &
          GTOL < ZERO .OR. XTOL < ZERO .OR. STPMIN < ZERO &
          .OR. STPMAX < STPMIN ) RETURN
!
!     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
!     AND CHECK THAT S IS A DESCENT DIRECTION.
!
      DGINIT = DOT_PRODUCT( G , S )

      IF (DGINIT > ZERO) then
        WRITE(LP,'(a)') ' THE SEARCH DIRECTION IS NOT A DESCENT'
        RETURN
      ENDIF
!
!     INITIALIZE LOCAL VARIABLES.
!


      BRACKT = .FALSE.
      STAGE1 = .TRUE.
      NFEV = 0
      FINIT = F
      DGTEST = FTOL * DGINIT
      WA(:) = X(:)


!
!     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
!     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
!     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
!     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
!     THE INTERVAL OF UNCERTAINTY.
!     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
!     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
!
      STX = ZERO
      FX = FINIT
      DGX = DGINIT
      STY = ZERO
      FY = FINIT
      DGY = DGINIT
!
!     START OF ITERATION.
!
   30 CONTINUE
!
!     SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
!     TO THE PRESENT INTERVAL OF UNCERTAINTY.
!
      IF (BRACKT) THEN
         STMIN = MIN(STX,STY)
         STMAX = MAX(STX,STY)
      ELSE
         STMIN = STX
         STMAX = STP + XTRAPF*(STP - STX)
      END IF
!
!     FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
!
      STP = MAX(STPMIN,STP)
      STP = MIN(STP,STPMAX)
!
!     IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
!     STP BE THE LOWEST POINT OBTAINED SO FAR.
!
      IF ((BRACKT .AND. (STP <= STMIN .OR. STP >= STMAX)) &
         .OR. NFEV >= MAXFEV-1 .OR. INFOC  ==  0 &
         .OR. (BRACKT .AND. STMAX-STMIN <= XTOL*STMAX)) STP = STX

!
!     EVALUATE THE FUNCTION AND GRADIENT AT STP
!     AND COMPUTE THE DIRECTIONAL DERIVATIVE.
!     We return to main program to obtain F and G.
!
      X(:) = WA(:) + STP * S(:)
      INFO = -1

      RETURN
!
   45 INFO = 0
      NFEV = NFEV + 1
      DG = SUM( G(:) * S(:) )
      FTEST1 = FINIT + STP * DGTEST
!
!     TEST FOR CONVERGENCE.
!
      IF ((BRACKT .AND. (STP <= STMIN .OR. STP >= STMAX)) &
         .OR. INFOC  ==  0) INFO = 6
      IF (STP  ==  STPMAX .AND. &
          F <= FTEST1 .AND. DG <= DGTEST) INFO = 5
      IF (STP  ==  STPMIN .AND.  &
          (F > FTEST1 .OR. DG >= DGTEST)) INFO = 4
      IF (NFEV >= MAXFEV) INFO = 3
      IF (BRACKT .AND. STMAX-STMIN <= XTOL*STMAX) INFO = 2
      IF (F <= FTEST1 .AND. ABS(DG) <= GTOL*(-DGINIT)) INFO = 1
!
!     CHECK FOR TERMINATION.
!
      IF (INFO /= 0) RETURN
!
!     IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
!     FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
!
      IF (STAGE1 .AND. F <= FTEST1 .AND. &
          DG >= MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE.
!
!     A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
!     WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
!     FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
!     DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
!     OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
!
      IF (STAGE1 .AND. F <= FX .AND. F > FTEST1) THEN
!
!        DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
!
         FM = F - STP * DGTEST
         FXM = FX - STX * DGTEST
         FYM = FY - STY * DGTEST
         DGM = DG - DGTEST
         DGXM = DGX - DGTEST
         DGYM = DGY - DGTEST
!
!        CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!        AND TO COMPUTE THE NEW STEP.
!
         CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,BRACKT,STMIN,STMAX,INFOC)
!
!        RESET THE FUNCTION AND GRADIENT VALUES FOR F.
!
         FX = FXM + STX * DGTEST
         FY = FYM + STY * DGTEST
         DGX = DGXM + DGTEST
         DGY = DGYM + DGTEST
      ELSE
!
!        CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!        AND TO COMPUTE THE NEW STEP.
!
         CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG, BRACKT,STMIN,STMAX,INFOC)
      END IF
!
!     FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
!     INTERVAL OF UNCERTAINTY.
!
      IF (BRACKT) THEN
         IF (ABS(STY-STX) >= P66 * WIDTH1) STP = STX + P5 * (STY - STX)
         WIDTH1 = WIDTH
         WIDTH = ABS(STY-STX)
      END IF
!
!     END OF ITERATION.
!
      GO TO 30

      END SUBROUTINE MCSRCH


      SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,STPMIN,STPMAX,INFO)
      use iso_c_binding
      implicit none
      integer(C_INT) INFO
      real(C_DOUBLE) STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
      logical(C_BOOL),intent(inout) :: BRACKT
      logical(C_BOOL) BOUND
!
!     SUBROUTINE MCSTEP
!
!     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
!     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
!     A MINIMIZER OF THE FUNCTION.
!
!     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
!     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
!     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
!     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
!     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
!     WITH ENDPOINTS STX AND STY.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
!                        STPMIN,STPMAX,INFO)
!
!     WHERE
!
!       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
!         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
!         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
!         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
!         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
!
!       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
!         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
!         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
!         UPDATED APPROPRIATELY.
!
!       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
!         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
!         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
!         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
!
!       BRACKT IS A logical(C_BOOL) VARIABLE WHICH SPECIFIES IF A MINIMIZER
!         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
!         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
!         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
!
!       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
!         AND UPPER BOUNDS FOR THE STEP.
!
!       INFO IS AN integer(C_INT) OUTPUT VARIABLE SET AS FOLLOWS:
!         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
!         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
!         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!     JORGE J. MORE', DAVID J. THUENTE
!
      real(C_DOUBLE) GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA

      INFO = 0
!
!     CHECK THE INPUT PARAMETERS FOR ERRORS.
!
      IF ((BRACKT .AND. (STP <= MIN(STX,STY) .OR. &
          STP >= MAX(STX,STY))) .OR.  &
          DX*(STP-STX) >= 0.0 .OR. STPMAX < STPMIN) RETURN
!
!     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
!
      SGND = DP*(DX/ABS(DX))
!
!     FIRST CASE. A HIGHER FUNCTION VALUE.
!     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
!     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
!     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
!
      IF (FP > FX) THEN
         INFO = 1
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP < STX) GAMMA = -GAMMA
         P = (GAMMA - DX) + THETA
         Q = ((GAMMA - DX) + GAMMA) + DP
         R = P/Q
         STPC = STX + R*(STP - STX)
         STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
         IF (ABS(STPC-STX) < ABS(STPQ-STX)) THEN
            STPF = STPC
         ELSE
           STPF = STPC + (STPQ - STPC)/2
         END IF
         BRACKT = .TRUE.
!
!     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
!     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
!     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
!     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
!
      ELSE IF (SGND < 0.0) THEN
         INFO = 2
         BOUND = .FALSE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP > STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = ((GAMMA - DP) + GAMMA) + DX
         R = P/Q
         STPC = STP + R*(STX - STP)
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (ABS(STPC-STP) > ABS(STPQ-STP)) THEN
            STPF = STPC
         ELSE
            STPF = STPQ
         END IF
         BRACKT = .TRUE.
!
!     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
!     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
!     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
!     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
!     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
!     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
!     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
!     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
!
      ELSE IF (ABS(DP) < ABS(DX)) THEN
         INFO = 3
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
!
!        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
!        TO INFINITY IN THE DIRECTION OF THE STEP.
!
         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
         IF (STP > STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = (GAMMA + (DX - DP)) + GAMMA
         R = P/Q
         IF (R < 0.0 .AND. GAMMA .NE. 0.0) THEN
            STPC = STP + R*(STX - STP)
         ELSE IF (STP > STX) THEN
            STPC = STPMAX
         ELSE
            STPC = STPMIN
         END IF
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (BRACKT) THEN
            IF (ABS(STP-STPC) < ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
            END IF
         ELSE
            IF (ABS(STP-STPC) > ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
            END IF
         END IF
!
!     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
!     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
!     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
!     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
!
      ELSE
         INFO = 4
         BOUND = .FALSE.
         IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP > STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
         ELSE IF (STP > STX) THEN
            STPF = STPMAX
         ELSE
            STPF = STPMIN
         END IF
      END IF
!
!     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
!     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
!
      IF (FP > FX) THEN
         STY = STP
         FY = FP
         DY = DP
      ELSE
         IF (SGND < 0.0) THEN
            STY = STX
            FY = FX
            DY = DX
         END IF
         STX = STP
         FX = FP
         DX = DP
      END IF
!
!     COMPUTE THE NEW STEP AND SAFEGUARD IT.
!
      STPF = MIN(STPMAX,STPF)
      STPF = MAX(STPMIN,STPF)
      STP = STPF
      IF (BRACKT .AND. BOUND) THEN
         IF (STY > STX) THEN
            STP = MIN(STX+0.66*(STY-STX),STP)
         ELSE
            STP = MAX(STX+0.66*(STY-STX),STP)
         END IF
      END IF


      END SUBROUTINE MCSTEP

