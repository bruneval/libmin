#define XTOL 1.0e-17
#define DEBUG .FALSE.
#define UNITDEBUG 101

!     ----------------------------------------------------------------------
!     This part of file contains the LBFGS algorithm and supporting routines
!     ----------------------------------------------------------------------
      subroutine LBFGS(N,M,X,F,G,DIAG,EPS,W,IFLAG,  &
                       GTOL,STPMIN,STPMAX,STP,ITER, &
                       INFO, NFEV,                  &
                       LINE_DGINIT,LINE_FINIT,      &
                       LINE_STX,LINE_FX,LINE_DGX,   &
                       LINE_STY,LINE_FY,LINE_DGY,   &
                       LINE_STMIN,LINE_STMAX,       &
                       LINE_BRACKT,LINE_STAGE1,LINE_INFOC) BIND(C)
      use, intrinsic :: iso_c_binding
      use, intrinsic :: iso_fortran_env, only : ERROR_UNIT
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
      POINT = MAX( 0 , MOD(ITER-1,M) )
      NPT = POINT * N
      ITER  = ITER + 1
      BOUND = MIN( ITER-1 , M)
      GNORM = NORM2(G(:))
!
!     TERMINATION TEST
!     ----------------
!
      XNORM = NORM2(X)
      XNORM = MAX(1.0D0,XNORM)
      IF ( GNORM / XNORM < EPS ) THEN
        IFLAG = 0
        RETURN
      ENDIF


      !
      ! Entering the subroutine with a new position and gradient
      ! or entering for the first time ever
      if( IFLAG == 1 ) then
        call MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,MAXFEV,INFO,NFEV, &
                    DIAG,GTOL,STPMIN,STPMAX,LINE_DGINIT,LINE_FINIT, &
                    LINE_STX,LINE_FX,LINE_DGX, &
                    LINE_STY,LINE_FY,LINE_DGY, &
                    LINE_STMIN,LINE_STMAX, &
                    LINE_BRACKT,LINE_STAGE1,LINE_INFOC)
!
!       COMPUTE THE NEW STEP AND GRADIENT CHANGE 
!
        NPT = POINT * N
        W(ISPT+NPT+1:ISPT+NPT+N) = STP * W(ISPT+NPT+1:ISPT+NPT+N)
        W(IYPT+NPT+1:IYPT+NPT+N) = G(1:N) - W(1:N)

      else
        W(ISPT+1:ISPT+N) = -G(1:N) * DIAG(1:N)
      endif

      INFO  = 0


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


      if( DEBUG ) then
        write(UNITDEBUG,*) 'INFO before second call to MCSRCH',INFO
      endif
      CALL MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,MAXFEV,INFO,NFEV, &
                  DIAG,GTOL,STPMIN,STPMAX,LINE_DGINIT,LINE_FINIT, &
                  LINE_STX,LINE_FX,LINE_DGX, &
                  LINE_STY,LINE_FY,LINE_DGY, &
                  LINE_STMIN,LINE_STMAX, &
                  LINE_BRACKT,LINE_STAGE1,LINE_INFOC)
      if( DEBUG ) then
        write(UNITDEBUG,*) 'INFO after  second call to MCSRCH',INFO
      endif

      IF (INFO  ==  -1) THEN
        IFLAG = 1
        RETURN
      else
        IFLAG = -1
        WRITE(ERROR_UNIT,'(/,a,i4,/,a,/,a,i4)') ' IFLAG= ',IFLAG,  &
                                   ' LINE SEARCH FAILED ', &
                                   ' ERROR RETURNED FROM LINE SEARCH ',INFO
        RETURN
      ENDIF


!
!     ------------------------------------------------------------
!     END OF MAIN ITERATION LOOP. ERROR EXITS.
!     ------------------------------------------------------------
!

      END SUBROUTINE LBFGS

