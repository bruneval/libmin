!=====================
!
!=====================
program main_test
 use,intrinsic :: iso_c_binding
 implicit none

 type(C_PTR)    :: plan
 integer(C_INT) :: ndim
 integer(C_INT) :: nhist
 real(C_DOUBLE) :: tolerance
 real(C_DOUBLE),allocatable :: xx(:)
 real(C_DOUBLE),allocatable :: gradf(:)
 real(C_DOUBLE),allocatable :: diag(:)
 real(C_DOUBLE)             :: ff
 integer(C_INT)             :: info

 integer :: idim
 integer :: iter
 integer,parameter :: niter=20

#include "libmin.f03"
!=====================

! ndim  = 1
 ndim  = 3
 nhist = 5
 tolerance = 1.0e-7
 allocate(xx(ndim),gradf(ndim))

! plan = libmin_init(ndim,nhist,tolerance)
 allocate(diag(ndim))
 diag(1) = 1.0 / 1.00 
 diag(2) = 1.0 / 0.25 
 diag(3) = 1.0 / 6.00 
 plan = libmin_init_diag(ndim,nhist,tolerance,diag)
 deallocate(diag)

 do idim=1,ndim
   xx(idim)    =  1.11258d0 + 1.5415 * idim
 enddo


 do iter=1,niter
   ff    = eval_function(ndim,xx)
   gradf = eval_gradient(ndim,xx)

   write(*,*) ' ===== iteration:',iter,' Function value: ',ff
   write(*,*) 'Input x:    ',xx(:)
   write(*,*) 'Input grad: ',gradf(:)

   info = libmin_execute(plan,xx,ff,gradf)

   write(*,*) 'Output x:   ',xx(:)

   if( info < 0 ) then
     write(*,*) 'Problem happened',info
     stop 'STOP'
   endif

   if( info == 0 ) then
     write(*,*) 'Convergence reached'
     write(*,*) 'Minarg f:',xx(:)
     exit
   endif
 enddo


 call libmin_destroy(plan)

 deallocate(xx,gradf)

contains

function eval_function(ndim,xx) result(ff)
 integer(C_INT),intent(in) :: ndim
 real(C_DOUBLE),intent(in) :: xx(ndim)
 real(C_DOUBLE)            :: ff
!====

 ff = 0.500_C_DOUBLE * ( xx(1) - 1.000_C_DOUBLE )**2  &
    + 0.125_C_DOUBLE * ( xx(2) - 3.000_C_DOUBLE )**2  &
   +  3.000_C_DOUBLE * ( xx(3) + 2.000_C_DOUBLE )**2
! ff = 0.50 * ( xx(1) - 0.500 )**2  + 0.0263 * ( xx(1) - 0.500)**4
 
end function eval_function

function eval_gradient(ndim,xx) result(gradf)
 integer(C_INT),intent(in) :: ndim
 real(C_DOUBLE),intent(in) :: xx(ndim)
 real(C_DOUBLE)            :: gradf(ndim)
!====
 
 gradf(1) = 0.5000_C_DOUBLE * ( xx(1) - 1.000_C_DOUBLE )  * 2.0_C_DOUBLE  
 gradf(2) = 0.1250_C_DOUBLE * ( xx(2) - 3.000_C_DOUBLE )  * 2.0_C_DOUBLE  
 gradf(3) = 3.0000_C_DOUBLE * ( xx(3) + 2.000_C_DOUBLE )  * 2.0_C_DOUBLE  

! gradf(1) = ( xx(1) - 0.500 ) + 4.0 * 0.0263 * ( xx(1) - 0.500)**3
 
end function eval_gradient


end program main_test

