

 interface

   type(C_PTR) function libmin_init(ndim,nhess) BIND(C)
     use,intrinsic :: iso_c_binding
     implicit none
     integer(C_INT),intent(in),value :: ndim
     integer(C_INT),intent(in),value :: nhess
   end function libmin_init

   type(C_PTR) function libmin_init_diag(ndim,nhess,diag) BIND(C)
     use,intrinsic :: iso_c_binding
     implicit none
     integer(C_INT),intent(in),value :: ndim
     integer(C_INT),intent(in),value :: nhess
     real(C_DOUBLE),intent(in)       :: diag(*)
   end function libmin_init_diag

   integer(C_INT) function libmin_execute(plan,xx,ff,gradf) BIND(C)
     use,intrinsic :: iso_c_binding
     implicit none
     type(C_PTR), value           :: plan
     real(C_DOUBLE),intent(inout) :: xx(*)
     real(C_DOUBLE), value        :: ff
     real(C_DOUBLE),intent(in)    :: gradf(*)
   end function libmin_execute

   subroutine libmin_destroy(plan) BIND(C)
     use,intrinsic :: iso_c_binding
     implicit none
     type(C_PTR), value        :: plan
   end subroutine libmin_destroy

 end interface

