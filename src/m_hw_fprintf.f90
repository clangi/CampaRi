module m_hw_fprintf
  interface
    function handwrapped_fprint(thename, &
      i1,i3,setis1,cutv1,distv1,ivec21,invvec1,cutv12, &
      vv11,vv21,vv31,vv12,vv22,vv32,vv13,vv23,vv33,ii1,jj1,kk1,&
      print_newline) bind(C, name='fprintf_wrapper')
      use, INTRINSIC:: ISO_C_BINDING
      implicit none
      integer(c_int) :: handwrapped_fprint, print_newline
      character(c_char) :: thename(*)
      integer(c_int) :: i1,i3,setis1,cutv1,ivec21,&
      invvec1,cutv12,vv11,vv21,vv31,vv12,vv22,vv32,vv13,vv23,vv33,ii1,jj1,kk1
      real(c_double) :: distv1
    end function handwrapped_fprint
  end interface
end module m_hw_fprintf
