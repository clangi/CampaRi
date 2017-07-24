! -----------------------------------------------------------------------------
! This function aims to reduce the fringe regions by collapsing the tree leaves
! on their branches. In this way we have less grouping of fringe regions
!
subroutine check_netcdf_installation(wanting_r_backend)
  use gutenberg
  use m_variables_gen
  implicit none
  logical, intent(inout) :: wanting_r_backend
#ifdef LINK_NETCDF
  if(wanting_r_backend) then
    call sl()
    call spr('------------------------------------------------------------------')
    call spr('ATTENTION: Even if CampaRi was installed using netcdf support, but')
    call spr('you selected to use the R data management system. Only R option ')
    call spr('will be followed (or both in the case of mst_from_trj).')
    call spr('------------------------------------------------------------------')
  end if
#else
  if(.not.wanting_r_backend) then
    call sl()
    call spr('------------------------------------------------------------------')
    call spr('ATTENTION: Even if CampaRi was installed without the netcdf support')
    call spr('the user tried to use the netcdf dumping functionality. This run ')
    call spr('will follow the usual flow without netcdf. If you want to use it ')
    call spr('install the full version of the package. R will be tried as backend.')
    call spr('------------------------------------------------------------------')
    call spr('Please be sure about the modality of analysis')
    wanting_r_backend = .true.
  end if
#endif
end
