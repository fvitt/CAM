!----------------------------------------------------------------------
! soa_gaex_diag_mod --- history diagnostics of SOA condensation and
!                       evaporation by gas-aerosol exchange, shared by
!                       the modal (MAM) and sectional (CARMA) aerosol models
!
! Fields (kg/kg/s):
!   qcon_gaex, qevap_gaex              total over all soa gas species
!   qcon<grp>_gaex, qevap<grp>_gaex    total over the species of each soa class,
!                                      when there is more than one class
!   qcon_gaexNN, qevap_gaexNN          soa gas species NN
!
! The soa class of a soa gas species is its name without the leading 'SOAG'
! and the trailing volatility bin digits, e.g., SOAGbb3 --> bb
!----------------------------------------------------------------------
module soa_gaex_diag_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols

  implicit none
  private

  public :: soa_gaex_diag_t

  integer, parameter :: grp_len = 16

  type :: soa_gaex_diag_t
     private
     integer :: nsoa = -1                          ! number of soa gas species
     integer :: ngrp = -1                          ! number of soa classes
     integer :: nfld_ = -1                         ! number of fields of each kind (cond or evap)
     integer :: ifld_grp0 = -1                     ! field index offset of the class fields
     integer :: ifld_spc0 = -1                     ! field index offset of the species fields
     integer, allocatable :: grp(:)                ! soa class of each soa gas species
     character(len=grp_len), allocatable :: grpname(:)
   contains
     procedure :: init
     procedure :: nfld
     procedure :: accum
     procedure :: write_fields
  end type soa_gaex_diag_t

contains

!----------------------------------------------------------------------
subroutine init(self, soag_names, history_aerosol, add_default_species)

!
! sets the soa classes from the soa gas species names and defines the history fields
!

   use cam_history, only: addfld, add_default

   class(soa_gaex_diag_t), intent(inout) :: self
   character(len=*), intent(in) :: soag_names(:)     ! soa gas species names
   logical,          intent(in) :: history_aerosol   ! add the fields to the default history output
   logical,          intent(in) :: add_default_species ! also add the species fields to the default output

   character(len=*), parameter :: unit = 'kg/kg/s'

   integer :: g, j
   character(len=grp_len) :: gname
   character(len=2)       :: outsoa
   character(len=64)      :: fieldname
   character(len=128)     :: long_name

   self%nsoa = size(soag_names)

   allocate( self%grp(self%nsoa), self%grpname(self%nsoa) )
   self%ngrp = 0
   do j = 1, self%nsoa
      gname = soa_class(soag_names(j))
      self%grp(j) = 0
      do g = 1, self%ngrp
         if (self%grpname(g) == gname) self%grp(j) = g
      end do
      if (self%grp(j) == 0) then
         self%ngrp = self%ngrp + 1
         self%grpname(self%ngrp) = gname
         self%grp(j) = self%ngrp
      end if
   end do

   ! field 1 is the total, then the class fields (if more than one class), then the species fields
   self%ifld_grp0 = 1
   if (self%ngrp > 1) then
      self%ifld_spc0 = self%ifld_grp0 + self%ngrp
   else
      self%ifld_spc0 = self%ifld_grp0
   end if
   self%nfld_ = self%ifld_spc0 + self%nsoa

   call addfld('qcon_gaex', (/'lev'/), 'A', unit, '3D fields for SOA condensation')
   call addfld('qevap_gaex', (/'lev'/), 'A', unit, '3D fields for SOA evaporation')
   if (history_aerosol) then
      call add_default('qcon_gaex', 1, ' ')
      call add_default('qevap_gaex', 1, ' ')
   end if

   if (self%ngrp > 1) then
      do g = 1, self%ngrp
         fieldname = 'qcon'//trim(self%grpname(g))//'_gaex'
         long_name = '3D fields for '//trim(class_description(self%grpname(g)))//' SOA condensation'
         call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
         if (history_aerosol) call add_default(fieldname, 1, ' ')
         fieldname = 'qevap'//trim(self%grpname(g))//'_gaex'
         long_name = '3D fields for '//trim(class_description(self%grpname(g)))//' SOA evaporation'
         call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
         if (history_aerosol) call add_default(fieldname, 1, ' ')
      end do
   end if

   do j = 1, self%nsoa
      write (outsoa, "(I2.2)") j
      fieldname = 'qcon_gaex'//outsoa
      long_name = '3D fields for SOA condensation of '//trim(soag_names(j))
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
      if (history_aerosol .and. add_default_species) call add_default(fieldname, 1, ' ')
      fieldname = 'qevap_gaex'//outsoa
      long_name = '3D fields for SOA evaporation of '//trim(soag_names(j))
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
      if (history_aerosol .and. add_default_species) call add_default(fieldname, 1, ' ')
   end do

end subroutine init

!----------------------------------------------------------------------
pure integer function nfld(self)

! number of fields of each kind (condensation or evaporation),
! the last dimension of the qcon and qevap arrays

   class(soa_gaex_diag_t), intent(in) :: self

   nfld = self%nfld_

end function nfld

!----------------------------------------------------------------------
subroutine accum(self, i, k, jsoa, dq, qcon, qevap)

!
! adds the soa aerosol tendency (kg/kg/s) from soa gas species jsoa at (i,k)
! to the condensation or evaporation fields
!

   class(soa_gaex_diag_t), intent(in) :: self
   integer,  intent(in)    :: i                  ! column index
   integer,  intent(in)    :: k                  ! level index
   integer,  intent(in)    :: jsoa               ! soa gas species index
   real(r8), intent(in)    :: dq                 ! soa aerosol tendency (kg/kg/s)
   real(r8), intent(inout) :: qcon(:,:,:)        ! condensation fields (pcols,pver,nfld)
   real(r8), intent(inout) :: qevap(:,:,:)       ! evaporation fields (pcols,pver,nfld)

   integer :: ig, is

   ig = self%ifld_grp0 + self%grp(jsoa)
   is = self%ifld_spc0 + jsoa

   if (dq >= 0.0_r8) then
      qcon(i,k,1) = qcon(i,k,1) + dq
      if (self%ngrp > 1) qcon(i,k,ig) = qcon(i,k,ig) + dq
      qcon(i,k,is) = qcon(i,k,is) + dq
   else if (dq < 0.0_r8) then
      qevap(i,k,1) = qevap(i,k,1) + dq
      if (self%ngrp > 1) qevap(i,k,ig) = qevap(i,k,ig) + dq
      qevap(i,k,is) = qevap(i,k,is) + dq
   end if

end subroutine accum

!----------------------------------------------------------------------
subroutine write_fields(self, lchnk, qcon, qevap)

!
! writes the condensation and evaporation fields to the history files
!

   use cam_history, only: outfld

   class(soa_gaex_diag_t), intent(in) :: self
   integer,  intent(in) :: lchnk                 ! chunk identifier
   real(r8), intent(in) :: qcon(:,:,:)           ! condensation fields (pcols,pver,nfld)
   real(r8), intent(in) :: qevap(:,:,:)          ! evaporation fields (pcols,pver,nfld)

   integer :: g, j
   character(len=2)  :: outsoa
   character(len=64) :: fieldname

   call outfld('qcon_gaex', qcon(:,:,1), pcols, lchnk)
   call outfld('qevap_gaex', qevap(:,:,1), pcols, lchnk)

   if (self%ngrp > 1) then
      do g = 1, self%ngrp
         fieldname = 'qcon'//trim(self%grpname(g))//'_gaex'
         call outfld(fieldname, qcon(:,:,self%ifld_grp0+g), pcols, lchnk)
         fieldname = 'qevap'//trim(self%grpname(g))//'_gaex'
         call outfld(fieldname, qevap(:,:,self%ifld_grp0+g), pcols, lchnk)
      end do
   end if

   do j = 1, self%nsoa
      write (outsoa, "(I2.2)") j
      call outfld('qcon_gaex'//outsoa, qcon(:,:,self%ifld_spc0+j), pcols, lchnk)
      call outfld('qevap_gaex'//outsoa, qevap(:,:,self%ifld_spc0+j), pcols, lchnk)
   end do

end subroutine write_fields

!----------------------------------------------------------------------
pure function soa_class(name) result(gname)

! soa class of a soa gas species: its name without the leading 'SOAG'
! and the trailing volatility bin digits

   character(len=*), intent(in) :: name
   character(len=grp_len) :: gname

   integer :: n

   n = len_trim(name)
   do while (n > 4)
      if (index('0123456789', name(n:n)) == 0) exit
      n = n - 1
   end do
   if (n > 4) then
      gname = name(5:n)
   else
      gname = ' '
   end if

end function soa_class

!----------------------------------------------------------------------
pure function class_description(gname) result(desc)

! description of a soa class, used in the field long names

   character(len=*), intent(in) :: gname
   character(len=32) :: desc

   select case (trim(gname))
   case ('ff')
      desc = 'Fossil'
   case ('bb')
      desc = 'Biomass'
   case ('bg')
      desc = 'Biogenic'
   case default
      desc = gname
   end select

end function class_description

end module soa_gaex_diag_mod
