module json_loader
  use const_props_mod, only: const_props_type
  use json_module, only: json_file, json_value, json_core
  use ccpp_kinds, only: kind_phys
  use cam_abortutils,  only: endrun
  use spmd_utils,      only: masterproc
  use cam_logfile,     only: iulog

  implicit none
  
contains

  subroutine json_loader_read( jsonfile, cnst_info, ncnst, nrxtn, nphot, jnames )
    
    character(len=*), intent(in) :: jsonfile
    type(const_props_type), allocatable :: cnst_info(:)
    integer, intent(out) :: ncnst, nrxtn, nphot
    character(len=16), allocatable, intent(out) :: jnames(:)
    
    ! local vars
    type(json_file) :: json       !! the JSON structure read from the file
    type(json_value),pointer :: p !! a pointer for low-level manipulations
    type(json_core) :: core       !! factory for manipulating `json_value` pointers
    type(json_value),pointer :: child 
    type(json_value),pointer :: child2
    type(json_value),pointer :: child3
    character(len=:),allocatable :: name
    
    logical :: found
    integer :: i, n, nsections
    character(len=:),allocatable :: string
    real(kind_phys) :: rval
    
    call json%initialize()

    if (masterproc) write(iulog,*) 'json_loader_read Load the file :'//jsonfile

    call json%load_file(filename = jsonfile)

!    call json%print_file()

    call core%initialize()
    call json%get(p) ! get root

    call core%get_child(p,child)

    nsections = core%count(child)

    !write(*,*)  'nsections : ', nsections

    do i = 1,nsections

       call core%get_child(child,i,child2)
       call core%info(child2,name=name)

       molecules: if (name=='molecules') then

          !write(*,*)  'Read obj data : '//name
          ncnst = core%count(child2)
          !write(*,*)  '  ncnst : ', ncnst

          allocate( cnst_info(ncnst) )

          do n = 1,ncnst
             call core%get_child(child2, n, child3, found)
             if (found) then
                call core%get(child3,'moleculename',string)
                call cnst_info(n)%set_name(string)
                deallocate(string)

                call core%get(child3,'formula',string)
                call cnst_info(n)%set_desc(string)
                deallocate(string)

                call core%get(child3,'molecular_weight',string)
                read( string, * ) rval
                deallocate(string)
                call cnst_info(n)%set_wght(rval)
             else
                !write(*,*) 'json_loader_read ERROR: Did not find molecules child ',n
                call endrun('json_loader_read ERROR: Did not find molecules child ')
             endif
          enddo

       end if molecules
       
       photolysis: if (name=='photolysis') then
          nphot = core%count(child2)
          !write(*,*)  '  nphot : ', nphot

          allocate( jnames(nphot) )
          do n = 1,nphot
             call core%get_child(child2, n, child3, found)
             if (found) then
                call core%get(child3,'rate',string)
                 jnames(n) = string
             else
                !write(*,*) 'json_loader_read ERROR: Did not find photolysis child ',n
                call endrun('json_loader_read ERROR: Did not find photolysis child ')
             endif
          end do

       end if photolysis
       
       reactions: if (name=='reactions') then
          nrxtn = core%count(child2)
          !write(*,*)  '  nrxtn : ', nrxtn
       end if reactions
       
    enddo

    call json%destroy()

  end subroutine json_loader_read

end module json_loader
