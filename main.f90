program main

  use, intrinsic :: iso_fortran_env
  use solid_volume

  implicit none

  ! Variables =================================================================
  ! xyz -----------------------------------------------------------------------
  logical :: xyz_loaded
  character(120) :: xyz_file ! currently loaded xyz file
  integer :: xyz_a ! number of atoms
  character(2), dimension(:), allocatable :: xyz_e   ! elements
  real(REAL64), dimension(:,:), allocatable :: xyz_c ! coordinates
  integer :: xyz_center ! possible values:
                        !     0 -> center not set
                        !    -1 -> arbitrary center set
                        ! n > 0 -> n-th atom set as center
  type(point) :: arbitrary_xyz_center
  real(REAL64), dimension(:), allocatable :: xyz_dist_from_center
  integer, dimension(:), allocatable :: xyz_sorted_by_dist
  integer :: nuc_in_ball
  ! generic -------------------------------------------------------------------
  character(120) :: pname ! program name
  integer :: i
  integer :: err_n
  character(120) :: err_msg
  ! shell ---------------------------------------------------------------------
  character(*), parameter :: shell_prompt_end = "~> "
  character(120) :: shell_prompt
  logical :: shell_exit
  character(120) :: shell_buff
  character(120) :: shell_cmd

  ! Program ===================================================================
  ! init ----------------------------------------------------------------------
  shell_prompt = "(file:-)"
  xyz_loaded = .false.
  xyz_file = ""
  xyz_a = 0
  xyz_center = 0
  nuc_in_ball = 0
  arbitrary_xyz_center%x = 0.0_REAL64
  arbitrary_xyz_center%y = 0.0_REAL64
  arbitrary_xyz_center%z = 0.0_REAL64

  ! argument parsing ----------------------------------------------------------
  call get_command_argument(0,pname)

  ! interactive shell ---------------------------------------------------------
  call shell_banner()

  shell_exit = .false.
  do while (shell_exit.eqv..false.)
    write(*,'(A)',advance="no") trim(shell_prompt)//shell_prompt_end
    read(input_unit,'(A120)',iostat=err_n,iomsg=err_msg) shell_buff
!    write(*,*) "(",err_n,": ",trim(shell_buff),")"
    if (err_n == iostat_end) then
      write(*,'(A)') "exit"
      shell_cmd = "exit"
    else
      if (len_trim(shell_buff) == 0) cycle
      read(shell_buff,*) shell_cmd
    end if

    selectcase (trim(shell_cmd))
    case ("ball")
      call write_xyz_ball(shell_buff,nuc_in_ball)
    case ("mergeball")
      call mergeball(shell_buff)
    case ("center")
      call set_xyz_center(shell_buff)
    case ("centerxyz")
      call set_arbitrary_xyz_center(shell_buff,arbitrary_xyz_center)
    case ("discard")
      call discard_xyz(xyz_file,xyz_loaded)
      shell_prompt = "(file:-)"
    case ("exit","quit")
      shell_exit = .true.
    case ("help","?")
      call shell_help()
    case ("info")
      call write_info()
    case ("list")
      if (xyz_loaded.eqv..true.) then
        do i = 1, xyz_a
          write(*,'(X,I5,X,A2,3(2X,F10.4))') &
            i, xyz_e(i), xyz_c(i,1), xyz_c(i,2), xyz_c(i,3)
        end do
      else
        write(*,*) "No xyz file loaded"
      end if
    case ("listball")
      call list_ball(nuc_in_ball)
    case ("listsort")
      call write_sorted_xyz(shell_buff)
    case ("ls")
      call execute_command_line(trim(shell_buff)//" --color=auto")
    case ("molden")
      call execute_command_line(trim(shell_buff))
    case ("read")
      call read_xyz(trim(shell_buff),xyz_file,xyz_loaded)
      if (len_trim(xyz_file) == 0) then
        shell_prompt = "(file:-)"
      else
        shell_prompt = "(file:"//trim(xyz_file)//")"
      end if
    case default
      write(*,*) "Invalid shell command: ",trim(shell_cmd)
    end select
  end do

contains !=====================================================================

subroutine help(pname)

  character(*), intent(in) :: pname

  write(*,*) "Usage: ",trim(pname)

end subroutine help

!==============================================================================

subroutine shell_banner()

end subroutine shell_banner

!==============================================================================

subroutine shell_help()

  write(*,*)
  write(*,*) "Commands:"
  write(*,*) "  help, ?                      print this help"
  write(*,*) "  exit, quit                   terminate the program"
  write(*,*) "  read <xyz_file>              read and load in memory the content of xyz_file"
  write(*,*) "  discard                      discard xyz file from memory"
  write(*,*) "  center <n>                   set the n-th atom of the xyz structure as center"
  write(*,*) "  centerxyz <x> <y> <z>        set the center on the given x, y, z coordinates"
  write(*,*) "  ball <radius> [fname]        write xyz ball of a given radius"
  write(*,*) "  mergeball <b1> <b2> [out]    merge the contents of ball files b1 and b2,"
  write(*,*) "                               writing the unique nuclei in the out file"
  write(*,*) "  list                         list xyz nuclei as they appear in the xyz file"
  write(*,*) "  listball                     list xyz nuclei in the last generated ball"
  write(*,*) "                                 at increasing distance from center"
  write(*,*) "  listsort [n]                 list xyz nuclei at increasing distance from center;"
  write(*,*) "                                 if n is present, list only the first n nuclei"
  write(*,*) "  info                         write crystal analysis info"
  write(*,*) "  ls [args]                    list files"
  write(*,*) "  molden [args]                run molden"
  write(*,*)

end subroutine shell_help

!==============================================================================

subroutine polite_reminder()

  write(*,*) "Please check the help for the correct usage of this command"

end subroutine polite_reminder

!==============================================================================

subroutine set_xyz_center(str)

  character(*), intent(in) :: str
  character(120) :: s1
  
  integer :: err_n
  character(120) :: err_msg

  if (xyz_loaded.eqv..false.) then
    write(*,*) "No xyz file loaded"
  else
    read(str,*,iostat=err_n,iomsg=err_msg) s1, xyz_center
    if (err_n /= 0) then
      write(*,*) err_msg
      xyz_center = 0
    else
      if (xyz_center > xyz_a) then
        write(*,*) "Error: center out of range"
        xyz_center = 0
      else
        call sort_xyz_center()
      end if
    end if
  end if

end subroutine set_xyz_center

!==============================================================================

subroutine set_arbitrary_xyz_center(str_in,p)

  character(*), intent(in) :: str_in
  type(point), intent(out) :: p

  character(120) :: str
  real(REAL64) :: x
  real(REAL64) :: y
  real(REAL64) :: z
  integer :: err_n
  character(120) :: err_msg

  if (xyz_loaded.eqv..false.) then
    write(*,*) "No xyz file loaded"
    return
  end if

  read(str_in,*,iostat=err_n,iomsg=err_msg) str, x, y, z
  if (err_n /= 0) then
    write(*,*) "Error: "//trim(err_msg)
    call polite_reminder()
    xyz_center = 0
  else
    p%x = x
    p%y = y
    p%z = z
    xyz_center = -1
    call sort_xyz_center()
  end if

end subroutine set_arbitrary_xyz_center

!==============================================================================

subroutine sort_xyz_center()

  character(*), parameter :: my_name = "sort_xyz_center"
  logical, dimension(:), allocatable :: mask
  integer :: i
  integer :: m
  real(REAL64) :: cx
  real(REAL64) :: cy
  real(REAL64) :: cz
  integer :: err_n
  character(120) :: err_msg

  ! allocation ----------------------------------------------------------------
  if (allocated(xyz_dist_from_center)) then
    deallocate(xyz_dist_from_center)
  end if

  if (allocated(xyz_sorted_by_dist)) then
    deallocate(xyz_sorted_by_dist)
  end if

  allocate(xyz_dist_from_center(xyz_a),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    write(*,*) "Cannot allocate xyz_dist_from_center"
    stop 1
  end if

  allocate(xyz_sorted_by_dist(xyz_a),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    write(*,*) "Cannot allocate xyz_sorted_by_dist"
    stop 1
  end if

  allocate(mask(xyz_a),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    write(*,*) "Cannot allocate mask"
    stop 1
  end if

  ! set center coordinates ----------------------------------------------------
  if (xyz_center > 0) then
    cx = xyz_c(xyz_center,1)
    cy = xyz_c(xyz_center,2)
    cz = xyz_c(xyz_center,3)
  else if (xyz_center == -1) then
    cx = arbitrary_xyz_center%x
    cy = arbitrary_xyz_center%y
    cz = arbitrary_xyz_center%z
  else
    write(*,*) "Fatal error: "//my_name//": this shouldn't have happened"
    write(*,*) "Please send a bug report"
    stop 1
  end if

  ! compute distances ---------------------------------------------------------
  do i = 1, xyz_a
    xyz_dist_from_center(i) = sqrt((xyz_c(i,1)-cx)**2 + &
                                   (xyz_c(i,2)-cy)**2 + &
                                   (xyz_c(i,3)-cz)**2)
  end do

  ! sorting -------------------------------------------------------------------
  mask = .true.

  do i = 1, xyz_a
    m = minloc(xyz_dist_from_center,1,mask=mask)
    xyz_sorted_by_dist(i) = m
    mask(m) = .false.
  end do

  ! finalize ------------------------------------------------------------------
  deallocate(mask)

end subroutine sort_xyz_center

!==============================================================================

subroutine list_ball(nuc_in_ball)

  integer, intent(in) :: nuc_in_ball

  integer :: i
  integer :: s

  if (xyz_center == 0) then
    write(*,*) "Please specify a center first"
    return
  end if

  write(*,*) "   i      n Nuc      x           y           z           dist"
  do i = 1, nuc_in_ball
    s = xyz_sorted_by_dist(i)
    write(*,'(I5,": ",I5,X,A2,4(2X,F10.4))') &
      i,s,xyz_e(s),xyz_c(s,1),xyz_c(s,2),xyz_c(s,3),xyz_dist_from_center(s)
  end do

end subroutine list_ball

!==============================================================================

subroutine write_sorted_xyz(str_in)

  character(*), intent(in) :: str_in

  character(120) :: str
  integer :: i
  integer :: i_max
  integer :: s

  if (xyz_center == 0) then
    write(*,*) "Please specify a center first"
    return
  end if

  read(str_in,*,iostat=err_n,iomsg=err_msg) str, i_max
  if (err_n /= 0) then
    i_max = xyz_a
  end if

  if (i_max < 0) then
    i_max = xyz_a
  else if (i_max > xyz_a) then
    i_max = xyz_a
  end if

  write(*,*) "   i      n Nuc      x           y           z           dist"
  do i = 1, i_max
    s = xyz_sorted_by_dist(i)
    write(*,'(I5,": ",I5,X,A2,4(2X,F10.4))') &
      i,s,xyz_e(s),xyz_c(s,1),xyz_c(s,2),xyz_c(s,3),xyz_dist_from_center(s)
  end do

end subroutine write_sorted_xyz

!==============================================================================

subroutine write_info()

  character(*), parameter :: my_name = "write_info"

  ! xyz file info
  if (xyz_loaded) then
    write(*,*) "xyz file loaded: "//trim(xyz_file)
  else
    write(*,*) "No xyz file loaded"
  end if

  ! nuclei info
  write(*,*) "Nuclei: ",xyz_a

  ! center info
  if (xyz_center == 0) then
    write(*,*) "Center not set"
  else if (xyz_center == -1) then
    write(*,*) "Arbitrary center on: ", &
      arbitrary_xyz_center%x, arbitrary_xyz_center%y, arbitrary_xyz_center%z
  else if (xyz_center > 0) then
    write(*,*) "Center on atom: ",xyz_center
  else
    write(*,*) "Fatal error: "//my_name//": this shouldn't have happened"
    write(*,*) "Please send a bug report"
    stop 1
  end if

end subroutine write_info

!==============================================================================

subroutine write_xyz_ball(shell_buff,nuc_in_ball)

  character(*), intent(in) :: shell_buff
  integer, intent(out) :: nuc_in_ball
  
  integer, parameter :: fnumb = 101
  character(120) :: fname
  integer :: i
  integer :: s
  integer :: counter
  character(120) :: str
  real(REAL64) :: rad
  integer :: err_n
  character(120) :: err_msg

  ! intial check --------------------------------------------------------------
  if (xyz_center == 0) then
    write(*,*) "Please select a center first"
    return
  end if

  ! read ball radius ----------------------------------------------------------
  rad = 0.0
  fname = ""
  read(shell_buff,*,iostat=err_n,iomsg=err_msg) str, rad, fname

  if (rad <= 0.0) then
    write(*,*) "Please specify a positive ball radius"
    return
  end if

  if (len_trim(fname) == 0) then
    fname = "ball.xyz"
  end if

  ! write file ----------------------------------------------------------------
  open(unit=fnumb,file=fname,status="replace",action="write",&
    iostat=err_n,iomsg=err_msg)
  if (err_n /= 0) then
    write(*,*) err_msg
    stop 1
  end if

  counter = 0
  do i = 1, xyz_a
    s = xyz_sorted_by_dist(i)
    if (xyz_dist_from_center(s) > rad) exit
    counter = counter + 1
  end do

  nuc_in_ball = counter

  write(fnumb,*) counter
  write(fnumb,*)
  do i = 1, counter
    s = xyz_sorted_by_dist(i)
    write(fnumb,*) xyz_e(s), xyz_c(s,1), xyz_c(s,2), xyz_c(s,3)
  end do

  close(unit=fnumb,iostat=err_n,iomsg=err_msg)
  if (err_n /= 0) then
    write(*,*) err_msg
    stop 1
  end if

end subroutine write_xyz_ball

!==============================================================================

subroutine read_xyz(shell_buff,xyz_file,xyz_loaded)

  character(*), intent(in) :: shell_buff
  character(*), intent(inout) :: xyz_file
  logical, intent(inout) :: xyz_loaded
  integer, parameter :: fnumb = 100
  character(120) :: s1
  integer :: err_n
  character(120) :: err_msg

  if (xyz_loaded) then
    call discard_xyz(xyz_file,xyz_loaded)
  end if

  read(shell_buff,*,iostat=err_n,iomsg=err_msg) s1, xyz_file

  if (len_trim(xyz_file) == 0) then
    write(*,*) "Error: specify a xyz file"
    call polite_reminder()
    return
  end if

!  write(*,*) "Opening file: ",xyz_file

  open(unit=fnumb,file=trim(xyz_file),status="old",action="read",&
       iostat=err_n,iomsg=err_msg)
  if (err_n /= 0) then
    write(*,*) err_msg
    xyz_file = ""
    return
  end if

  read(fnumb,*) xyz_a
  read(fnumb,*)

  allocate(xyz_e(xyz_a),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    write(*,*) err_msg
    stop 1
  end if

  allocate(xyz_c(xyz_a,3),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    write(*,*) err_msg
    stop 1
  end if

  do i = 1, xyz_a
    read(fnumb,*) xyz_e(i), xyz_c(i,1), xyz_c(i,2), xyz_c(i,3)
  end do

!  call write_xyz(output_unit,xyz_a,xyz_e,xyz_c)

  close(unit=fnumb,iostat=err_n,iomsg=err_msg)
  if (err_n /= 0) then
    write(*,*) err_msg
    stop 1
  end if

  xyz_loaded = .true.

end subroutine read_xyz

!==============================================================================

subroutine discard_xyz(xyz_file,xyz_loaded)

  character(*), intent(out) :: xyz_file
  logical, intent(out) :: xyz_loaded

  call deallocate_xyz()

  xyz_file = ""
  xyz_loaded = .false.

end subroutine discard_xyz

!==============================================================================

subroutine deallocate_xyz()

  integer :: err_n
  character(120) :: err_msg

  xyz_a = 0
  xyz_center = 0
  arbitrary_xyz_center%x = 0.0_REAL64
  arbitrary_xyz_center%y = 0.0_REAL64
  arbitrary_xyz_center%z = 0.0_REAL64

  if (allocated(xyz_e)) then
    deallocate(xyz_e,stat=err_n,errmsg=err_msg)
    if (err_n /= 0) then
      write(*,*) "Error: "//trim(err_msg)
      stop 1
    end if
  end if

  if (allocated(xyz_c)) then
    deallocate(xyz_c,stat=err_n,errmsg=err_msg)
    if (err_n /= 0) then
      write(*,*) "Error: "//trim(err_msg)
      stop 1
    end if
  end if

end subroutine deallocate_xyz

!==============================================================================

subroutine write_xyz(fnumb,a,e,c)

  integer, intent(in) :: fnumb
  integer, intent(in) :: a
  character(2), dimension(:), intent(in) :: e
  real(REAL64), dimension(:,:), intent(in) :: c

  integer :: i

  write(fnumb,*) a
  write(fnumb,*)
  do i = 1, a
    write(fnumb,*) e(i), c(i,1), c(i,2), c(i,3)
  end do

end subroutine write_xyz

!==============================================================================

subroutine mergeball(shell_buff)

  use, intrinsic :: iso_c_binding

  interface
    integer(kind=C_INT) function mergeball_kernel(b1_c,b2_c,f_out_c) bind(C)
      use, intrinsic :: iso_c_binding
      character(kind=C_CHAR), dimension(120) :: b1_c
      character(kind=C_CHAR), dimension(120) :: b2_c
      character(kind=C_CHAR), dimension(120) :: f_out_c
    end function mergeball_kernel
  end interface

  character(*), intent(in) :: shell_buff

  character(120) :: cmd_name
  character(120) :: b1
  character(120) :: b2
  character(120) :: f_out
  integer :: i
  integer :: i_max
  integer :: err_n
  character(120) :: err_msg

  integer(kind=C_INT) :: stat
  character(kind=C_CHAR), dimension(120) :: b1_c
  character(kind=C_CHAR), dimension(120) :: b2_c
  character(kind=C_CHAR), dimension(120) :: f_out_c

  b1 = ""
  b2 = ""
  f_out = ""
  read(shell_buff,*,iostat=err_n,iomsg=err_msg) cmd_name, b1, b2, f_out

  if (len_trim(b1) == 0) then
    write(*,*) "usage: "//trim(cmd_name)//" <ball1.xyz> <ball2.xyz> [output.xyz]"
    return
  end if

  if (len_trim(b2) == 0) then
    write(*,*) "usage: "//trim(cmd_name)//" <ball1.xyz> <ball2.xyz> [output.xyz]"
    return
  end if

  if (len_trim(f_out) == 0) then
    f_out = "mergeball.xyz"
  end if

!  write(*,*) "b1    = "//trim(b1)
!  write(*,*) "b2    = "//trim(b2)
!  write(*,*) "f_out = "//trim(f_out)

  i_max = len_trim(b1)
  do i = 1, i_max
    b1_c(i) = b1(i:i)
  end do
  b1_c(i_max + 1) = C_NULL_CHAR

  i_max = len_trim(b2)
  do i = 1, i_max
    b2_c(i) = b2(i:i)
  end do
  b2_c(i_max + 1) = C_NULL_CHAR

  i_max = len_trim(f_out)
  do i = 1, i_max
    f_out_c(i) = f_out(i:i)
  end do
  f_out_c(i_max + 1) = C_NULL_CHAR

  stat = mergeball_kernel(b1_c,b2_c,f_out_c)

end subroutine mergeball

end program main
