program main

  use iso_fortran_env

  ! Variables =================================================================
  ! xyz -----------------------------------------------------------------------
  logical :: xyz_loaded
  character(120) :: xyz_file ! currently loaded xyz file
  integer :: xyz_a ! number of atoms
  character(2), dimension(:), allocatable :: xyz_e   ! elements
  real(REAL64), dimension(:,:), allocatable :: xyz_c ! coordinates
  integer :: xyz_center ! center
  real(REAL64), dimension(:), allocatable :: xyz_dist_from_center
  integer, dimension(:), allocatable :: xyz_sorted_by_dist
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
      call write_xyz_ball(shell_buff)
    case ("center")
      call set_xyz_center(shell_buff)
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
    case ("listsort")
      call write_sorted_xyz()
    case ("ls")
      call execute_command_line(trim(shell_buff)//" --color=auto")
    case ("molden")
      call execute_command_line(trim(shell_buff))
    case ("read")
      call read_xyz(trim(shell_buff),xyz_file,xyz_loaded)
      shell_prompt = "(file:"//trim(xyz_file)//")"
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
  write(*,*) "  ball radius [fname]     write xyz ball of a given radius"
  write(*,*) "  center                  select center of xyz structure"
  write(*,*) "  discard                 discard xyz file from memory"
  write(*,*) "  exit, quit              terminate the program"
  write(*,*) "  help, ?                 print this help"
  write(*,*) "  info                    write crystal analysis info"
  write(*,*) "  list                    list xyz nuclei"
  write(*,*) "  listsort                list xyz nuclei at increasing distance from center"
  write(*,*) "  ls [args]               list files"
  write(*,*) "  molden [args]           run molden"
  write(*,*) "  read xyz_file           read and load in memory the content of xyz_file"
  write(*,*)

end subroutine shell_help

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

subroutine sort_xyz_center()

  logical, dimension(:), allocatable :: mask
  integer :: i
  integer :: m
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

  ! compute distances ---------------------------------------------------------
  do i = 1, xyz_a
    xyz_dist_from_center(i) = sqrt((xyz_c(i,1)-xyz_c(xyz_center,1))**2 + &
                                   (xyz_c(i,2)-xyz_c(xyz_center,2))**2 + &
                                   (xyz_c(i,3)-xyz_c(xyz_center,3))**2)
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

subroutine write_sorted_xyz()

  integer :: i
  integer :: s

  if (xyz_center == 0) then
    write(*,*) "Please specify a center first"
    return
  end if

  write(*,*) "  i     Nuc        x           y           z           dist"
  do i = 1, xyz_a
    s = xyz_sorted_by_dist(i)
    write(*,'(I4,": ",I4,X,A2,4(2X,F10.4))') &
      i,s,xyz_e(s),xyz_c(s,1),xyz_c(s,2),xyz_c(s,3),xyz_dist_from_center(s)
  end do

end subroutine write_sorted_xyz

!==============================================================================

subroutine write_info()

  write(*,*) "Nuclei: ",xyz_a
  write(*,*) "Center: ",xyz_center

end subroutine write_info

!==============================================================================

subroutine write_xyz_ball(shell_buff)

  character(*), intent(in) :: shell_buff
  
  integer, parameter :: fnumb = 101
  character(120) :: fname
  integer :: i
  integer :: s
  integer :: counter
  character(120) :: str
  real(REAL64) :: rad
  integer :: err_n
  character(120) :: err_msg

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

  read(shell_buff,*,iostat=err_n,iomsg=err_msg) s1, xyz_file

  if (len_trim(xyz_file) == 0) then
    write(*,*) "Please specify a xyz file"
  else
    if (xyz_loaded.eqv..true.) call deallocate_xyz()

!    write(*,*) "Opening file: ",xyz_file

    open(unit=fnumb,file=trim(xyz_file),status="old",action="read",&
         iostat=err_n,iomsg=err_msg)
    if (err_n /= 0) then
      write(*,*) err_msg
    else
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

  !    call write_xyz(output_unit,xyz_a,xyz_e,xyz_c)

      close(unit=fnumb,iostat=err_n,iomsg=err_msg)
      if (err_n /= 0) then
        write(*,*) err_msg
        stop 1
      end if
      
      xyz_loaded = .true.
    end if
  end if

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

  if (allocated(xyz_e)) then
    deallocate(xyz_e,stat=err_n,errmsg=err_msg)
    if (err_n /= 0) then
      write(*,*) err_msg
      stop 1
    end if
  end if

  if (allocated(xyz_c)) then
    deallocate(xyz_c,stat=err_n,errmsg=err_msg)
    if (err_n /= 0) then
      write(*,*) err_msg
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

end program main
