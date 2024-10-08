!=======================================================================

! Diagnostic information output during run
!
! author: Tony Craig

      module icedrv_system

      use icedrv_kinds
      use icedrv_constants, only: nu_diag
      use icedrv_state, only: aice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      private
      public :: icedrv_system_abort, &
                icedrv_system_flush

!=======================================================================

      contains

!=======================================================================
! prints error information prior to aborting

      subroutine icedrv_system_abort(icell, istep, string, file, line)

      integer (kind=int_kind), intent(in), optional :: &
         icell       , & ! indices of grid cell where model aborts
         istep       , & ! time step number
         line            ! line number

      character (len=*), intent(in), optional :: string, file

      ! local variables

      character(len=*), parameter :: subname='(icedrv_system_abort)'

      write(nu_diag,*) ' '

      call icepack_warnings_flush(nu_diag)

      write(nu_diag,*) ' '
      write(nu_diag,*) subname,' ABORTED: '
      if (present(file))   write (nu_diag,*) subname,' called from ',trim(file)
      if (present(line))   write (nu_diag,*) subname,' line number ',line
      if (present(istep))  write (nu_diag,*) subname,' istep =', istep
      if (present(icell))  write (nu_diag,*) subname,' i, aice =', icell, aice(icell)
      if (present(string)) write (nu_diag,*) subname,' string = ',trim(string)
      call icedrv_system_flush(nu_diag)
      stop

      end subroutine icedrv_system_abort

!=======================================================================
! flushes iunit IO buffer

      subroutine icedrv_system_flush(iunit)

      integer (kind=int_kind), intent(in) :: &
         iunit        ! unit number to flush

      ! local variables

      character(len=*), parameter :: subname='(icedrv_system_flush)'

#ifndef NO_F2003
      flush(iunit)
#endif

      end subroutine icedrv_system_flush

!=======================================================================

      end module icedrv_system

!=======================================================================
