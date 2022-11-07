!==========================================================================================!
!==========================================================================================!
!   module hurricane_coms                                                                  !
!   This module contains variables used to control hurricanes.                             !
!------------------------------------------------------------------------------------------!
module hurricane_coms
   use ed_max_dims, only : str_len      & ! intent(in)
                         , max_hurricanes  ! ! intent(in)
   implicit none


   !> Whether to include hurricanes. 0 = no; 1 = yes. If hurricanes are included, the
   !> storm schedule must be included in a text file.
   integer :: include_hurricanes

   !----- Hurricane specification list ----------------------------------------------------!
   !> A structure for holding the info for a single hurricane
   type hurrtime
       !> Year hurricane occurs
       integer :: year
       !> Month hurricane occurs, 1-12
       integer :: month
       !> Between 0 and 1; 1 is a max severity storm
       real    :: severity
   end type hurrtime

   !> Filename containing the hurricane list
   character(len=str_len) :: hurricane_db

   !> List of scheduled hurricanes
   type(hurrtime), dimension(max_hurricanes) :: hurricane_db_list

   !> Total length of the hurricane schedule
   integer :: hurricane_db_list_len

   !> Counter for cohorts off allometry
   integer :: coh_off_allom

   !> Filename of detailed hurricane report, if the user asks for it
   character(len=str_len) :: hurricane_report
   !=======================================================================================!
   !=======================================================================================!

end module hurricane_coms
!==========================================================================================!
!==========================================================================================!
