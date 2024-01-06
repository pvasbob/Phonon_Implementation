
program x
integer :: i
character(len=10)  :: file_id
logical :: file_exists
do i = 1, 4
   write(file_id,'(i0)') i
   inquire(file='Output'//trim(adjustl(file_id))//'.txt',exist=file_exists)
  !inquire(file='Output2.txt',exist=file_exists)
   if(file_exists) then
     open(11,file='Output'//trim(adjustl(file_id))//'.txt')
     write(11,*) i, 'file_exist'
     close(11)
   else
     exit
   end if
   write(*,*) file_exists
end do

end program x
