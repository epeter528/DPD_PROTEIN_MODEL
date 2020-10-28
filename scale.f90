program w

implicit none

character(len=1)   :: char1
character(len=100) :: filename

real(kind=8) :: scaling,x,y

integer :: p,i


write(*,*) 'file please'

read(*,*) filename

write(*,*) 'scaling'

read(*,*) scaling

open(unit=2,file='output.xvg')
open(unit=1,file=trim(filename))

p = 0

do 


   read(1,'(a1)') char1

   if(char1 .ne. '#' .and. char1 .ne. '@') then

      exit

   endif

   p = p + 1


enddo

close(unit=1)

open(unit=1,file=trim(filename))


do i=1,p 

   read(1,*)

enddo

do 


   read(1,*) x, y


   write(2,*) x,y*scaling


enddo


close(unit=1)
close(unit=2)

end program w
