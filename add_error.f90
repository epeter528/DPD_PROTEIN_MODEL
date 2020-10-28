program add_error

implicit none


real(kind=8) :: x,y

character(len=100) :: output,filename

integer :: i

write(*,*) 'give me the filename'

read(*,*) filename
read(*,*) output

open(unit=1,file=trim(filename))
open(unit=2,file=trim(output))

do

   read(1,*,end=2,err=2) x,y
   if(mod(i,50) == 0) then

   write(2,*) x,y,0.4
 
   endif

   i = i + 1

enddo

2 continue

close(unit=1)
close(unit=2)


end program add_error
