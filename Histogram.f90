program histogram

implicit none

integer :: i,j,k,l,m,n,o,p


real(kind=8) , dimension(100000) :: x,y

character(len=100) :: filename,output

character(len=1)   :: char1

real(kind=8)   , dimension(1000000) :: hisval_x, hisval_y 

real(kind=8)       :: range_x, range_y

integer   :: number_of_rangers

write(*,*) 'Filename ?'

read(*,*) filename

write(*,*) 'Output ?'

read(*,*) output

write(*,*) 'lowest range'

read(*,*) range_x

write(*,*) 'largest range'

read(*,*) range_y

open(unit=1,file=trim(filename))

do
 
    read(1,'(a1)') char1
    
    if(char1 .ne. '@' .and. char1 .ne. '#') then
    
       exit
       
    endif
    
enddo

i = 1

do 

    read(1,*,end=1,err=1) x(i),y(i)
    
    i = i + 1
    
enddo

1 continue

close(unit=1)

number_of_rangers = i


do i = 1,10000 

      hisval_x(i) = range_x
      
      do k = 1,number_of_rangers
      
          if( hisval_x(i-1) .le. y(k) .and. hisval_x(i) .ge. y(k)) then
          
              hisval_y(i) = hisval_y(i) + 1     
          
          endif
      
      enddo

      range_x = range_x + (range_y - range_x)/10000.0

enddo

open(unit=1,file=trim(output))

do i=1,10000

   if(hisval_y(i) .ne. 0) then

   write(1,*) hisval_x(i), hisval_y(i)

   endif
   
enddo

close(unit=1)

end program histogram


