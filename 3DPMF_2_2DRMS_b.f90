program pmf22d


implicit none

integer :: i,j,k,l,m,n

real(kind=8),dimension(10000,10000) :: prob, rmsd, rgyr

real(kind=8),dimension(100000)       :: probmax

real(kind=8) :: maxi,probber

maxi = 0

open(unit=1,file='normFEL.dat')
open(unit=2,file='normFEL2.dat')

k = 1

do

    do i=1,100
    
        read(1,'(f5.3,1x,f6.4,1x,f4.2)',end=1,err=1) rgyr(i,k),rmsd(i,k),prob(i,k)
   
        probber = probber + prob(i,k)
 
    enddo

    read(1,*,end=1,err=1)
    
    k = k + 1
    
enddo

1 continue

close(unit=1)

do l=1,k-1  
     
     do i=1,100
         
         write(2,'(f5.3,1x,f6.4,1x,f8.6)') rgyr(i,l),rmsd(i,l),prob(i,l)/probber
        
         if(maxi .lt. prob(i,l)/probber) then

            maxi = prob(i,l)/probber

        endif
 
     enddo
     
     write(2,*)
     
 enddo

write(*,*) maxi

close(unit=2)


end program pmf22d
