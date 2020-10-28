program pmf22d


implicit none

integer :: i,j,k,l,m,n

real(kind=8),dimension(10000,10000) :: prob, rmsd, rgyr

real(kind=8),dimension(100000)       :: probmax

open(unit=1,file='normFEL.dat')
open(unit=2,file='PMF2DRMSD.dat')

k = 1

do

    do i=1,100
    
        read(1,'(f5.3,1x,f6.4,1x,f4.2)',end=1,err=1) rgyr(i,k),rmsd(i,k),prob(i,k)
    
    enddo

    read(1,*,end=1,err=1)
    
    k = k + 1
    
enddo

1 continue

close(unit=1)

do l=1,100

   probmax(l) = 0

   do i=1,k
    
        if(prob(l,i) .gt. probmax(l)) then
    
           probmax(l) = prob(l,i)
    
        endif

   enddo

enddo   

do l=1,100

if(probmax(l) .ne. 0) then

   write(2,*) rmsd(l,1), -log(probmax(l)) 

endif   
   
enddo

close(unit=2)


end program pmf22d
