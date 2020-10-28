program writegro

implicit none

integer :: natoms

integer :: i,j,k,l,m,n

integer :: time

real(kind=8) :: x,y,z

real(kind=8) :: x1,y1,z1

integer      :: resnum,atomnum

character(len=5) :: restype,atomtype

character(len=100) :: xyzfile,grofile

resnum = 1
atomnum = 1

write(*,*) 'give me the file name xyz'
read(*,*) xyzfile

write(*,*) 'give me the name gro'
read(*,*) grofile

open(unit=1,file=trim(xyzfile))
open(unit=2,file=trim(grofile))
open(unit=3,file='chain.top')

read(1,*) 
read(1,*) natoms 

write(2,*) 
write(2,*) natoms

do i = 1,natoms

   read(1,'(i2,3f7.2)') j,x,y,z

   if(j == 1) then 

      write(3,*) 'ALA 1'

   endif

   if(j == 12) then

      write(3,*) 'SOL 1'

   endif

   if(j == 16) then

      write(3,*) 'CL 1'

   endif


   if(j .ge. 1 .and. j .lt. 12) then

      atomtype = '    C'
      restype  = 'ALA  '

   endif

   if(j .ge. 12 .and. j .le. 14) then

      atomtype = '    O'
      restype  = 'SOL  '

   endif

   if(j .ge. 15 .and. j .le. 20) then

      atomtype = '    C'
      restype  = 'ALA  '

   endif

   atomnum = i

   write(2,'(i5,a5,a5,i5,f7.2,f7.2,f7.2)') resnum,restype,atomtype,atomnum,x/10.0,y/10.0,z/10.0

    if(mod(i,3)==0) then
  
     resnum = resnum + 1
  
     endif


enddo

write(2,*) ' 12.5 12.5 12.5 '

close(unit=1)
close(unit=2)

end program writegro
