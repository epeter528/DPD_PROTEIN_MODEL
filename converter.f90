program converter

implicit none

integer :: i,j,k,l

integer :: natoms1,natoms2

character(len=100) :: inputfile,inputfile2,outputfile

integer :: resnum,resnum2,atomnum,atomnum2

real(kind=8) :: x,x2,y,y2,z,z2

real(kind=8) :: box_x,box_y,box_z

character(len=5) :: atomtype,atomtype2,restype,restype2

resnum2 = 1
atomnum2 = 1

write(*,*) 'input 1'
read(*,*) inputfile

write(*,*) 'output'
read(*,*) outputfile


open(unit=1,file=trim(inputfile))
open(unit=3,file=trim(outputfile))

read(1,*)
read(1,*) natoms1

write(3,*)
write(3,*) natoms1

do i=1,natoms1

   read(1,'(i5,a5,a5,i5,3f8.3)') resnum,restype,atomtype,atomnum,x,y,z

   x = x*0.561
   y = y*0.561
   z = z*0.561

   if(atomtype == '   CA') then

     write(3,'(i5,a5,a5,i5,3f8.3)') resnum2,'GLY  ','    N',atomnum2,x,y,z
     write(3,'(i5,a5,a5,i5,3f8.3)') resnum2,'GLY  ','    H',atomnum2+1,x-0.075,y-0.051,z-0.041
     write(3,'(i5,a5,a5,i5,3f8.3)') resnum2,'GLY  ','   CA',atomnum2+2,x-0.03,y+0.117,z+0.082
     write(3,'(i5,a5,a5,i5,3f8.3)') resnum2,'GLY  ','    C',atomnum2+3,x+0.013,y+0.096,z+0.226
     write(3,'(i5,a5,a5,i5,3f8.3)') resnum2,'GLY  ','    O',atomnum2+4,x+0.071,y+0.186,z+0.289

    atomnum2 = atomnum2 + 5 
    resnum2 = resnum2 + 1

   endif

!   10GLY     N1  114   0.766   0.937   0.682
!   10GLY      H  115   0.691   0.886   0.641
!   10GLY     CA  116   0.736   1.054   0.764
!   10GLY      C  117   0.779   1.033   0.908
!   10GLY      O  118   0.837   1.123   0.971

enddo

read(1,*) box_x,box_y,box_z
write(3,*) box_x,box_y,box_z

close(unit=1)
close(unit=3)

end program converter

