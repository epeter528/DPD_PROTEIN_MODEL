program scale

implicit none

integer :: i,natoms, counter 

integer :: frames_mod,atomnum, resnum

character(len=5) :: restype, atomtype

real               :: coord_x, coord_y ,coord_z, scaling

real               :: box_x, box_y, box_z

character(len=100) :: outputfile,filename

write(*,*) 'filename of gro-trajectory'

read(*,*) filename

write(*,*) 'pdb output ?'

read(*,*) outputfile

write(*,*) 'scaling ?'

read(*,*) scaling

write(*,*) 'frames mod ?'

read(*,*) frames_mod

counter = 0

open(unit=1,file=trim(filename))
open(unit=2,file=trim(outputfile))

do 

read(1,*,end=2,err=2)
read(1,*,end=2,err=2) natoms

counter = counter + 1

if(mod(counter,frames_mod) == 0) then

write(2,'(a6,i4)') 'MODEL ', counter

endif

do i=1,natoms

    read(1,'(i5,a5,a5,i5,3f8.3)',end=1,err=1) resnum, restype, atomtype, atomnum, &
coord_x,coord_y,coord_z

    if(atomtype(4:4) == 'C') then

       atomtype(1:1) = 'C'
       atomtype(4:4) = ' '

    endif

    if(atomtype(5:5) == 'C') then

       atomtype(1:1) = 'C'

    endif

    if(atomtype(5:5) == 'N') then

       atomtype(1:1) = 'N'

    endif

    if(atomtype(5:5) == 'O') then

       atomtype(1:1) = 'O'

    endif

    if(atomtype(5:5) == 'A') then

       atomtype(2:2) = 'A'

    endif

    if(atomtype(5:5) == 'B') then

       atomtype(2:2) = 'B'

    endif

   if(mod(counter,frames_mod) == 0)  then   

    write(2,'( A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2 )') 'ATOM  ',atomnum,trim(atomtype),' ', trim(restype),& 
'A',resnum,' ',  coord_x* scaling*10.0, coord_y * scaling*10.0, coord_z * scaling*10.0,0.0,0.0,'  ',' C'

   endif

enddo

1 continue

read(1,*,end=2,err=2) box_x,box_y,box_z

if(mod(counter,frames_mod) == 0) then

write(2,'(a3)') 'TER'
write(2,'(a6)') 'ENDMDL'

endif

enddo

2 continue

close(unit=1)
close(unit=2)


end program scale
