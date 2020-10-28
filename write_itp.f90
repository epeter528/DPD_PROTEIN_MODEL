program bonds

implicit none

character(len=7) :: char7

integer  :: k,i,n,dummy2,dummy3,dummy4

integer  :: natoms, hem_int, ident_prot_max,ident_prot_min,num_bonds

real     :: dummy

integer :: count_prot, count_prot2 

real,dimension(1000000) :: x,y,z

integer, dimension(100000) :: atomtype,restype,bondtype,bonded_i,bonded_k

integer, dimension(100000) :: angletype,angle_i,angle_k,angle_j

integer, dimension(100000) :: dihedraltype,dih_i,dih_k,dih_j,dih_l

real   :: d_x,d_y,d_z,d_tot

 count_prot  = 0
 count_prot2 = 0
 
open(unit=1,file='output.txt')
open(unit=3,file='owow.ndx')

do k=1,50

   read(1,'(a7)') char7
   
!   write(*,*) char7
   
   if(char7(3:7) == 'Atoms') then
   
!      write(*,*) 'here'
   
      read(1,*) 
!      write(*,*) char7
      
      i = 1
      
      do 
      
         read(1,'(i10,i10,i10,5x,f4.1,f10.4,f10.4,f10.4,i5,i5,i5)',end=1,err=1) n,restype(i),atomtype(i),&
         dummy,x(i),y(i),z(i),dummy2,dummy3,dummy4
  !       write(*,'(i10,i10,i10,5x,f4.1,f10.4,f10.4,f10.4,i5,i5,i5)') i,restype,atomtype,dummy,x(i),&
  !       y(i),z(i),dummy2,dummy3,dummy4
   
  !        read(1,*) char7
  !        write(*,*) char7
  
         i = i + 1
   
      enddo
      
   1 continue
 
       natoms = i - 1
 
         i = 1
 
      do
      
        read(1,'(i10,i10,i10,i10)',end=2,err=2) n,bondtype(i),bonded_i(i),bonded_k(i) 
   !     write(*,'(i10,i10,i10,i10)') i,bondtype(i),bonded_i(i),bonded_k(i)
   
        d_x = x(bonded_i(i)) - x(bonded_k(i))
        d_y = y(bonded_i(i)) - y(bonded_k(i))
        d_z = z(bonded_i(i)) - z(bonded_k(i))        
    
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)
        
        i = i + 1
    
      enddo
    
    2 continue
    
       num_bonds = i-1
    
        i = 1
    
        read(1,*)
    
      do
      
        read(1,'(i10,i10,i10,i10,i10)',end=3,err=3) n,angletype(i),angle_i(i),angle_k(i),angle_j(i) 
  !      write(*,'(i10,i10,i10,i10,i10)')  i,angletype(i),angle_i(i),angle_k(i),angle_j(i) 
    
        d_x = x(angle_i(i)) - x(angle_j(i))
        d_y = y(angle_i(i)) - y(angle_j(i))
        d_z = z(angle_i(i)) - z(angle_j(i))
        
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)       
        
        d_x = x(angle_i(i)) - x(angle_k(i))
        d_y = y(angle_i(i)) - y(angle_k(i))
        d_z = z(angle_i(i)) - z(angle_k(i))
        
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)

        i = i + 1
    
      enddo
    
    3 continue
    
       read(1,*)
       
       i = 1
       
      do 
      
        read(1,'(i10,i10,i10,i10,i10,i10)',end=4,err=4) n, dihedraltype(i), dih_i(i),dih_j(i),dih_k(i),dih_l(i)      
    !    write(*,'(i10,i10,i10,i10,i10,i10)') n, dihedraltype(i), dih_i(i),dih_j(i),dih_k(i),dih_l(i) 

        d_x = x(dih_i(i)) - x(dih_j(i))
        d_y = y(dih_i(i)) - y(dih_j(i))
        d_z = z(dih_i(i)) - z(dih_j(i))
        
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)       

        d_x = x(dih_i(i)) - x(dih_k(i))
        d_y = y(dih_i(i)) - y(dih_k(i))
        d_z = z(dih_i(i)) - z(dih_k(i))
        
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)
              
        d_x = x(dih_i(i)) - x(dih_l(i))
        d_y = y(dih_i(i)) - y(dih_l(i))
        d_z = z(dih_i(i)) - z(dih_l(i))
        
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)
        
             
        i = i + 1
      
      enddo
    
          
    endif
 
enddo

4 continue

close(unit=1)

open(unit=1,file='topol.top')

write(1,'(a39)') '#include "amber99sb.ff/forcefield.itp"'


do i=1,natoms

    if(atomtype(i) .ge. 17 .and. atomtype(i) .le. 20) then
    
       hem_int = 1
       exit
   
    else
    
       hem_int = 0
       
    endif
    
enddo    

if(hem_int == 1) then

write(1,'(a16)') '[ moleculetype ]'
write(1,'(a10)') 'HEM  3'
write(1,'(a9)') '[ atoms ]'

i = 1

do n=1,natoms

   if(atomtype(n) .ge. 17 .and. atomtype(i) .le. 20) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '   C ',i,' HEM ','   C ',i,'  0.00 ','  1.00 '   
  i = i + 1   
   endif     

 
   
enddo

write(1,'(a9)') '[ bonds ]'

endif

i = 1

write(1,'(a16)') '[ moleculetype ]'
write(1,'(a10)') 'Protein  3'
write(1,'(a9)') '[ atoms ]'
!      1         CT      1    ACE    CH3      1    -0.3662      12.01
do n=1,natoms

   if(atomtype(n) == 1 .and. count_prot == 0) then
   
      ident_prot_min = n
      
      count_prot = count_prot + 1
    
   endif

   if(atomtype(n) == 1) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '   C ',i,' ALA ','   N ',i,'  0.00 ','  1.00 '
   i = i + 1
   endif  
   if(atomtype(n) == 2) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '   C ',i,' ALA ','  CA ',i,'  0.00 ','  1.00 '
   i = i + 1
   endif     
   if(atomtype(n) == 3) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '   O ',i,' ALA ','  N  ',i,' -0.25 ','  1.00 '
   i = i + 1
   endif     
   if(atomtype(n) == 4) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '   O ',i,' ALA ','  C  ',i,'  0.25 ','  1.00 '
   i = i + 1
   endif     
   
   if(atomtype(n) .ge. 5 .and. atomtype(n) .le. 11) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '   C ',i,' ALA ','  CB ',i,'  0.00 ','  1.00 '
   i = i + 1
   endif   
   if(atomtype(n) == 15) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '   C ',i,' ALA ','   O ',i,'  0.00 ','  1.00 '
   i = i + 1
   endif     
   
   if(atomtype(n) == 12 .and. count_prot2 == 0) then
   
      ident_prot_max = n-1
      
      count_prot2 = count_prot2 + 1
      
   endif 
  
     
enddo

do i=1,num_bonds

   do n=1,natoms
   
       if(n == bonded_i(i)) then
       
          bonded_i(i) = bonded_i(i) - ident_prot_min + 1
          
       endif
       if(n == bonded_k(i)) then
       
          bonded_k(i) = bonded_k(i) - ident_prot_min + 1
          
       endif       
       
    enddo
    
enddo

write(*,*) ident_prot_min, ident_prot_max

write(1,'(a9)') '[ bonds ]'

do i=1,num_bonds

    if(bonded_i(i) .ge. 1 .and. bonded_i(i) .le. ident_prot_max-ident_prot_min &
    .and. bonded_k(i) .ge. 1 .and. bonded_k(i) .le. ident_prot_max-ident_prot_min) then
    
        write(1,'(i10,i10,i10)') bonded_i(i), bonded_k(i), 1
     
    endif

enddo


i = 1

write(1,'(a16)') '[ moleculetype ]'
write(1,'(a20)') 'non-protein  3'
write(1,'(a9)') '[ atoms ]'

do n=1,natoms

   if(atomtype(n) == 12) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '  OW ',i,' SOL ','  OW ',i,'  0.00 ','  1.00 '   
   i = i + 1  
     write(3,*) n
 
   endif   
   if(atomtype(n) == 13) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '   H ',i,' SOL ','  H ',i,' -0.75 ','  1.00 '   
   i = i + 1   
   endif
   if(atomtype(n) == 14) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '   H ',i,' SOL ','  H ',i,'  0.75 ','  1.00 '   
   i = i + 1   
   endif
   if(atomtype(n) == 16) then
   
     write(1,'(i10,a5,i10,a5,a5,i10,a7,a7)') i, '  NA ',i,'  NA ','  NA ',i,'  1.00 ','  1.00 '   
   i = i + 1   
   endif   
   

   
enddo

write(1,'(a10)') '[ system ]' 
write(1,'(a13)')'[ molecules ]'

write(1,*) 'HEM 1 '
write(1,*) 'Protein 1'
write(1,*) 'non-protein 1'

close(unit=1)
close(unit=3)

end program bonds
