program bonds

implicit none

character(len=7) :: char7

integer  :: k,i,n, restype,atomtype,dummy2,dummy3,dummy4

real     :: dummy

real,dimension(1000000) :: x,y,z

integer, dimension(100000) :: bondtype,bonded_i,bonded_k

integer, dimension(100000) :: angletype,angle_i,angle_k,angle_j

integer, dimension(100000) :: dihedraltype,dih_i,dih_k,dih_j,dih_l

real   :: d_x,d_y,d_z,d_tot

open(unit=1,file='output.txt')

do k=1,50

   read(1,'(a7)') char7
   
!   write(*,*) char7
   
   if(char7(3:7) == 'Atoms') then
   
!      write(*,*) 'here'
   
      read(1,*) 
!      write(*,*) char7
      
      i = 1
      
      do 
      
         read(1,'(i10,i10,i10,5x,f4.1,f10.4,f10.4,f10.4,i5,i5,i5)',end=1,err=1) n,restype,atomtype,&
         dummy,x(i),y(i),z(i),dummy2,dummy3,dummy4
  !       write(*,'(i10,i10,i10,5x,f4.1,f10.4,f10.4,f10.4,i5,i5,i5)') i,restype,atomtype,dummy,x(i),&
  !       y(i),z(i),dummy2,dummy3,dummy4
   
  !        read(1,*) char7
  !        write(*,*) char7
  
         i = i + 1
   
      enddo
      
   1 continue
 
         i = 1
 
      do
      
        read(1,'(i10,i10,i10,i10)',end=2,err=2) n,bondtype(i),bonded_i(i),bonded_k(i) 
   !     write(*,'(i10,i10,i10,i10)') i,bondtype(i),bonded_i(i),bonded_k(i)
   
        d_x = x(bonded_i(i)) - x(bonded_k(i))
        d_y = y(bonded_i(i)) - y(bonded_k(i))
        d_z = z(bonded_i(i)) - z(bonded_k(i))        
    
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)
        
        if(d_tot .gt. 2.0) then
        
           write(*,*) bonded_i(i),bonded_k(i),d_tot, '! too large bond '
           
        endif
        if(d_tot == 0.0) then
        
           write(*,*) bonded_i(i),bonded_k(i),d_tot, '! 0 bond '
           
        endif  
        
        i = i + 1
    
      enddo
    
    2 continue
    
        i = 1
    
        read(1,*)
    
      do
      
        read(1,'(i10,i10,i10,i10,i10)',end=3,err=3) n,angletype(i),angle_i(i),angle_k(i),angle_j(i) 
  !      write(*,'(i10,i10,i10,i10,i10)')  i,angletype(i),angle_i(i),angle_k(i),angle_j(i) 
    
        d_x = x(angle_i(i)) - x(angle_j(i))
        d_y = y(angle_i(i)) - y(angle_j(i))
        d_z = z(angle_i(i)) - z(angle_j(i))
        
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)
        
        if(d_tot .gt. 2.0) then
        
           write(*,*) angle_i(i),angle_j(i),d_tot, '! too large angle '
           
        endif    
        if(d_tot == 0) then
        
           write(*,*) angle_i(i),angle_j(i),d_tot, '! 0 angle '
           
        endif         
        
        d_x = x(angle_i(i)) - x(angle_k(i))
        d_y = y(angle_i(i)) - y(angle_k(i))
        d_z = z(angle_i(i)) - z(angle_k(i))
        
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)
        
        if(d_tot .gt. 2.0) then
        
           write(*,*) angle_i(i),angle_k(i),d_tot, '! too large angle '
           
        endif         
        
        if(d_tot == 0.0) then
        
           write(*,*) angle_i(i),angle_k(i),d_tot, '! 0 angle '
           
        endif            
        
    
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
        
        if(d_tot .gt. 2.0) then
        
           write(*,*) dih_i(i),dih_j(i),dih_k(i),dih_l(i),d_tot, '! too large dihedral '
           
        endif
        if(d_tot == 0) then
        
           write(*,*) dih_i(i),dih_j(i),dih_k(i),dih_l(i),d_tot, '! 0 dihedral '
           
        endif        

        d_x = x(dih_i(i)) - x(dih_k(i))
        d_y = y(dih_i(i)) - y(dih_k(i))
        d_z = z(dih_i(i)) - z(dih_k(i))
        
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)
        
        if(d_tot .gt. 2.0) then
        
           write(*,*) dih_i(i),dih_j(i),dih_k(i),dih_l(i),d_tot, '! too large dihedral '
           
        endif
         if(d_tot == 0) then
        
           write(*,*) dih_i(i),dih_j(i),dih_k(i),dih_l(i),d_tot, '! 0 dihedral '
           
        endif       
        d_x = x(dih_i(i)) - x(dih_l(i))
        d_y = y(dih_i(i)) - y(dih_l(i))
        d_z = z(dih_i(i)) - z(dih_l(i))
        
        d_tot = sqrt(d_x**2+d_y**2+d_z**2)
        
        if(d_tot .gt. 2.0) then
        
           write(*,*) dih_i(i),dih_j(i),dih_k(i),dih_l(i),d_tot, '! too large dihedral '
           
        endif
        if(d_tot == 0) then
        
           write(*,*) dih_i(i),dih_j(i),dih_k(i),dih_l(i),d_tot, '! 0 dihedral '
           
        endif        
        
        
        i = i + 1
      
      enddo
    
          
    endif
 
enddo

4 continue

close(unit=1)

end program bonds