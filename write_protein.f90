program write_membrane

implicit none

integer :: i,j,k,l,m,n

real    :: box_l,ran

integer  :: seed

integer  :: num_bond, num_angle

real,dimension(50000) :: x_bl,y_bl,z_bl,x,y,z

integer  :: num_dihedral,number_of_dppc,no_water

real     :: charge,scale_y,x_box,y_box,z_box

integer, dimension(50000) :: angled5,angled2,angled3,angled4,dihedral,dihedral2,dihedral3,dihedral4

integer, dimension(50000) :: angletype5,angletype2,angletype3,angletype4,dihedraltype2,dihedraltype3,dihedraltype4

integer, dimension(50000) :: dihedraltype,bonded2,bonded, angled, bondtype2,bondtype, angletype, atomtype, restype

character(len=50),dimension(50000) :: dihedralchar,bondedchar2,bondedchar, anglechar

character(len=50),dimension(50000) :: dihedralchar2,dihedralchar3,dihedralchar4

character(len=50),dimension(50000) :: anglechar5,anglechar2,anglechar3,anglechar4

integer    :: num_sol,number_tot_dppc
integer    :: p
integer    :: restype_pres

real       :: diff_x,diff_y,diff_z

real       :: range_z,d_tot,scale_x

integer    :: npos,moddar,modmin,number_of_aminos

character(len=1) , dimension(10000) :: seq

seed    = 9999

call random_seed(seed)

npos    = 0

box_l = 15 ! cubic box
num_bond  = 0
num_angle = 0

open(unit=1,file='seq.txt')

i = 1

do 

  read(1,'(a1)',end=22,err=22) seq(i)

  write(*,*) seq(i)

  i = i + 1

enddo

22 continue

close(unit=1)

number_of_aminos = i

! now build the molecule

p = 1

write(*,*) number_of_aminos

do i=1,number_of_aminos



   if(seq(i) == 'Z') then

      x(p) = abs(box_l/2.0)

      y(p) = abs(box_l/2.0)
      z(p) = 0 
      atomtype(p) = 1
      restype(p) = 1
  
      p = p + 1

   endif

   if(seq(i) == 'S') then

    do k = 1,5
 
     if(k==1) then  
 
      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p 
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1) 
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2)
      y(p) = y(p-2) + 0.2 ! pol charge2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      angled(p)   = p
      write(anglechar(p) ,'(i10,i10,i10)') p,p-2,p-1
      bondtype(p) = 2
      angletype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif
     if(k==4) then

      x(p) = x(p-3) + 0.33 ! sidechain
      y(p) = y(p-3)
      z(p) = z(p-3)
      atomtype(p) = 5
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -3, p-4
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==5) then

      x(p) = x(p-4) 
      y(p) = y(p-4)
      z(p) = z(p-4) + 0.33 ! b ackbone
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -4
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -4, p-5
      angletype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

    enddo

   endif

   if(seq(i) == 'L' .or. seq(i) == 'I' .or. seq(i) == 'A') then

    do k = 1,5

     if(k==1) then

      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2)
      y(p) = y(p-2) + 0.2 ! pol charge2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      angled(p)   = p
      write(anglechar(p) ,'(i10,i10,i10)') p,p-2,p-1
      bondtype(p) = 2
      angletype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif
     if(k==4) then

      x(p) = x(p-3) + 0.33 ! apolar sidechain
      y(p) = y(p-3)
      z(p) = z(p-3)
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -3, p-4
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==5) then

      x(p) = x(p-4)
      y(p) = y(p-4)
      z(p) = z(p-4) + 0.33 ! backbone
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -4
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -4, p-5
      angletype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

    enddo

   endif

   if(seq(i) == 'Y') then

    do k = 1,8

     if(k==1) then

      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2)
      y(p) = y(p-2) + 0.2 ! pol charge2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      angled(p)   = p
      write(anglechar(p) ,'(i10,i10,i10)') p,p-2,p-1
      bondtype(p) = 2
      angletype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif
     if(k==4) then

      x(p) = x(p-3) + 0.33 ! apolar sidechain
      y(p) = y(p-3)
      z(p) = z(p-3)
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -3, p-4
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==5) then

      x(p) = x(p-1) + 0.33 ! apolar sidechain
      y(p) = y(p-1)
      z(p) = z(p-1)
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p+1
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==6) then

      x(p) = x(p-1) + 0.33 ! apolar sidechain
      y(p) = y(p-1)
      z(p) = z(p-1)
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p+1
      angletype(p) = 3
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif
     if(k==7) then

      x(p) = x(p-2)  ! polar sidechain
      y(p) = y(p-2) + 0.33
      z(p) = z(p-2)
      atomtype(p) = 7
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      bonded2(p)   = p 
      write(bondedchar2(p),'(i10,i10)') p , p  -2
      bondtype2(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p-2
      angletype(p) = 3
      restype(p) = 1
      num_bond = num_bond + 2
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==8) then

      x(p) = x(p-7)
      y(p) = y(p-7)
      z(p) = z(p-7) + 0.33 ! backbone
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -7
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -7, p-8
      angletype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

    enddo

   endif

   if(seq(i) == 'R') then

    do k = 1,8

     if(k==1) then

      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2)
      y(p) = y(p-2) + 0.2 ! pol charge2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      angled(p)   = p
      write(anglechar(p) ,'(i10,i10,i10)') p,p-2,p-1
      bondtype(p) = 2
      angletype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif
     if(k==4) then

      x(p) = x(p-3) + 0.33 ! apolar sidechain
      y(p) = y(p-3)
      z(p) = z(p-3)
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -3, p-4
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==5) then

      x(p) = x(p-1) + 0.33 ! apolar sidechain
      y(p) = y(p-1)
      z(p) = z(p-1)
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p+1
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif
     if(k==6) then

      x(p) = x(p-1)  ! apolar sidechain
      y(p) = y(p-1) + 0.33
      z(p) = z(p-1)
      atomtype(p) = 15
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p+1
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==7) then

      x(p) = x(p-2)  ! apolar sidechain
      y(p) = y(p-2) - 0.33
      z(p) = z(p-2)
      atomtype(p) = 15
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p-2
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==8) then

      x(p) = x(p-7)
      y(p) = y(p-7)
      z(p) = z(p-7) + 0.33 ! backbone
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -7
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -7, p-8
      angletype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

    enddo

   endif


   if(seq(i) == 'T') then

    do k = 1,5
 
     if(k==1) then  
 
      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p 
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1) 
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2)
      y(p) = y(p-2) + 0.2 ! pol charge2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      angled(p)   = p
      write(anglechar(p),'(i10,i10,i10)') p , p -2,p -1
      angletype(p) = 1 
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif
     if(k==4) then

      x(p) = x(p-3) + 0.33 ! sidechain
      y(p) = y(p-3)
      z(p) = z(p-3)
      atomtype(p) = 7
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -3, p-4
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==5) then

      x(p) = x(p-4) 
      y(p) = y(p-4)
      z(p) = z(p-4) + 0.33 ! backbone
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -4
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -4, p-5
      angletype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

    enddo


   endif
   

   if(seq(i) == 'W') then

    do k=1,8

      if(k==1) then

      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2)
      y(p) = y(p-2) + 0.2 ! pol charge2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      angled(p) = p
      write(anglechar(p),'(i10,i10,i10)') p,p-2,p-1
      angletype(p) = 1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

     if(k==4) then

      x(p) = x(p-3) + 0.33 ! sidechain
      y(p) = y(p-3)
      z(p) = z(p-3)
      atomtype(p) = 8
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -3, p-4
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==5) then

      x(p) = x(p-1) + 0.33 ! sidechain
      y(p) = y(p-1)
      z(p) = z(p-1)
      atomtype(p) = 8
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p-2
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==6) then

      x(p) = x(p-1)  ! sidechain (2)
      y(p) = y(p-1) + 0.33
      z(p) = z(p-1)
      atomtype(p) = 8
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p +1
      angletype(p) = 3
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==7) then

      x(p) = x(p-2)  ! sidechain (3)
      y(p) = y(p-2)
      z(p) = z(p-2) + 0.33
      atomtype(p) = 8
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p -2 ! 3rd angle type
      angletype(p) = 3
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif


     if(k==8) then

      x(p) = x(p-7)
      y(p) = y(p-7)
      z(p) = z(p-7) + 0.33 ! backbone
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -7
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -7, p-8
      angletype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

    enddo

   endif

   if(seq(i) == 'E' .or. seq(i) == 'D') then

    do k = 1,8

     if(k==1) then

      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2)
      y(p) = y(p-2) + 0.2 ! pol charge2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      angled(p)   = p
      write(anglechar(p) ,'(i10,i10,i10)') p , p-2,p-1
      angletype(p) = 1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif
     if(k==4) then

      x(p) = x(p-3) + 0.33 ! sidechain
      y(p) = y(p-3)
      z(p) = z(p-3)
      atomtype(p) = 8 ! uncharged sidechain
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -3, p-4
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif
     if(k==5) then

      x(p) = x(p-1) + 0.33 ! sidechain
      y(p) = y(p-1)
      z(p) = z(p-1)
      atomtype(p) = 8 ! uncharged sidechain
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -3, p-4
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif


     if(k==6) then

      x(p) = x(p-1) + 0.33 ! sidechain
      y(p) = y(p-1)
      z(p) = z(p-1)
      atomtype(p) = 9 ! charged sidechain
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p-2
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif


     if(k==7) then

      x(p) = x(p-2)  ! sidechain
      y(p) = y(p-2) + 0.33
      z(p) = z(p-2)
      atomtype(p) = 9 ! charged sidechain
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -2, p-1
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==8) then

      x(p) = x(p-5)
      y(p) = y(p-5)
      z(p) = z(p-5) + 0.33 ! backbone
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -7
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -6, p-7
      angletype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

    enddo


   endif


   if(seq(i) == 'G') then

   do k=1,4

   if(k==1) then

      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2)
      y(p) = y(p-2) + 0.2 ! pol charge2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      write(anglechar(p),'(i10,i10,i10)') p , p -2 , p -1
      angled(p) = p
      angletype(p) = 1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif
     if(k == 4) then
     
     x(p) = x(p-3)
     y(p) = y(p-3) 
     z(p) = z(p-3) + 0.33
     atomtype(p) = 6
     bondtype(p) = 1
     bonded(p)   = p
     write(bondedchar(p),'(i10,i10)') p,p-3 
     angled(p)   = p
     angletype(p) = 2
     restype(p)   = 1
     write(anglechar(p),'(i10,i10,i10)') p,p-3,p-4
     p = p + 1
     num_bond = num_bond + 1
     num_angle = num_angle + 1

     endif

    enddo

   endif

   if(seq(i) == 'P') then

   do k=1,5

   if(k==1) then

      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2) 
      y(p) = y(p-2) + 0.2 ! pol charge 1
      z(p) = z(p-2)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif

     if(k==4) then

      x(p) = x(p-3) + 0.15
      y(p) = y(p-3) + 0.15 ! pol charge2
      z(p) = z(p-3)
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      write(anglechar(p),'(i10,i10,i10)') p , p -3 , p + 1
      bonded2(p) = p
      write(bondedchar2(p),'(i10,i10)') p , p + 1
      bondtype2(p) = 1
      angled(p) = p
      angletype(p) = 3
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 2
      num_angle = num_angle + 1

     endif

     if(k==5) then

      x(p) = x(p-4)
      y(p) = y(p-4) ! pol charge 1
      z(p) = z(p-4) + 0.33
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -4
      bondtype(p) = 1
      restype(p) = 1
      angled(p)  = p 
      angletype(p) = 2
      write(anglechar(p),'(i10,i10,i10)') p , p -4 , p - 5
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

    enddo

   endif

if(seq(i) == 'N' .or. seq(i) == 'Q') then

   do k=1,8

   if(k==1) then

      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2)
      y(p) = y(p-2) + 0.2 ! pol charge2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      write(anglechar(p),'(i10,i10,i10)') p , p -2 , p -1
      angled(p) = p
      angletype(p) = 1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif
     if(k==4) then

      x(p) = x(p-3) + 0.33 ! sidechain
      y(p) = y(p-3)
      z(p) = z(p-3)
      atomtype(p) = 8
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -3, p-4
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif
     if(k==5) then

      x(p) = x(p-1) + 0.33 ! sidechain
      y(p) = y(p-1)
      z(p) = z(p-1)
      atomtype(p) = 8
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p-2
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==6) then

      x(p) = x(p-1)  ! pol
      y(p) = y(p-1) - 0.2
      z(p) = z(p-1)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p+1
      angletype(p) = 1
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif
     if(k==7) then

      x(p) = x(p-2)  ! pol
      y(p) = y(p-2) + 0.2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      bondtype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1

      p = p + 1

     endif


     if(k==8) then

      x(p) = x(p-7)
      y(p) = y(p-7)
      z(p) = z(p-7) + 0.33 ! backbone
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -7
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -7, p-8
      angletype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

    enddo

   endif

   if(seq(i) == 'K') then

    do k=1,7

      if(k==1) then

      x(p) = x(p-1)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33 ! backbone central
      atomtype(p) = 2
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p - 1
      bondtype(p) = 1
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==2) then

      x(p) = x(p-1)
      y(p) = y(p-1) - 0.2 ! pol charge 1
      z(p) = z(p-1)
      atomtype(p) = 3
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1

     endif
     if(k==3) then

      x(p) = x(p-2)
      y(p) = y(p-2) + 0.2 ! pol charge2
      z(p) = z(p-2)
      atomtype(p) = 4
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -2
      write(anglechar(p),'(i10,i10,i10)') p , p -2 , p -1
      angled(p) = p
      angletype(p) = 1
      bondtype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif
     if(k==4) then

      x(p) = x(p-3) + 0.33 ! sidechain
      y(p) = y(p-3)
      z(p) = z(p-3)
      atomtype(p) = 8
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -3
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -3, p-4
      angletype(p) = 2
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==5) then

      x(p) = x(p-1)  ! sidechain (2)
      y(p) = y(p-1) + 0.33
      z(p) = z(p-1)
      atomtype(p) = 8
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p +1
      angletype(p) = 1
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif

     if(k==6) then

      x(p) = x(p-1)  ! sidechain (3)
      y(p) = y(p-1) + 0.33
      z(p) = z(p-1) 
      atomtype(p) = 10
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1 ! head charge positive
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p -2 ! 3rd angle type
      angletype(p) = 1
      restype(p) = 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

      p = p + 1

     endif


     if(k==7) then

      x(p) = x(p-6)
      y(p) = y(p-6)
      z(p) = z(p-6) + 0.33 ! backbone
      atomtype(p) = 6
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -6
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -6, p-7
      angletype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

     endif

    enddo


   endif


   if(seq(i) == 'X') then

      x(p) = x(p-1)  ! sidechain (3)
      y(p) = y(p-1)
      z(p) = z(p-1) + 0.33
      atomtype(p) = 11
      bonded(p)   = p
      write(bondedchar(p),'(i10,i10)') p , p  -1
      bondtype(p) = 1
      angled(p)  = p
      write(anglechar(p),'(i10,i10,i10)') p , p  -1, p-8
      angletype(p) = 2
      restype(p) = 1
      p = p + 1
      num_bond = num_bond + 1
      num_angle = num_angle + 1

   endif



enddo

num_dihedral = 0

number_of_dppc = p - 1 

do i=1,number_of_dppc

    if(atomtype(i) == 3) then

       do k=i+1,i+10

           if(atomtype(k) == 3) then

              dihedral(i) = i
              dihedraltype(i) = 1 ! 0 degrees

              write(dihedralchar(i),'(i4,i4,i4,i4)') i,i-1,k-1,k
             
              num_dihedral = num_dihedral + 1

              exit

           endif

       enddo

    endif

    if(atomtype(i) == 4) then


       do k=i+1,i+10

           if(atomtype(k) == 4) then

              dihedral2(i) = i
              dihedraltype2(i) = 1 ! 0 degrees

              write(dihedralchar2(i),'(i4,i4,i4,i4)') i,i-1,k-1,k
             
              num_dihedral = num_dihedral + 1

              exit

           endif

       enddo

    endif


    if(atomtype(i) == 3) then

       do k=i+2,i+10

           if(atomtype(k) == 4) then

              dihedral3(i) = i
              dihedraltype3(i) = 2 ! 180 degrees

              write(dihedralchar3(i),'(i4,i4,i4,i4)') i,i-1,k-1,k
             
              num_dihedral = num_dihedral + 1

              exit

           endif

       enddo

    endif

    if(atomtype(i) == 4) then


       do k=i+2,i+10

           if(atomtype(k) == 3) then

              dihedral4(i) = i
              dihedraltype4(i) = 2 ! 180 degrees

              write(dihedralchar4(i),'(i4,i4,i4,i4)') i,i-1,k-1,k

              num_dihedral = num_dihedral + 1

              exit

           endif

       enddo

    endif

    if(atomtype(i) == 4) then


              angled5(i) = i
              angletype5(i) = 5 ! 90 degrees

              write(anglechar5(i),'(i10,i10,i10)') i,i-2,i+1

              num_angle = num_angle + 1

    endif


    if(atomtype(i) == 3) then


              angled2(i) = i
              angletype2(i) = 5 ! 90 degrees

              write(anglechar2(i),'(i10,i10,i10)') i,i-1,i+2

              num_angle = num_angle + 1

    endif


    if(atomtype(i) == 3) then


           do k=i+1,i+10

             if(atomtype(k) == 6) then

                   angled3(i) = i
                   angletype3(i) = 6 ! 109 degrees

                  write(anglechar3(i),'(i10,i10,i10)') i,i-1,k

              num_angle = num_angle + 1

             exit

             endif

           enddo

    endif


    if(atomtype(i) == 4) then


           do k=i+1,i+10

             if(atomtype(k) == 6) then

                   angled4(i) = i
                   angletype4(i) = 7 ! 71 degrees

                  write(anglechar4(i),'(i10,i10,i10)') i,i-2,k

              num_angle = num_angle + 1

             exit

             endif

           enddo

    endif


enddo

restype_pres = 1

restype_pres = restype_pres + 1

write(*,*) number_of_dppc,'number of dppc'
write(*,*) restype_pres,'restype pres'

! now fill the rest of the box with water

x_box = 0
y_box = 0
z_box = 0
num_sol = 0
n       = 1

do

!call random_number(ran)

! check if molecule will be too close to membrane

!write(*,*) 'number_tot_dppc',number_of_dppc

     do k=1,number_of_dppc 

        diff_x = x_box - x(k)
        diff_y = y_box - y(k)
        diff_z = z_box - z(k)

        d_tot  = sqrt( diff_x**2 + diff_y**2 + diff_z**2)

        if(d_tot .lt. 0.3) then

           no_water = 1

           write(*,'(a8,2x,f7.4,2x,f7.4,2x,f7.4,2x,f7.4)') 'no_water', x_box,x(k),y(k),z(k)
 
           exit

        else

          no_water = 0

        endif

     enddo

     if(no_water == 0) then

     do i=1,3

       x(p) = x_box
       y(p) = y_box
       z(p) = z_box
       restype(p) = restype_pres
write(*,*) restype(p),'restype'


       if(i == 1) then

          atomtype(p) = 12
          num_sol = num_sol + 1

       endif

       if(i == 2) then

          atomtype(p) = 13

       endif

       if(i == 3) then

          atomtype(p) = 14

       endif

     if(i == 1) then

        bonded(p) = p
        write(bondedchar(p),'(i10,i10)') p,p+1
        bondtype(p) = 3
        num_bond = num_bond + 1

     endif

     if(i == 2) then

        bonded(p) = p
        write(bondedchar(p),'(i10,i10)') p-1,p+1
        bondtype(p) = 3
        num_bond = num_bond + 1

     endif

     if(i == 3) then

        angled(p) = p
        write(anglechar(p),'(i10,i10,i10)') p-1,p-2,p
        angletype(p) = 4
        num_angle = num_angle + 1

     endif

       write(*,*) p

       p = p + 1

     if(i == 1) then

        x_box = x_box + 0.3
        
     endif

   enddo

   endif

   x_box = x_box + 0.6

 !  write(*,*) x_box, y_box, z_box

   restype_pres = restype_pres + 1

     if( x_box .ge. box_l) then

         y_box = y_box + 0.6
         x_box = 0

     endif

     if( y_box .ge. box_l) then

         z_box = z_box + 0.6
         y_box = 0

    endif


   if(z_box .ge. box_l) then

      exit

   endif

enddo


! write everyting into the data file

open(unit=1,file='output.txt')
open(unit=2,file='chain.xyz')

write(1,*)
write(1,*)
write(1,*) p-1 , ' atoms ' 
write(1,*) num_bond, ' bonds '
write(1,*) num_angle, ' angles '
write(1,*) num_dihedral, ' dihedrals '
write(1,*)
write(1,*)
write(1,*) ' 15    atom types'
write(1,*) ' 3    bond types'
write(1,*) ' 7    angle types'
write(1,*) ' 2    dihedral types ' 
write(1,*) 
write(1,*) '0 ', box_l, ' xlo xhi '
write(1,*) '0 ', box_l, ' ylo yhi '
write(1,*) '0 ', box_l, ' zlo zhi '
write(1,*)
write(1,*) ' Masses '
write(1,*)
write(1,*) ' 1   1.00000'
write(1,*) ' 2   1.00000'
write(1,*) ' 3   1.00000'
write(1,*) ' 4   1.00000'
write(1,*) ' 5   1.00000'
write(1,*) ' 6   1.00000'
write(1,*) ' 7   1.00000'
write(1,*) ' 8   1.00000'
write(1,*) ' 9   1.00000'
write(1,*) ' 10   1.00000'
write(1,*) ' 11   1.00000'
write(1,*) ' 12   1.00000'
write(1,*) ' 13   1.00000'
write(1,*) ' 14   1.00000'
write(1,*) ' 15   1.00000'
write(1,*)
write(1,*) ' Atoms '
write(1,*) 
write(2,*) 
write(2,*) p-1

do i=1, p-1

         write(1,'(i10,i10,i10,5x,f4.1,f10.4,f10.4,f10.4,i5,i5,i5)') i,restype(i),atomtype(i),0.0,x(i),y(i),z(i),0,0,0
         write(2,'(i2,3f7.2)') atomtype(i),x(i),y(i),z(i)


enddo
write(1,*)
write(1,*) 'Bonds'
write(1,*)

k = 1

do i=1,p-1

         if(i == bonded(i)) then

                 write(1,'(i10,i10,a20)') k,bondtype(i),bondedchar(i) 

                 k = k + 1

         endif


         if(i == bonded2(i)) then

                 write(1,'(i10,i10,a20)') k,bondtype2(i),bondedchar2(i)

                 k = k + 1

         endif


enddo
write(1,*)
write(1,*) 'Angles '
write(1,*)
k = 1


do i=1,p-1

         if(i == angled(i)) then

                 write(1,'(i10,i10,a30)') k,angletype(i),anglechar(i)             

                 k = k + 1

         endif


         if(i == angled2(i)) then

                 write(1,'(i10,i10,a30)') k,angletype2(i),anglechar2(i)

                 k = k + 1

         endif

         if(i == angled3(i)) then

                 write(1,'(i10,i10,a30)') k,angletype3(i),anglechar3(i)

                 k = k + 1

         endif

         if(i == angled4(i)) then

                 write(1,'(i10,i10,a30)') k,angletype4(i),anglechar4(i)

                 k = k + 1

         endif


         if(i == angled5(i)) then

                 write(1,'(i10,i10,a30)') k,angletype5(i),anglechar5(i)

                 k = k + 1

         endif

enddo
write(1,*)
write(1,*) 'Dihedrals '
write(1,*)
k = 1


do i=1,p-1

         if(i == dihedral(i)) then

                 write(1,'(i5,i5,a18)') k,dihedraltype(i),dihedralchar(i)

                 k = k + 1

         endif


         if(i == dihedral2(i)) then

                 write(1,'(i5,i5,a18)') k,dihedraltype2(i),dihedralchar2(i)

                 k = k + 1

         endif

         if(i == dihedral3(i)) then

                 write(1,'(i5,i5,a18)') k,dihedraltype3(i),dihedralchar3(i)

                 k = k + 1

         endif

         if(i == dihedral4(i)) then

                 write(1,'(i5,i5,a18)') k,dihedraltype4(i),dihedralchar4(i)

                 k = k + 1

         endif

enddo


close(unit=1)
close(unit=2)


end program write_membrane
