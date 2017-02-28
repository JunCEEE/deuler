

program random 
implicit none 
	integer :: i 
	DOUBLE PRECISION, PARAMETER :: pi = 3.1415926536
	DOUBLE PRECISION :: x,y,z,lx,ly,lz,latx,laty,latz
        DOUBLE PRECISION :: alfa,beta,gama,ang
	INTEGER :: numgrain,seed 
	DOUBLE PRECISION, DIMENSION(3) :: l,m,n
	DOUBLE PRECISION :: phi1,phi2,phi3,mod1,mod2,mod3,dx1,dx2,dx3
	OPEN(7,FILE='info.in')
	OPEN(6,FILE='input.txt')
	READ(7,*)lx,ly,lz,numgrain
	READ(7,*)latx,laty,latz
!	WRITE(6,*)'AA'
	WRITE(6,*)latx,' # lattice parameter'
	WRITE(6,72)numgrain,lx,ly,lz !,'# number of grains, box-x, box-y, box-z'
   72  format(i5,3(1x,f10.5)," # number of grains, box-x, box-y, box-z")
	WRITE(6,*)'0.0 1.0'	
	
	call init_random_seed() ! 系统根据日期和时间随机地提供种子 
	do i=1,numgrain 
		call random_number (x)
		call random_number (y)	 ! 每次的随机数就都不一样了 
                call random_number (z)
                call random_number (alfa)
                call random_number (beta)
                call random_number (gama)
		ang=acos(beta*2.0-1)
		phi1=alfa*360
		phi2=ang/pi*180
		phi3=gama*360
		l(1)=cos(phi2)*sin(phi1)
		m(1)=sin(phi2)*sin(phi1)
		n(1)=cos(phi1)
	l(2)=(l(1)*l(1)+(m(1)*m(1)+n(1)*n(1))*cos(phi3))*(-1.0)*sin(phi2)+(l(1)*m(1)*(1-cos(phi3))-n(1)*sin(phi3))*cos(phi2)
	m(2)=(l(1)*m(1)*(1-cos(phi3))+n(1)*sin(phi3))*(-1.0)*sin(phi2)+(m(1)*m(1)+(l(1)*l(1)+n(1)*n(1))*cos(phi3))*cos(phi2)
	n(2)=(l(1)*n(1)*(1-cos(phi3))-m(1)*sin(phi3))*(-1.0)*sin(phi2)+(m(1)*n(1)*(1-cos(phi3))+l(1)*sin(phi3))*cos(phi2)
	l(3)=m(1)*n(2)-n(1)*m(2)
	m(3)=-1.0*l(1)*n(2)+n(1)*l(2)
	n(3)=l(1)*m(2)-m(1)*l(2)
	mod1=l(1)*l(1)+m(1)*m(1)+n(1)*n(1)
	mod2=l(2)*l(2)+m(2)*m(2)+n(2)*n(2)
	mod3=l(3)*l(3)+m(3)*m(3)+n(3)*n(3)
	dx1=l(1)*l(2)+m(1)*m(2)+n(1)*n(2)
	dx2=l(3)*l(2)+m(3)*m(2)+n(3)*n(2)
	dx3=l(1)*l(3)+m(1)*m(3)+n(1)*n(3)
		write(6,2002) x*lx,y*ly,z*lz,phi1,phi2,phi3
!		write(6,2003)l(1),m(1),n(1),l(2),m(2),n(2),l(3),m(3),n(3)
!		write(6,2004)	mod1,mod2,mod3,dx1,dx2,dx3
  2002   format(6(1x,f11.6))
  2003   format(9(1x,f11.6))
  2004   format(6(1x,f11.4))
	end do 
end program random

subroutine init_random_seed()
	integer :: i,n,clock
	integer,dimension(:),allocatable :: seed
	call random_seed(size=n)
	allocate(seed(n))
	call system_clock(count=clock)
	seed=clock+37*(/(i-1,i=1,n)/)
	call random_seed(PUT=seed)
	deallocate(seed)
end subroutine
