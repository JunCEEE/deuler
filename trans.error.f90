program trans
  implicit none
  double precision pi,lat,b0(3,3)
  double precision x,y,z,phi1,phi2,phi3
  double precision cosa1,cosa2,cosa3
  double precision cosb1,cosb2,cosb3
  double precision cosr1,cosr2,cosr3
  double precision,allocatable,dimension(:) :: thida,phi,kosin
  double precision sint,sinab,cosab
  double precision cos1,cos2,sin1,sin2
  integer,allocatable,dimension(:) :: BC
  integer i, ngrain

    pi=atan(1.0)*4.0
    open(101,file='input.txt')
    open(102,file='ebsd.dat')   

    read(101,*)lat
    read(101,*)ngrain,b0(1,1),b0(2,2),b0(3,3)
    read(101,*)

    allocate(thida(ngrain),phi(ngrain),kosin(ngrain),BC(ngrain))

    do i=1,ngrain
     read(101,*)x,y,z,phi1,phi2,phi3
     cosa1=cos(phi2)*sin(phi1)
     cosb1=sin(phi2)*sin(phi1)
     cosr1=cos(phi1)
     cosa2=(cosa1*cosa1+(cosa2*cosa2+cosa3*cosa3)*cos(phi3))*(-1.0)*sin(phi2)+(cosa1*cosa2*(1-cos(phi3))-cosa3*sin(phi3))*cos(phi2);
     cosb2=(cosa1*cosa2*(1-cos(phi3))+cosa3*sin(phi3))*(-1.0)*sin(phi2)+(cosa2*cosa2+(cosa1*cosa1+cosa3*cosa3)*cos(phi3))*cos(phi2);
     cosr2=(cosa1*cosa3*(1-cos(phi3))-cosa2*sin(phi3))*(-1.0)*sin(phi2)+(cosa2*cosa3*(1-cos(phi3))+cosa1*sin(phi3))*cos(phi2);
     cosa3=cosa2*cosb3-cosa3*cosb2;
     cosb3=-1.0*cosa1*cosb3+cosa3*cosb1;
     cosr3=cosa1*cosb2-cosa2*cosb1;

     if(cosa1>1.0) cosa1=1.0; if(cosa1<-1.0) cosa1=-1.0;
     if(cosa2>1.0) cosa2=1.0; if(cosa2<-1.0) cosa2=-1.0;
     if(cosa3>1.0) cosa3=1.0; if(cosa3<-1.0) cosa3=-1.0;
     if(cosb1>1.0) cosb1=1.0; if(cosb1<-1.0) cosb1=-1.0;
     if(cosb2>1.0) cosb2=1.0; if(cosb2<-1.0) cosb2=-1.0;
     if(cosb3>1.0) cosb3=1.0; if(cosb3<-1.0) cosb3=-1.0; 
     if(cosr1>1.0) cosr1=1.0; if(cosr1<-1.0) cosr1=-1.0;
     if(cosr2>1.0) cosr2=1.0; if(cosr2<-1.0) cosr2=-1.0;
     if(cosr3>1.0) cosr3=1.0; if(cosr3<-1.0) cosr3=-1.0; 

       thida(i)=acos(cosr3)
       sint=sin(thida(i))
      if(abs(sint)<5e-3) then

       cosab=cosa1
       sinab=cosb1
        if(sinab>0.0 .and. cosab>0.0) then
         phi(i)=0.0
         kosin(i)=acos(abs(cosa1))
        endif
        if(sinab>0.0 .and. cosab<0.0) then
         phi(i)=0.0
         kosin(i)=pi-acos(abs(cosa1))
        endif
        if(sinab<0.0 .and. cosab<0.0) then
         phi(i)=0.0
         kosin(i)=pi+acos(abs(cosa1))
        endif
        if(sinab<0.0 .and. cosab>0.0) then
         phi(i)=0.0
         kosin(i)=pi*2.0-acos(abs(cosa1))
        endif
        if(abs(sinab)<5e-3) then
         if(cosab>0.0) then
          phi(i)=0.0
          kosin(i)=0.0
         else
          phi(i)=0.0
          kosin(i)=pi
         endif
        endif
       if(abs(cosab)<5e-3) then
        if(sinab>0.0) then
         phi(i)=0.0
         kosin(i)=0.5*pi
        else
         phi(i)=0.0
         kosin(i)=1.5*pi
        endif
       endif

      else

       if(abs(cosa3/sint)<1.0) then
         phi(i)=asin(abs(cosa3/sint))
         kosin(i)=asin(abs(cosr1/sint))
         sin1=cosa3/sint
         cos1=-cosb3/sint
         sin2=cosr1/sint
         cos2=cosr2/sint
         if(sin1>0.0 .and. cos1>0.0) phi(i)=phi(i)
         if(sin1>0.0 .and. cos1<0.0) phi(i)=pi-phi(i)
         if(sin1<0.0 .and. cos1<0.0) phi(i)=pi+phi(i)
         if(sin1<0.0 .and. cos1>0.0) phi(i)=2.0*pi-phi(i)
         if(abs(cos1)<5e-3) then
          if(sin1>0.0) then
           phi(i)=0.5*pi
          else
           phi(i)=1.5*pi
          endif
         endif
         if(abs(sin1)<5e-3) then
          if(cos1>0.0) then
           phi(i)=0.0
          else
           phi(i)=pi
          endif
         endif
         if(sin2>0.0 .and. cos2>0.0) kosin(i)=kosin(i)
         if(sin2>0.0 .and. cos2<0.0) kosin(i)=pi-kosin(i)
         if(sin2<0.0 .and. cos2<0.0) kosin(i)=pi+kosin(i)
         if(sin2<0.0 .and. cos2>0.0) kosin(i)=2.0*pi-kosin(i)       
         if(abs(cos2)<5e-3) then
          if(sin2>0.0) then
           kosin(i)=0.5*pi
          else
           kosin(i)=1.5*pi
          endif
         endif
         if(abs(sin2)<5e-3) then
          if(cos2>0.0) then
           kosin(i)=0.0
          else
           kosin(i)=pi
          endif
         endif
       endif

       if(abs(cosa3/sint)>1.0) then
         if(cosa3*sint>0.0) phi(i)=0.5*pi
         if(cosa3*sint<0.0) phi(i)=-0.5*pi
       endif
       if(abs(cosr1/sint)>1.0) then
         if(cosr1*sint>0.0) kosin(i)=0.5*pi
         if(cosr1*sint<0.0) kosin(i)=1.5*pi  
       endif 
       if(abs(abs(cosa3/sint)-1.0)<5e-3) then
         if(cosa3*sint>0.0) phi(i)=0.5*pi
         if(cosa3*sint<0.0) phi(i)=1.5*pi
       endif
       if(abs(abs(cosr1/sint)-1.0)<5e-3) then
         if(cosr1*sint>0.0) kosin(i)=0.5*pi
         if(cosr1*sint<0.0) kosin(i)=1.5*pi
       endif

      endif           

    thida(i)=thida(i)*180/pi
    phi(i)=phi(i)*180/pi
    kosin(i)=kosin(i)*180/pi
    BC(i)=255

     write(102,2001)x,y,z,phi(i),thida(i),kosin(i),BC(i)
2001 format('0',6(1x,f10.5),1x,i5)

    enddo
        close(101)
        close(102)
end program
