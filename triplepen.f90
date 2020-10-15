program Doublependulum


real::dt,x1,y1,t0,tf,x_10,y_10,t,x_20,y_20,x2,y2,om1,om2,om_10,om_20,x3,y3,x_30,y_30,&
      f1,f2,g1,g2,k_11,k_21,k_31,k_41,l_11,l_21,l_31,l_41,k_12,k_22,k_32,k_42,l_12,l_22,l_32,l_42,&
      theta_10,theta_20,rad1,rad2,thetarad_10,thetarad_20,theta_30,om_30,rad3,thetarad_30,omega_10,omega_20,omega_30,om3,&
      k_13,l_13,k_23,l_23,k_33,l_33,k_43,l_43
real,parameter::pi=acos(-1.0),l1=3.0,l2=3.0,l3=3.0,xor=0.0,yor=0.0,g=9.81,m1=2.0,m2=2.0,m3=2.0
integer::n,i



print*,"give the initial time and final time"
read*,t0,tf
print*,"give the number of interval"
read*,n
!print*,"give the length of the pendulum "
!read*,R
print*,"give the initial value of theta1(in deg),theta2,theta3,omega1,omega2 and omega3 "
read*,theta_10,theta_20,theta_30,omega_10,omega_20,omega_30



thetarad_10=(pi*theta_10/180)
thetarad_20=(pi*theta_20/180)
thetarad_30=(pi*theta_30/180)

open(30,file="threepen.dat")
!open(35,file="modeuler")
dt=(tf-t0)/n
print*,"dt=",dt
om1=omega_10
om2=omega_20
om3=omega_30
rad1=thetarad_10
rad2=thetarad_20
rad3=thetarad_30

t=t0

x_10=xor+l1*sin(rad1)
y_10=yor-l1*cos(rad1)

x_20=x_10+l2*sin(rad2)
y_20=y_10-l2*cos(rad2)

x_30=x_20+l3*sin(rad3)
y_30=y_20-l3*cos(rad3)

x1=x_10
y1=y_10

x2=x_20
y2=y_20

x3=x_30
y3=y_30

do i=0,n
!==================================================================================================
k_11=dt*f1(rad1,rad2,rad3,om1,om2,om3,t)
l_11=dt*g1(rad1,rad2,rad3,om1,om2,om3,t)
k_21=dt*f1(rad1+k_11*0.5,rad2+k_11*0.5,rad3+k_11*0.5,om1+l_11*0.5,om2+l_11*0.5,om3+l_11*0.5,t+dt*0.5)
l_21=dt*g1(rad1+k_11*0.5,rad2+k_11*0.5,rad3+k_11*0.5,om1+l_11*0.5,om2+l_11*0.5,om3+l_11*0.5,t+dt*0.5)
k_31=dt*f1(rad1+k_21*0.5,rad2+k_21*0.5,rad3+k_21*0.5,om1+l_21*0.5,om2+l_21*0.5,om3+l_21*0.5,t+dt*0.5)
l_31=dt*g1(rad1+k_21*0.5,rad2+k_21*0.5,rad3+k_21*0.5,om1+l_21*0.5,om2+l_21*0.5,om3+l_21*0.5,t+dt*0.5)
k_41=dt*f1(rad1+k_31*0.5,rad2+k_31*0.5,rad3+k_31*0.5,om1+l_31*0.5,om2+l_31*0.5,om3+l_31*0.5,t+dt*0.5)
l_41=dt*g1(rad1+k_31*0.5,rad2+k_31*0.5,rad3+k_31*0.5,om1+l_31*0.5,om2+l_31*0.5,om3+l_31*0.5,t+dt*0.5)
!===================================================================================================

k_12=dt*f2(rad1,rad2,rad3,om1,om2,om3,t)
l_12=dt*g2(rad1,rad2,rad3,om1,om2,om3,t)
k_22=dt*f2(rad1+k_12*0.5,rad2+k_12*0.5,rad3+k_12*0.5,om1+l_12*0.5,om2+l_12*0.5,om3+l_12*0.5,t+dt*0.5)
l_22=dt*g2(rad1+k_12*0.5,rad2+k_12*0.5,rad3+k_12*0.5,om1+l_12*0.5,om2+l_12*0.5,om3+l_12*0.5,t+dt*0.5)
k_32=dt*f2(rad1+k_22*0.5,rad2+k_22*0.5,rad3+k_22*0.5,om1+l_22*0.5,om2+l_22*0.5,om3+l_22*0.5,t+dt*0.5)
l_32=dt*g2(rad1+k_22*0.5,rad2+k_22*0.5,rad3+k_22*0.5,om1+l_22*0.5,om2+l_22*0.5,om3+l_22*0.5,t+dt*0.5)
k_42=dt*f2(rad1+k_32*0.5,rad2+k_32*0.5,rad3+k_32*0.5,om1+l_32*0.5,om2+l_32*0.5,om3+l_32*0.5,t+dt*0.5)
l_42=dt*g2(rad1+k_32*0.5,rad2+k_32*0.5,rad3+k_32*0.5,om1+l_32*0.5,om2+l_32*0.5,om3+l_32*0.5,t+dt*0.5)
!========================================================================================================

k_13=dt*f3(rad1,rad2,rad3,om1,om2,om3,t)
l_13=dt*g3(rad1,rad2,rad3,om1,om2,om3,t)
k_23=dt*f3(rad1+k_13*0.5,rad2+k_13*0.5,rad3+k_13*0.5,om1+l_13*0.5,om2+l_13*0.5,om3+l_13*0.5,t+dt*0.5)
l_23=dt*g3(rad1+k_13*0.5,rad2+k_13*0.5,rad3+k_13*0.5,om1+l_13*0.5,om2+l_13*0.5,om3+l_13*0.5,t+dt*0.5)
k_33=dt*f3(rad1+k_23*0.5,rad2+k_23*0.5,rad3+k_23*0.5,om1+l_23*0.5,om2+l_23*0.5,om3+l_23*0.5,t+dt*0.5)
l_33=dt*g3(rad1+k_23*0.5,rad2+k_23*0.5,rad3+k_23*0.5,om1+l_23*0.5,om2+l_23*0.5,om3+l_23*0.5,t+dt*0.5)
k_43=dt*f3(rad1+k_33*0.5,rad2+k_33*0.5,rad3+k_33*0.5,om1+l_33*0.5,om2+l_33*0.5,om3+l_33*0.5,t+dt*0.5)
l_43=dt*g3(rad1+k_33*0.5,rad2+k_33*0.5,rad3+k_33*0.5,om1+l_33*0.5,om2+l_33*0.5,om3+l_33*0.5,t+dt*0.5)

!=========================================================================================================




write(30,*)0,0,x1,y1,x2,y2,x3,y3,rad1,om1,rad2,om2
!===================================================================================================
t=i*dt+t0
!========================================================================
rad1=rad1+(k_11+2*k_21+2*k_31+k_41)/6.0
om1=om1+(l_11+2*l_11+2*l_31+l_41)/6.0
!=========================================================================
rad2=rad2+(k_12+2*k_22+2*k_32+k_42)/6.0
om2=om2+(l_12+2*l_12+2*l_32+l_42)/6.0
!=========================================================================
rad3=rad3+(k_13+2*k_23+2*k_33+k_43)/6.0
om3=om3+(l_13+2*l_13+2*l_33+l_43)/6.0
!=========================================================================
x1=xor+l1*sin(rad1)
y1=yor-l1*cos(rad1)

x2=xor+l1*sin(rad1)+l2*sin(rad2)
y2=yor-l1*cos(rad1)-l2*cos(rad2)

x3=x2+l3*sin(rad3)
y3=y2-l3*cos(rad3)

!=======================================================================
!print*,x,y
!write(30,*)x,y,t

end do

end program
!close(35)
!call system ('gnuplot -p euler_plot.plt')
!call system ('gnuplot -p modeuler_plot.plt')

!---------------------------------------------------------------------
real function f1(rad1,rad2,rad3,om1,om2,om3,t)
implicit none

real,intent(in)::rad1,rad2,rad3,om1,om2,om3,t
f1=om1
end function f1
!=======================================================================================================

real function f2(rad1,rad2,rad3,om1,om2,om3,t)
implicit none

real,intent(in)::rad1,rad2,rad3,om1,om2,om3,t
f2=om2
end function f2

!========================================================================================================

real function f3(rad1,rad2,rad3,om1,om2,om3,t)
implicit none

real,intent(in)::rad1,rad2,rad3,om1,om2,om3,t
f3=om3
end function f3
!========================================================================================================
real function g1(rad1,rad2,rad3,om1,om2,om3,t)
implicit none
real,intent(in)::rad1,rad2,rad3,om1,om2,om3,t
real,parameter::pi=acos(-1.0),l1=3.0,l2=3.0,l3=3.0,xor=0.0,yor=0.0,g=9.81,m1=2.0,m2=2.0,m3=2.0

g1=(2*rad2-3*rad1)*g/l1


end function g1
!============================================++++++++======================================================
real function g2(rad1,rad2,rad3,om1,om2,om3,t)
implicit none
real,intent(in)::rad1,rad2,rad3,om1,om2,om3,t
real,parameter::pi=acos(-1.0),l1=3.0,l2=3.0,l3=3.0,xor=0.0,yor=0.0,g=9.81,m1=2.0,m2=2.0,m3=2.0

g2=(3*rad1-4*rad2+rad3)*g/l1

end function g2
!===================================================================================================================================
real function g3(rad1,rad2,rad3,om1,om2,om3,t)
implicit none
real,intent(in)::rad1,rad2,rad3,om1,om2,om3,t
real,parameter::pi=acos(-1.0),l1=3.0,l2=3.0,l3=3.0,xor=0.0,yor=0.0,g=9.81,m1=2.0,m2=2.0,m3=2.0

g3=2*(rad2-rad3)*g/l1


end function g3
