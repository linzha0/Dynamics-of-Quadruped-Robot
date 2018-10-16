%use Runge-Kutta to simulate the walk status 

close all
clear all

global g m2 J2 m3 J3 m4 J4 m5 J5 m6 J6 m7 J7 m8 J8 m9 J9 m10 J10 
global L2 L3 L4 L5 L6 L7 L8 L9 L10
global M Bi Cq Qe gamma_i M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11 M12 M13

%define time step
h=0.001;

%define number of time steps
imax=100;

%define masses and lengths
L2=0.24;
L3=0.21;
L4=0.65;
L5=0.21;
L6=0.24;
L7=0.24;
L8=0.21;
L9=0.21;
L10=0.24;

m2=1.25;
m3=0.625;
m4=28.5;
m5=0.625;
m6=1.25;
m7=1.25;
m8=0.625;
m9=0.625;
m10=1.25;

J2=0.0092;
J3=0.0023;
J4=1.06;
J5=0.0023;
J6=0.0092;
J7=0.0092;
J8=0.0023;
J9=0.0023;
J10=0.0092;

g=9.81;

%define initial conditions for independent variables

y(1,1)=0.5593;
y(2,1)=0.3250;
y(3,1)=5.7238;
y(4,1)=2.4064;
y(5,1)=0.5593;
y(6,1)=5.7238;
y(7,1)=3.8768;
y(8,1)=0;
y(9,1)=0;
y(10,1)=0;
y(11,1)=0;
y(12,1)=0;
y(13,1)=0;
y(14,1)=0;

%define initial conditions for dependent varibales
%link 2
yd(3,1)=2.4064;
yd(1,1)=(L2/2)*cos(yd(3,1));
yd(2,1)=(L2/2)*sin(yd(3,1));

%link 3
yd(4,1)=(L2)*cos(yd(3,1))+(L3/2)*cos(y(1,1));
yd(5,1)=(L2)*sin(yd(3,1))+(L3/2)*sin(y(1,1));

%link 4
yd(7,1)=0;
yd(6,1)=(L2)*sin(yd(3,1))+(L3)*sin(y(1,1))+(L4/2)*sin(yd(7,1));

%link 5
yd(8,1)=(L2)*cos(yd(3,1))+(L3)*cos(y(1,1))+L4*cos(yd(7,1))+(L5/2)*cos(y(3,1));
yd(9,1)=(L2)*sin(yd(3,1))+(L3)*sin(y(1,1))+L4*sin(yd(7,1))+(L5/2)*sin(y(3,1));

%link 6
yd(12,1)=3.8768;
yd(10,1)=(L2)*cos(yd(3,1))+(L3)*cos(y(1,1))+L4*cos(yd(7,1))+L5*cos(y(3,1))+(L6/2)*cos(yd(12,1));
yd(11,1)=(L2)*sin(yd(3,1))+(L3)*sin(y(1,1))+L4*sin(yd(7,1))+L5*sin(y(3,1))+(L6/2)*sin(yd(12,1));

%link 7
yd(13,1)=(L2/2)*cos(yd(3,1));
yd(14,1)=(L2/2)*sin(yd(3,1));

%link 8
yd(15,1)=(L2)*cos(yd(3,1))+(L3/2)*cos(y(1,1));
yd(16,1)=(L2)*sin(yd(3,1))+(L3/2)*sin(y(1,1));

%link 9
yd(17,1)=(L2)*cos(yd(3,1))+(L3)*cos(y(1,1))+L4*cos(yd(7,1))+(L5/2)*cos(y(3,1));
yd(18,1)=(L2)*sin(yd(3,1))+(L3)*sin(y(1,1))+L4*sin(yd(7,1))+(L5/2)*sin(y(3,1));

%link 10
yd(19,1)=(L2)*cos(yd(3,1))+(L3)*cos(y(1,1))+L4*cos(yd(7,1))+L5*cos(y(3,1))+(L6/2)*cos(yd(12,1));
yd(20,1)=(L2)*sin(yd(3,1))+(L3)*sin(y(1,1))+L4*sin(yd(7,1))+L5*sin(y(3,1))+(L6/2)*sin(yd(12,1));

% inputs of torque

M9=100;
M10=100;
M11=0;
M12=-50;
M13=165;



%solve the differential equations
for i=1:imax
    t(i)=(i-1)*h;
    
    M1=-50*(yd(3,i)-y(1,i)-pi/6);
    M2=50*(yd(3,i)-y(1,i)-pi/6);
    M3=-50*(y(3,i)-yd(12,i)-pi/6);
    M4=50*(y(3,i)-yd(12,i)-pi/6); 
    M5=-50*(y(4,i)-y(5,i)-pi/6);
    M6=50*(y(4,i)-y(5,i)-pi/6);
    M7=-50*(y(6,i)-y(7,i)-pi/6);
    M8=50*(y(6,i)-y(7,i)-pi/6);

 
    f1=h*FR(y(:,i),yd(:,i),t(i));
    f2=h*FR(y(:,i)+0.5*f1,yd(:,i),t(i)+0.5*h);
    f3=h*FR(y(:,i)+0.5*f2,yd(:,i),t(i)+0.5*h);
    f4=h*FR(y(:,i)+f3,yd(:,i),t(i)+h);
    y(:,i+1)=y(:,i)+(f1+2*f2+2*f3+f4)/6;

    R2x=yd(1,i);
    R2y=yd(2,i);
    th2=yd(3,i);
    R3x=yd(4,i);
    R3y=yd(5,i);
    th3=y(1,i+1);
    R4x=y(2,i+1);
    R4y=yd(6,i);
    th4=yd(7,i); 
    R5x=yd(8,i);
    R5y=yd(9,i);
    th5=y(3,i+1);
    R6x=yd(10,i);
    R6y=yd(11,i);
    th6=yd(12,i);
    R7x=yd(13,i);
    R7y=yd(14,i);
    th7=y(4,i+1);
    R8x=yd(15,i);
    R8y=yd(16,i);
    th8=y(5,i+1);
    R9x=yd(17,i);
    R9y=yd(18,i);
    th9=y(6,i+1);
    R10x=yd(19,i);
    R10y=yd(20,i);
    th10=y(7,i+1);
    
    C=[R2x-(L2/2)*cos(th2);
       R2y-(L2/2)*sin(th2);
       R2x+(L2/2)*cos(th2)-R3x+(L3/2)*cos(th3);
       R2y+(L2/2)*sin(th2)-R3y+(L3/2)*sin(th3);
       R3x+(L3/2)*cos(th3)-R4x+(L4/2)*cos(th4);
       R3y+(L3/2)*sin(th3)-R4y+(L4/2)*sin(th4);
       R4x+(L4/2)*cos(th4)-R5x+(L5/2)*cos(th5);
       R4y+(L4/2)*sin(th4)-R5y+(L5/2)*sin(th5);
       R5x+(L5/2)*cos(th5)-R6x+(L6/2)*cos(th6);
       R5y+(L5/2)*sin(th5)-R6y+(L6/2)*sin(th6);
       R6x+(L6/2)*cos(th6)-L4;
       R6y+(L6/2)*sin(th6);
       R7x+(L7/2)*cos(th7)-R8x+(L8/2)*cos(th8);
       R7y+(L7/2)*sin(th7)-R8y+(L8/2)*sin(th8);
       R8x+(L8/2)*cos(th8)-R4x+(L4/2)*cos(th4);
       R8y+(L8/2)*sin(th8)-R4y+(L4/2)*sin(th4);
       R4x+(L4/2)*cos(th4)-R9x+(L9/2)*cos(th9);
       R4y+(L4/2)*sin(th4)-R9y+(L9/2)*sin(th9);
       R9x+(L9/2)*cos(th9)-R10x+(L10/2)*cos(th10);
       R9y+(L9/2)*sin(th9)-R10y+(L10/2)*sin(th10)];
   
   
   Cqi =...
[                0,  0,                0,                0,                0,                0,                  0;
                0,  0,                0,                0,                0,                0,                  0;
 -(L3*sin(th3))/2,  0,                0,                0,                0,                0,                  0;
  (L3*cos(th3))/2,  0,                0,                0,                0,                0,                  0;
 -(L3*sin(th3))/2, -1,                0,                0,                0,                0,                  0;
  (L3*cos(th3))/2,  0,                0,                0,                0,                0,                  0;
                0,  1, -(L5*sin(th5))/2,                0,                0,                0,                  0;
                0,  0,  (L5*cos(th5))/2,                0,                0,                0,                  0;
                0,  0, -(L5*sin(th5))/2,                0,                0,                0,                  0;
                0,  0,  (L5*cos(th5))/2,                0,                0,                0,                  0;
                0,  0,                0,                0,                0,                0,                  0;
                0,  0,                0,                0,                0,                0,                  0;
                0,  0,                0, -(L7*sin(th7))/2, -(L8*sin(th8))/2,                0,                  0;
                0,  0,                0,  (L7*cos(th7))/2,  (L8*cos(th8))/2,                0,                  0;
                0, -1,                0,                0, -(L8*sin(th8))/2,                0,                  0;
                0,  0,                0,                0,  (L8*cos(th8))/2,                0,                  0;
                0,  1,                0,                0,                0, -(L9*sin(th9))/2,                  0;
                0,  0,                0,                0,                0,  (L9*cos(th9))/2,                  0;
                0,  0,                0,                0,                0, -(L9*sin(th9))/2, -(L10*sin(th10))/2;
                0,  0,                0,                0,                0,  (L9*cos(th9))/2,  (L10*cos(th10))/2];
Cqd =...
[ 1, 0,  (L2*sin(th2))/2,  0,  0,  0,                0,  0,  0,  0,  0,                0, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 1, -(L2*cos(th2))/2,  0,  0,  0,                0,  0,  0,  0,  0,                0, 0, 0,  0,  0,  0,  0,  0,  0;
 1, 0, -(L2*sin(th2))/2, -1,  0,  0,                0,  0,  0,  0,  0,                0, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 1,  (L2*cos(th2))/2,  0, -1,  0,                0,  0,  0,  0,  0,                0, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 0,                0,  1,  0,  0, -(L4*sin(th4))/2,  0,  0,  0,  0,                0, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 0,                0,  0,  1, -1,  (L4*cos(th4))/2,  0,  0,  0,  0,                0, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 0,                0,  0,  0,  0, -(L4*sin(th4))/2, -1,  0,  0,  0,                0, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 0,                0,  0,  0,  1,  (L4*cos(th4))/2,  0, -1,  0,  0,                0, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 0,                0,  0,  0,  0,                0,  1,  0, -1,  0, -(L6*sin(th6))/2, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 0,                0,  0,  0,  0,                0,  0,  1,  0, -1,  (L6*cos(th6))/2, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 0,                0,  0,  0,  0,                0,  0,  0,  1,  0, -(L6*sin(th6))/2, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 0,                0,  0,  0,  0,                0,  0,  0,  0,  1,  (L6*cos(th6))/2, 0, 0,  0,  0,  0,  0,  0,  0;
 0, 0,                0,  0,  0,  0,                0,  0,  0,  0,  0,                0, 1, 0, -1,  0,  0,  0,  0,  0;
 0, 0,                0,  0,  0,  0,                0,  0,  0,  0,  0,                0, 0, 1,  0, -1,  0,  0,  0,  0;
 0, 0,                0,  0,  0,  0, -(L4*sin(th4))/2,  0,  0,  0,  0,                0, 0, 0,  1,  0,  0,  0,  0,  0;
 0, 0,                0,  0,  0, -1,  (L4*cos(th4))/2,  0,  0,  0,  0,                0, 0, 0,  0,  1,  0,  0,  0,  0;
 0, 0,                0,  0,  0,  0, -(L4*sin(th4))/2,  0,  0,  0,  0,                0, 0, 0,  0,  0, -1,  0,  0,  0;
 0, 0,                0,  0,  0,  1,  (L4*cos(th4))/2,  0,  0,  0,  0,                0, 0, 0,  0,  0,  0, -1,  0,  0;
 0, 0,                0,  0,  0,  0,                0,  0,  0,  0,  0,                0, 0, 0,  0,  0,  1,  0, -1,  0;
 0, 0,                0,  0,  0,  0,                0,  0,  0,  0,  0,                0, 0, 0,  0,  0,  0,  1,  0, -1];

      %Generate the new values of the dependent variables
    yd(:,i+1)=-inv(Cqd)*(C+Cqi*(y(1:7,i+1)-y(1:7,i)))+yd(:,i);  


end
t(i+1)=t(i)+h;

%calculate joint locations

xa=0;
ya=0;
xb=yd(1,:)+(L2/2)*cos(yd(3,:));
yb=yd(2,:)+(L2/2)*sin(yd(3,:));
xc=yd(4,:)+(L3/2)*cos(y(1,:));
yc=yd(5,:)+(L3/2)*sin(y(1,:));
xe=y(2,:)+(L4/2)*cos(yd(7,:));
ye=yd(6,:)+(L4/2)*sin(yd(7,:));
xf=yd(8,:)+(L5/2)*cos(y(3,:));
yf=yd(9,:)+(L5/2)*sin(y(3,:));
xg=yd(13,:)+(L7/2)*cos(y(4,:));
yg=yd(14,:)+(L7/2)*sin(y(4,:));
xm=yd(13,:)-(L7/2)*cos(y(4,:));
ym=yd(14,:)-(L7/2)*sin(y(4,:));
xh=yd(17,:)+(L9/2)*cos(y(6,:));
yh=yd(18,:)+(L9/2)*sin(y(6,:));
xn=yd(19,:)+(L10/2)*cos(y(7,:));
yn=yd(20,:)+(L10/2)*sin(y(7,:));

%plot the links
figure
for i=1:imax
    %generate the links at every instant
    u(1,i)=0;
    v(1,i)=0;
    u(2,i)=xb(i);
    v(2,i)=yb(i);
    u(3,i)=xc(i);
    v(3,i)=yc(i);
    u(4,i)=xe(i);
    v(4,i)=ye(i);
    u(5,i)=xf(i);
    v(5,i)=yf(i);
    u(6,i)=0.65;
    v(6,i)=0;
    u(7,i)=xg(i);
    v(7,i)=yg(i);
    u(8,i)=xm(i);
    v(8,i)=ym(i);
    u(9,i)=xh(i);
    v(9,i)=yh(i);
    u(10,i)=xn(i);
    v(10,i)=yn(i);
    
    U1=[u(7,i);
        u(3,i)];
    V1=[v(7,i);
        v(3,i)];
    U2=[u(4,i);
        u(9,i)];
    V2=[v(4,i);
        v(9,i)];    
    
    
    
   if(i/4)-floor(i/4)==0
    plot(u(1:2,i),v(1:2,i),'b',u(2:3,i),v(2:3,i),'r',u(3:4,i),v(3:4,i),'g',u(4:5,i),v(4:5,i),'r',u(5:6,i),v(5:6,i),'b',u(7:8,i),v(7:8,i),'y',U1,V1,'k',U2,V2,'k',u(9:10,i),v(9:10,i),'y','LineWidth',2)
   end
   hold on
    
    
end

title('motion of th links')
axis equal
grid on











