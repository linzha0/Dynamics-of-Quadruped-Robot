
function [F]=FRRRRRR(y,yd,t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global g m2 J2 m3 J3 m4 J4 m5 J5 m6 J6 M1 M2 M3 M4 M9 M10 M13
global L2 L3 L4 L5 L6
global M Bi Qe gamma_i

%in put
th3=y(1);
R4x=y(2);
th5=y(3);
th3d=y(4);
R4xd=y(5);
th5d=y(6);

%label dependent and independent variables

R2x=yd(1,:);
R2y=yd(2,:);
th2=yd(3,:);
R3x=yd(4,:);
R3y=yd(5,:);
th3=y(1,:);
R4x=y(2,:);
R4y=yd(6,:);
th4=yd(7,:);
R5x=yd(8,:);
R5y=yd(9,:);
th5=y(3,:);
R6x=yd(10,:);
R6y=yd(11,:);
th6=yd(12,:);



%define variables and their time differentiation


qi=[th3;R4x;th5];
qi_dot=[th3d;R4xd;th5d];



Cqi =...
[                0,  0,                0;
                0,  0,                0;
 -(L3*sin(th3))/2,  0,                0;
  (L3*cos(th3))/2,  0,                0;
 -(L3*sin(th3))/2, -1,                0;
  (L3*cos(th3))/2,  0,                0;
                0,  1, -(L5*sin(th5))/2;
                0,  0,  (L5*cos(th5))/2;
                0,  0, -(L5*sin(th5))/2;
                0,  0,  (L5*cos(th5))/2;
                0,  0,                0;
                0,  0,                0];

Cqd=...
[ 1, 0,  (L2*sin(th2))/2,  0,  0,  0,                0,  0,  0,  0,  0,                0;
 0, 1, -(L2*cos(th2))/2,  0,  0,  0,                0,  0,  0,  0,  0,                0;
 1, 0, -(L2*sin(th2))/2, -1,  0,  0,                0,  0,  0,  0,  0,                0;
 0, 1,  (L2*cos(th2))/2,  0, -1,  0,                0,  0,  0,  0,  0,                0;
 0, 0,                0,  1,  0,  0, -(L4*sin(th4))/2,  0,  0,  0,  0,                0;
 0, 0,                0,  0,  1, -1,  (L4*cos(th4))/2,  0,  0,  0,  0,                0;
 0, 0,                0,  0,  0,  0, -(L4*sin(th4))/2, -1,  0,  0,  0,                0;
 0, 0,                0,  0,  0,  1,  (L4*cos(th4))/2,  0, -1,  0,  0,                0;
 0, 0,                0,  0,  0,  0,                0,  1,  0, -1,  0, -(L6*sin(th6))/2;
 0, 0,                0,  0,  0,  0,                0,  0,  1,  0, -1,  (L6*cos(th6))/2;
 0, 0,                0,  0,  0,  0,                0,  0,  0,  1,  0, -(L6*sin(th6))/2;
 0, 0,                0,  0,  0,  0,                0,  0,  0,  0,  1,  (L6*cos(th6))/2];


Bid=-inv(Cqd)*Cqi;

Bi=[Bid(1:5,:); 
    1 0 0;
    0 1 0;
    Bid(6:9,:);
    0 0 1;
    Bid(10:12,:)];

qd_dot=Bid*qi_dot;

th2d=qd_dot(3);
th4d=qd_dot(7);
th6d=qd_dot(12);

Cqd_qdot_qd =...
[ 0, 0,  (L2*th2d*cos(th2))/2, 0, 0, 0,                     0, 0, 0, 0, 0,                     0;
 0, 0,  (L2*th2d*sin(th2))/2, 0, 0, 0,                     0, 0, 0, 0, 0,                     0;
 0, 0, -(L2*th2d*cos(th2))/2, 0, 0, 0,                     0, 0, 0, 0, 0,                     0;
 0, 0, -(L2*th2d*sin(th2))/2, 0, 0, 0,                     0, 0, 0, 0, 0,                     0;
 0, 0,                     0, 0, 0, 0, -(L4*th4d*cos(th4))/2, 0, 0, 0, 0,                     0;
 0, 0,                     0, 0, 0, 0, -(L4*th4d*sin(th4))/2, 0, 0, 0, 0,                     0;
 0, 0,                     0, 0, 0, 0, -(L4*th4d*cos(th4))/2, 0, 0, 0, 0,                     0;
 0, 0,                     0, 0, 0, 0, -(L4*th4d*sin(th4))/2, 0, 0, 0, 0,                     0;
 0, 0,                     0, 0, 0, 0,                     0, 0, 0, 0, 0, -(L6*th6d*cos(th6))/2;
 0, 0,                     0, 0, 0, 0,                     0, 0, 0, 0, 0, -(L6*th6d*sin(th6))/2;
 0, 0,                     0, 0, 0, 0,                     0, 0, 0, 0, 0, -(L6*th6d*cos(th6))/2;
 0, 0,                     0, 0, 0, 0,                     0, 0, 0, 0, 0, -(L6*th6d*sin(th6))/2];

Cqi_qdot_qi =...
[                     0, 0,                     0;
                     0, 0,                     0;
 -(L3*th3d*cos(th3))/2, 0,                     0;
 -(L3*th3d*sin(th3))/2, 0,                     0;
 -(L3*th3d*cos(th3))/2, 0,                     0;
 -(L3*th3d*sin(th3))/2, 0,                     0;
                     0, 0, -(L5*th5d*cos(th5))/2;
                     0, 0, -(L5*th5d*sin(th5))/2;
                     0, 0, -(L5*th5d*cos(th5))/2;
                     0, 0, -(L5*th5d*sin(th5))/2;
                     0, 0,                     0;
                     0, 0,                     0];
gamma_id=inv(Cqd)*(Cqd_qdot_qd*Bid*qi_dot+Cqi_qdot_qi*qi_dot);

gamma_i=[gamma_id(1:5,:);
    0;
    0;
    gamma_id(6:9,:);
    0;
    gamma_id(10:12,:)];

M=[m2 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 m2 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 J2 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 m3 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 m3 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 J3 0 0 0 0 0 0 0 0 0; 
    0 0 0 0 0 0 m4 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 m4 0 0 0 0 0 0 0; 
    0 0 0 0 0 0 0 0 J4 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 m5 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 m5 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 J5 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 m6 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 m6 0;
    0 0 0 0 0 0 0 00 0 0 0 0 0 0 J6];

Qe=[0;
    -m2*g;
    M1; 
    0;
    -m3*g;
    M2+M9;
    0
    -m4*g; 
    M13;
    0;
    -m5*g;
    M3+M10;
    0;
    -m6*g;
    M4];

Mi=[transpose(Bi)*M*Bi];

Qi=transpose(Bi)*(Qe-M*gamma_i);

qidd=inv(Mi)*Qi;

F=[th3d;R4xd;th5d;qidd];
end













