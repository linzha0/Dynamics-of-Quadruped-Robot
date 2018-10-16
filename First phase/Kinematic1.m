%Derivations for Kinematics of RRRRRR Manipulator
clc;

close all;
clear all;
format compact

syms R2x R2y th2 R3x R3y th3 R4x R4y th4 R5x R5y th5 R6x R6y th6 L2 L3 L4 L5 L6
syms R2xd R2yd th2d R3xd R3yd th3d R4xd R4yd th4d R5xd R5yd th5d R6xd R6yd th6d
syms g m2 J2 m3 J3 m4 J4 m5 J5 m6 J6 M1 M2 M3 M4



%define variables and their time differentiation
q=[R2x;R2y;th2;R3x;R3y;th3;R4x;R4y;th4;R5x;R5y;th5;R6x;R6y;th6]; 
q_dot=[R2xd;R2yd;th2d;R3xd;R3yd;th3d;R4xd;R4yd;th4d;R5xd;R5yd;th5d;R6xd;R6yd;th6d];

%define independent and dependent variables
qi=[th3;R4x;th5]; 
qi_dot=[th3d;R4xd;th5d]; 
qd=[R2x;R2y;th2;R3x;R3y;R4y;th4;R5x;R5y;R6x;R6y;th6]; 
qd_dot=[R2xd;R2yd;th2d;R3xd;R3yd;R4yd;th4d;R5xd;R5yd;R6xd;R6yd;th6d];

%Kinematic constraint matrix
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
    R6y+(L6/2)*sin(th6)]
%Jacobian matrix 
Cq=jacobian(C,q)

%Split Jacobian to independent and deoendent submatrices
Cqi=[Cq(:,6) Cq(:,7) Cq(:,12)] 
Cqd=[Cq(:,1:5) Cq(:,8:11) Cq(:,13:15)]

%Bid relates the dependent velocity to indeopndent ones

Bid=(-inv(Cqd)*Cqi)    
Bi=[Bid(1:5,:); 
    1 0 0;
    0 1 0;
    Bid(6:9,:);
    0 0 1;
    Bid(10:12,:)]    
  
  
%Generate 2nd order acceleration terms
%multiply Cqd by dependent velocity matrix 
Cqd_qdot=Cqd*qd_dot

%dieffervtiate w.r.t. dependent variables 
Cqd_qdot_qd=simplify(jacobian(Cqd_qdot,qd))

%repeat the process for the indepndent variable
Cqi_qdot=Cqi*qi_dot 
Cqi_qdot_qi=simplify(jacobian(Cqi_qdot,qi))

% gamma_id=(inv(Cqd)*(Cqd_qdot_qd*Bid*qi_dot+Cqi_qdot_qi*qi_dot))
% 
% %expand gamma_id to include the indepndent variable coeffiecient at the
% gamma_i=[gamma_id(1:5,:);
%     0;
%     0;
%     gamma_id(6:9,:);
%     0;
%     gamma_id(10:12,:)]
% 
% %mass matrix (15x15)
% M=[m2 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%     0 m2 0 0 0 0 0 0 0 0 0 0 0 0 0;
%     0 0 J2 0 0 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 m3 0 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 0 m3 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 0 0 J3 0 0 0 0 0 0 0 0 0; 
%     0 0 0 0 0 0 m4 0 0 0 0 0 0 0 0;
%     0 0 0 0 0 0 0 m4 0 0 0 0 0 0 0; 
%     0 0 0 0 0 0 0 0 J4 0 0 0 0 0 0;
%     0 0 0 0 0 0 0 0 0 m5 0 0 0 0 0;
%     0 0 0 0 0 0 0 0 0 0 m5 0 0 0 0;
%     0 0 0 0 0 0 0 0 0 0 0 J5 0 0 0;
%     0 0 0 0 0 0 0 0 0 0 0 0 m6 0 0;
%     0 0 0 0 0 0 0 0 0 0 0 0 0 m6 0;
%     0 0 0 0 0 0 0 00 0 0 0 0 0 0 J6]
% 
% %External force matrix
% Qe=[0;
%     -m2*g;
%     M1; 
%     0;
%     -m3*g;
%     M2;
%     0
%     -m4*g; 
%     0;
%     0;
%     -m5*g;
%     M3;
%     0;
%     -m6*g;
%     M4]
% 
% %Generalized mass matrix
% Mi=[transpose(Bi)*M*Bi]
% 
% %Generalized force matrix
% Qi=transpose(Bi)*(Qe-M*gamma_i)
% 
% qidd=inv(Mi)*Qi







