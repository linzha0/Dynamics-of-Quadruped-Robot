%Derivations for Kinematics of big dog for second phase
clc;

close all;
clear all;
format compact

syms R2x R2y th2 R3x R3y th3 R4x R4y th4 R5x R5y th5 R6x R6y th6 R7x R7y th7 R8x R8y th8 R9x R9y th9 R10x R10y th10
syms R2xd R2yd th2d R3xd R3yd th3d R4xd R4yd th4d R5xd R5yd th5d R6xd R6yd th6d R7xd R7yd th7d R8xd R8yd th8d R9xd R9yd th9d R10xd R10yd th10d
syms g m2 J2 m3 J3 m4 J4 m5 J5 m6 J6 m7 J7 m8 J8 m9 J9 m10 J10   
syms L2 L3 L4 L5 L6 L7 L8 L9 L10 M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11 M12 M13

%define variables and their time differentiation
q=[R2x;R2y;th2;R3x;R3y;th3;R4x;R4y;th4;R5x;R5y;th5;R6x;R6y;th6;R7x;R7y;th7;R8x;R8y;th8;R9x;R9y;th9;R10x;R10y;th10]; 
q_dot=[R2xd;R2yd;th2d;R3xd;R3yd;th3d;R4xd;R4yd;th4d;R5xd;R5yd;th5d;R6xd;R6yd;th6d;R7xd;R7yd;th7d;R8xd;R8yd;th8d;R9xd;R9yd;th9d;R10xd;R10yd;th10d];

%define independent and dependent variables
qi=[th3;R4x;th5;th7;th8;th9;th10]; 
qi_dot=[th3d;R4xd;th5d;th7d;th8d;th9d;th10d]; 
qd=[R2x;R2y;th2;R3x;R3y;R4y;th4;R5x;R5y;R6x;R6y;th6;R7x;R7y;R8x;R8y;R9x;R9y;R10x;R10y]; 
qd_dot=[R2xd;R2yd;th2d;R3xd;R3yd;R4yd;th4d;R5xd;R5yd;R6xd;R6yd;th6d;R7xd;R7yd;R8xd;R8yd;R9xd;R9yd;R10xd;R10yd];

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
    R6y+(L6/2)*sin(th6);
    R7x+(L7/2)*cos(th7)-R8x+(L8/2)*cos(th8);
    R7y+(L7/2)*sin(th7)-R8y+(L8/2)*sin(th8);
    R8x+(L8/2)*cos(th8)-R4x+(L4/2)*cos(th4);
    R8y+(L8/2)*sin(th8)-R4y+(L4/2)*sin(th4);
    R4x+(L4/2)*cos(th4)-R9x+(L9/2)*cos(th9);
    R4y+(L4/2)*sin(th4)-R9y+(L9/2)*sin(th9);
    R9x+(L9/2)*cos(th9)-R10x+(L10/2)*cos(th10);
    R9y+(L9/2)*sin(th9)-R10y+(L10/2)*sin(th10)]

%Jacobian matrix 
Cq=jacobian(C,q)

%Split Jacobian to independent and deoendent submatrices
Cqi=[Cq(:,6) Cq(:,7) Cq(:,12) Cq(:,18) Cq(:,21) Cq(:,24) Cq(:,27)] 
Cqd=[Cq(:,1:5) Cq(:,8:11) Cq(:,13:17) Cq(:,19:20) Cq(:,22:23) Cq(:,25:26)]

%Bid relates the dependent velocity to indeopndent ones

Bid=(-inv(Cqd)*Cqi)    
Bi=[Bid(1:5,:); 
    1 0 0 0 0 0 0;
    0 1 0 0 0 0 0;
    Bid(6:9,:);
    0 0 1 0 0 0 0;
    Bid(10:14,:);
    0 0 0 1 0 0 0;
    Bid(15:16,:);
    0 0 0 0 1 0 0;
    Bid(17:18,:);
    0 0 0 0 0 1 0;
    Bid(19:20,:);
    0 0 0 0 0 0 1]
%Generate 2nd order acceleration terms
%multiply Cqd by dependent velocity matrix 
Cqd_qdot=Cqd*qd_dot

%dieffervtiate w.r.t. dependent variables 
Cqd_qdot_qd=simplify(jacobian(Cqd_qdot,qd))

%repeat the process for the indepndent variable
Cqi_qdot=Cqi*qi_dot 
Cqi_qdot_qi=simplify(jacobian(Cqi_qdot,qi))











