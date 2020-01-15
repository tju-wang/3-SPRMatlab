clear
clc
syms pi

syms alpiha beta gama fai x0 y0 z0 ux uy uz vx vy vz wx wy wz
syms k1_x k1_y k1_z k2_x k2_y k2_z k3_x k3_y k3_z  
syms a b k   %结构参数

% Rz_alpiha = [cos(alpiha),-sin(alpiha),0;sin(alpiha),cos(alpiha),0;0,0,1];
% Rx_beta = [1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
% Rz_gama = [cos(gama),-sin(gama),0;sin(gama),cos(gama),0;0,0,1];
% R = Rz_alpiha*Rx_beta*Rz_gama

R = [ux,vx,wx;uy,vy,wy;uz,vz,wz]
B1 = b*[cos(0);sin(0);0];
B2 = b*[cos(pi/3);sin(pi/3);0];
B3 = b*[cos(-pi/3);sin(-pi/3);0];
B = [B1,B2,B3]      %B 是 b(i)值 

A1 = a*R*[cos(0);sin(0);0];
A2 = a*R*[cos(pi/3);sin(pi/3);0];
A3 = a*R*[cos(-pi/3);sin(-pi/3);0];
A = [A1,A2,A3]      %A是 a(i)值

K_1 = [k1_x;k1_y;k1_z];
K_2 = [k2_x;k2_y;k2_z];
K_3 = [k3_x;k3_y;k3_z];
K_K = [K_1,K_2,K_3]     %偏距 k
K_K = 0
OBA = [x0 y0 z0].'   %由OB指向BA的列矢量
A_A = [A(:,1)+OBA,A(:,2)+OBA,A(:,3)+OBA]
B_B = B
C_C = B_B + K_K

CA = A_A-B_B
%求 ci
c1 = R*[cos(pi/2),sin(pi/2),0].';
c2 = R*[cos(pi/6),sin(pi/6),0].'
c3 = R*[cos(5*pi/6),sin(5*pi/6),0].';
Ci = [c1,c2,c3] 

CA(:,1).'
c1
J1 = CA(:,1).'.* c1.'
[x0,y0,z0] = solve(J1)
J2 = CA(:,2).'.*c2.'
J3 = CA(:,3).'.*c3.'
[x0,y0,z0] = solve(J2)
[x0,y0,z0] = solve(J3)





















