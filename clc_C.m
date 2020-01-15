clear
clc

syms alpha beta  gama X0 Y0 Z0
syms a b ux uy uz vx vy vz wx wy wz s %a b 分别为静动平台三角形外接圆半径  s为机构偏距
syms Xc1 Xc2 Xc3 Yc1 Yc2 Yc3 Zc1 Zc2 Zc3
RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)]
% RX_gama = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)]

%欧拉角
R = RX_alpha*RY_beta*RX_gama
A2_o = [-(3.^(1/2)/2)*a;1/2*a;0]
A3_o = [(3^(1/2)/2)*a;1/2*a;0]
A1_o = [0;-a;0]

B2 = [-(3^(1/2)/2)*b;1/2*b;0]
B3 = [(3^(1/2)/2)*b;1/2*b;0]
B1 = [0;-b;0]

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3))
% X0 = (b*uy*(ux-3*vy)+2*Z0*wx)/(2*wz)
% Y0 = (b*ux*(vy-ux)+2*b*vx*uy+2*Z0*wy)/(2*wz)

alpha = gama

Ao = [X0;Y0;Z0]
A1 = R*A1_o + Ao
A2 = R*A2_o + Ao
A3 = R*A3_o + Ao


C1 = [Xc1,Yc1,Zc1]
C2 = [Xc2,Yc2,Zc2]
C3 = [Xc3,Yc3,Zc3]

vA1C1 = vector(A1,C1) 
vA2C2 = vector(A2,C2)
vA3C3 = vector(A3,C3)

vB1C1 = vector(B1,C1) 
vB2C2 = vector(B2,C2)
vB3C3 = vector(B3,C3)

q1 = 287;  q2 = 323; q3 = 309;

norm(vA1C1) = q1;
norm(vA2C2) = q2;
norm(vA3C3) = q3;



