clear all
clc
syms alpha beta  gama X0 Y0 Z0 
syms a b ux uy uz vx vy vz wx wy wz s %a b 分别为静动平台三角形外接圆半径  s为机构偏距
RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)]
% RX_gama = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)]

%欧拉角
R = RX_alpha*RY_beta*RX_gama
A1_o = [-(3.^(1/2)/2)*a;1/2*a;0]
A2_o = [(3^(1/2)/2)*a;1/2*a;0]
A3_o = [0;-a;0]

B1 = [-(3^(1/2)/2)*b;1/2*b;0]
B2 = [(3^(1/2)/2)*b;1/2*b;0]
B3 = [0;-b;0]

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3))
X0 = (b*uy*(3*vy-ux)+2*Z0*wx)/(2*wz)
Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*wy)/(2*wz)

alpha = gama

Ao = [X0;Y0;Z0]
A1 = R*A1_o + Ao
A2 = R*A2_o + Ao
A3 = R*A3_o + Ao

% q1 = sqrt((B1-A1)'*(B1-A1)-s.^2)    %q表示 P副长度 + 球铰连接柱长度（两者始终平行）
% q2 = sqrt((B2-A2)'*(B2-A2)-s.^2)
% q3 = sqrt((B3-A3)'*(B3-A3)-s.^2)

f = [(B1-A1)'*(B1-A1)-s.^2;(B2-A2)'*(B2-A2)-s.^2;(B3-A3)'*(B3-A3)-s.^2]
argu = [alpha beta Z0]
J = []
J = jacobian(f,argu)
%计算C1 C2 C3坐标
syms c1X c1Y c1Z K1x K1y K1z
C1 = [c1X c1Y c1Z]
K1 = [K1x K1y K1z] %K1 为A1处的转动副轴线向量（方向）
K1x = 


