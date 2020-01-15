clear all
clc
tic  %��ʱ��ʼ
syms alpha beta  gama X0 Y0 Z0
syms a b ux uy uz vx vy vz wx wy wz s  %a b �ֱ�Ϊ����ƽ̨���������Բ�뾶  sΪ����ƫ��
global gX0 gY0 gZ0;
RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)];

%ŷ����
R = RX_alpha*RY_beta*RX_gama;
[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3))

%�޸�����ϵ  x����������  ��   y����������   z �������ɾ�ƽָ̨��ƽ̨
A1_o = [0;a;0];
A2_o = [-(3.^(1/2)/2)*a;-1/2*a;0];
A3_o = [(3^(1/2)/2)*a;-1/2*a;0];

B1 = [0;b;0];
B2 = [-(3^(1/2)/2)*b;-1/2*b;0];
B3 = [(3^(1/2)/2)*b;-1/2*b;0];

X0 = (b*uy*(3*vy-ux)+2*Z0*(vz*uy-vy*uz))/(2*(ux*vy-vx*uy))   %���Ľ�� 
Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*(uz*vx-vz*ux))/(2*(ux*vy-vx*uy))
alpha = gama
gX0 = -10;  gY0 = -10; gZ0 = 250;
fun = @solve_fun;
x0 = [0.1,0.1,0.1];
[ang_argu] = fsolve(fun,x0)

