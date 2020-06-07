%% 论文 利用3-SPR逆解显式关系 正解迭代
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
[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3))

%旧的坐标系
% A1_o = [0;-a;0]
% A2_o = [(3.^(1/2)/2)*a;1/2*a;0]
% A3_o = [-(3^(1/2)/2)*a;1/2*a;0]
% 
% B1 = [0;-b;0]
% B2 = [(3^(1/2)/2)*b;1/2*b;0]
% B3 = [-(3^(1/2)/2)*b;1/2*b;0]


%修改坐标系  x轴正方向向  右   y正方向向上   z 正方向由静平台指向动平台
A1_o = [0;a;0];
A2_o = [-(3.^(1/2)/2)*a;-1/2*a;0];
A3_o = [(3^(1/2)/2)*a;-1/2*a;0];

B1 = [0;b;0];
B2 = [-(3^(1/2)/2)*b;-1/2*b;0];
B3 = [(3^(1/2)/2)*b;-1/2*b;0];

% X0 = (b*uy*(ux-3*vy)+2*Z0*wx)/(2*wz)
% Y0 = (b*ux*(vy-ux)+2*b*vx*uy+2*Z0*wy)/(2*wz)

% X0 =  (b*uy*(3*vy-ux)+2*Z0*wx)/(2*wz)    %论文结果 
% Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*wy)/(2*wz)
X0 =  (b*uy*(3*vy-ux)+2*Z0*(vz*uy-vy*uz))/(2*(ux*vy-vx*uy))   %论文结果 
Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*(uz*vx-vz*ux))/(2*(ux*vy-vx*uy))

% X0 = 2*Z0*wx/(2*wz)
% Y0 = 2*Z0*wy/(2*wz)

alpha = gama

Ao = [X0;Y0;Z0]
A1 = R*A1_o + Ao
A2 = R*A2_o + Ao
A3 = R*A3_o + Ao


f = [(B1-A1)'*(B1-A1);(B2-A2)'*(B2-A2);(B3-A3)'*(B3-A3)]
argu = [alpha beta Z0]
J = []
J = jacobian(f,argu)

% J_1 = inv(J)
% 开始数值运算  结构常数等进行赋值
s = 62  %测量得 58.5mm + 6.5/2
alpha = 0.1;
beta = 0.1;
gama = alpha;
Z0 = 250;
%q1 =287;   q2 = 309;   q3 = 322;
%q1 =261;   q2 =261;   q3 = 261;

 pq = [273.87 261.49 258.50];
 q1 = pq(1);  q2 = pq(2);  q3 = pq(3);


a = 41.56 %动平台外接圆半径
b = 80
argu = eval(argu)

err = 1;  err_2 = 10;
Fi = eval([-q1^2+(B1-A1)'*(B1-A1);-q2^2+(B2-A2)'*(B2-A2);-q3^2+(B3-A3)'*(B3-A3)])
%  diff_argu' = - inv(eval(J))*Fi;
diff_argu = [0 0 0]
num = 0
while(err>1.0e-4 || err_2>1.0e-4 )
   argu = (argu + diff_argu);
   alpha = argu(1);
   beta = argu(2);
   Z0 = argu(3);
   gama = alpha;
   Fi = eval([-q1^2+(B1-A1)'*(B1-A1);-q2^2+(B2-A2)'*(B2-A2);-q3^2+(B3-A3)'*(B3-A3)]);
   diff_argu = (-inv(eval(J))*Fi)';
   err = diff_argu(1)*diff_argu(1)+diff_argu(2)*diff_argu(2)
   err_2 = (diff_argu(3)*diff_argu(3))
   num = num+1
end
num
argu
err
X0 = eval(X0)
Y0 = eval(Y0)

