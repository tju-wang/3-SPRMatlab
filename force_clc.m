clear
clc
syms Xb1 Yb1 Zb1 Xb2 Yb2 Zb2 Xb3 Yb3 Zb3 Gx Gy Gz Xg1 Yg1 Zg1 G2x G2y G2z Xg12 Yg12 Zg12
syms Fx_b1 Fy_b1 Fz_b1 Fx_b2 Fy_b2 Fz_b2 Fx_b3 Fy_b3 Fz_b3 
%定义各点坐标
b = 80;
Xb1 = 0; Yb1 = -b; Zb1 = 0;
Xb2 = sqrt(3)/2*b; Yb2 = 1/2*b; Zb2 = 0;
Xb3 = -sqrt(3)/2*b;Yb3 = 1/2*b;Zb3 = 0;
Xg1 = 40; Yg1 = 10; Zg1 = 12;
Xg12 = -6; Yg2 = 3; Zg12 = 4;

Gx = 0; Gy = 60; Gz= 0;
G2x = 0; G2y = 0; G2z = 0;

B1 = [Xb1,Yb1,Zb1]; 
B2 = [Xb2,Yb2,Zb2];
B3 = [Xb3,Yb3,Zb3];
G1 = [Xg1,Yg1,Zg1];
G2 = [Xg12,Yg2,Zg12];
O = [0,0,0];
%定义各点力的大小  方向
FB1 = [Fx_b1,Fy_b1,Fz_b1];
FB2 = [Fx_b2,Fy_b2,Fz_b2];
FB3 = [Fx_b3,Fy_b3,Fz_b3];
FG1 = [Gx,Gy,Gz]
FG2 = [G2x,G2y,G2z];
FG = FG1 + FG2;

%求力矩
TG1 =  cross(G1,FG)
TG2 = cross(G2,FG2)
TG = TG1 + TG2;
TB1 = cross(B1,FB1);
TB2 = cross(B2,FB2);
TB3 = cross(B3,FB3)

Fz_b1 = -2*Fz_b2;
Fz_b2 = Fz_b3;
TB1(1)
TB2(1)
TB3(1)
TG(1)

Fz_b2 = solve(80*Fz_b2 + 40*Fz_b2 + 40*Fz_b2 + TG(1),Fz_b2)
Fz_b3 = Fz_b2
Fz_b1 = -2*Fz_b2

Fy_b2 = -(FG(2))/3;
Fy_b3 = Fy_b2;
Fy_b1 = Fy_b2;

Fx_b2 = TG(3)/3/2/b
Fx_b1 = -2*Fx_b2
Fx_b3 = Fx_b2

FB1 = eval(FB1)
FB2 = eval(FB2)
FB3 = eval(FB3)



% [A,b] = equationsToMatrix(M,[Fz_b1 Fz_b2 Fz_b3 ])



% 3.10
% % 对B1点出取距
% Rb1_G = clc_libi(B1,G);
% Rb1_B2 = clc_libi(B1,B2);
% Rb1_B3 = clc_libi(B1,B3);
% Rb1_G2 = clc_libi(B1,G2);
% 
% TG_B1 = clc_liju(FG,Rb1_G);
% TB2_B1 = clc_liju(FB2,Rb1_B2);
% TB3_B1 = clc_liju(FB3,Rb1_B3);
% TG2_B1 = clc_liju(G2,Rb1_G2);
% 
% %对B2点处取距
% Rb2_G = clc_libi(B2,G);
% Rb2_B1 = clc_libi(B2,B1);
% Rb2_B3 = clc_libi(B2,B3);
% Rb2_G2 = clc_libi(B2,G2);
% 
% TG_B2 = clc_liju(FG,Rb2_G);
% TB1_B2 = clc_liju(FB1,Rb2_B1);
% TB3_B2 = clc_liju(FB3,Rb2_B3);
% TG2_B2 = clc_liju(FG2,Rb2_G2);
% 
% %对B3点处取距
% Rb3_G = clc_libi(B3,G);
% Rb3_B1 = clc_libi(B3,B1);
% Rb3_B2 = clc_libi(B3,B2);
% Rb3_G2 = clc_libi(B3,G2);
% 
% TG_B3 = clc_liju(FG,Rb3_G);
% TB1_B3 = clc_liju(FB1,Rb3_B1);
% TB2_B3 = clc_liju(FB2,Rb3_B2);
% TG2_B3 = clc_liju(G2,Rb3_G2);
% %9*9 矩阵解方程
% M = [TG_B1(1)+TB2_B1(1)+TB3_B1(1)+TG2_B1(1);
%      TG_B1(2)+TB2_B1(2)+TB3_B1(2)+TG2_B1(2);
%      TG_B1(3)+TB2_B1(3)+TB3_B1(3)+TG2_B1(3);
%      TG_B2(1)+TB1_B2(1)+TB3_B2(1)+TG2_B2(1);
%      TG_B2(2)+TB1_B2(2)+TB3_B2(2)+TG2_B2(2);
%      TG_B2(3)+TB1_B2(3)+TB3_B2(3)+TG2_B2(3);
%      TG_B3(1)+TB1_B3(1)+TB2_B3(1)+TG2_B3(1);
%      TG_B3(2)+TB1_B3(2)+TB2_B3(2)+TG2_B3(2);
%      TG_B3(3)+TB1_B3(3)+TB2_B3(3)+TG2_B3(3);
% %      Fx_b1+Fx_b2+Fx_b3;
% %      Fy_b1+Fy_b2+Fy_b3+Gy;
% %      Fz_b1+Fz_b2+Fz_b3;
% ]
% %将等式中的系数提取出来 作为系数矩阵
% [A,b] = equationsToMatrix(M,[Fx_b1 Fy_b1 Fz_b1 Fx_b2 Fy_b2 Fz_b2 Fx_b3 Fy_b3 Fz_b3 ])
% 
% format rat
% X = vpa(rref([A b]))    %解非齐次矩阵




