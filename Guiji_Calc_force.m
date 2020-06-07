%解算机构重力补偿 电机力
%条件  q1-q3   求解Fm1 - Fm3
%步骤1.求解3-SPR机构正解  求末端X Y Z坐标点  2.求解三个球铰受力(重力等效) 
%3.力平衡方程 求解P副受力  4.比较受力相对大小  乘系数 确定机构补偿重力的PWM上下限

%相较于motor_force.m文件  修正了B点坐标  修正求B点力矩 
%% 1.
%符号运算
clear all
clc
syms alpha beta  gama X0 Y0 Z0 real
syms a b ux uy uz vx vy vz wx wy wz s real%a b 分别为静动平台三角形外接圆半径  s为机构偏距
%C点坐标 P副长度q
syms Xc1 Xc2 Xc3 Yc1 Yc2 Yc3 Zc1 Zc2 Zc3 q1 q2 q3 real
RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
% RX_gama = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
% RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)];
RZ_gama = [cos(gama),-sin(gama),0;sin(gama),cos(gama),0;0,0,1];
%欧拉角
R = RX_alpha*RY_beta*RZ_gama;

%修改坐标系  x轴正方向向  右   y正方向向上   z 正方向由静平台指向动平台
A1_o = [0;a;0];
A2_o = [-(3.^(1/2)/2)*a;-1/2*a;0];
A3_o = [(3^(1/2)/2)*a;-1/2*a;0];

B1 = [0;b;20];
B2 = [-(3^(1/2)/2)*b;-1/2*b;20];
B3 = [(3^(1/2)/2)*b;-1/2*b;20];

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3));

Ao = [X0;Y0;Z0];
A1 = R*A1_o + Ao;
A2 = R*A2_o + Ao;
A3 = R*A3_o + Ao;
%定义用到的6条支链
r1 = A2' - A3';
r1 = r1/norm(r1);
r2 = A1' - A3';
r2 = r2/norm(r2);
r3 = A1' - A2';
r3 = r3/norm(r3);
L = 40;     %假设的R副长度的一半是20mm
A11 = A1' + r1*L;
A12 = A1' - r1*L;
A21 = A2' + r2*L;
A22 = A2' - r2*L;
A31 = A3' + r3*L;
A32 = A3' - r3*L;

%表达6条支链
syms q11 q12 q21 q22 q31 q32 q1 q2 q3 real

q11 = sqrt(q1^2 + L^2);   q12 = q11;
q21 = sqrt(q2^2 + L^2);   q22 = q21;
q31 = sqrt(q3^2 + L^2);   q32 = q31;

BA = [ B1'-A11;B1'-A12; B2'-A21;B2'-A22; B3'-A31;B3'-A32];
Fi = [ BA(1,:)*BA(1,:)'-q11^2;BA(2,:)*BA(2,:)'-q12^2; BA(3,:)*BA(3,:)'-q21^2;BA(4,:)*BA(4,:)'-q22^2; BA(5,:)*BA(5,:)'-q31^2;BA(6,:)*BA(6,:)'-q32^2 ];
Fii =[ BA(1,:)*BA(1,:)'-q11^2;BA(2,:)*BA(2,:)'-q12^2; BA(3,:)*BA(3,:)'-q21^2;BA(4,:)*BA(4,:)'-q22^2; BA(5,:)*BA(5,:)'-q31^2;BA(6,:)*BA(6,:)'-q32^2 ];

argu = [alpha beta gama X0 Y0 Z0];
J = [];
J = jacobian(Fi,argu);

%% 2.求C点坐标 球铰力的表达 符号运算
% 3.力平衡方程 求解P副受力   先求偏距与P副交点C1 C2 C3
%求C点的条件为 VBC vAC垂直 vAC模长已知 vBC与R副的轴线vAxi 垂直
C1 = [Xc1,Yc1,Zc1];
C2 = [Xc2,Yc2,Zc2];
C3 = [Xc3,Yc3,Zc3];

vA1C1 = vector(A1',C1);     % A1C1直线的向量
vA2C2 = vector(A2',C2);
vA3C3 = vector(A3',C3);

vB1C1 = vector(B1,C1);
vB2C2 = vector(B2,C2);
vB3C3 = vector(B3,C3);

vAxi_A1 = A2' - A3';
vAxi_A2 = A1' - A3';
vAxi_A3 = A1' - A2';
%% 3.数值运算 机构正解
% 开始数值运算  结构常数等进行赋值

syms Xc Yc Zc XB YB ZB XA YA ZA Xi Yi Zi real

s = 62;  %测量得 58.5mm + 6.5/2

a = 41.56; %动平台外接圆半径
b = 80;


load .\..\..\'3-SPR Kinematic'\MATLAB\calc_forceData\force.mat

alpha = 0.1;
beta = 0.15;
gama = 0.12;
Z0 = 280;
X0 = 10;
Y0 = 10;
argu = eval(argu);
Fi = eval(Fi);
for i = 1:size(result,1)
    q1 =result(i,2);   q2 =result(i,3);   q3 =result(i,4);

err = 1;

diff_argu = [0 0 0 0 0 0];
num = 0;
while(err>1.0e-6)
    argu = (argu + diff_argu);
    alpha = argu(1);
    beta = argu(2);
    gama = argu(3);
    X0 = argu(4);
    Y0 = argu(5);
    Z0 = argu(6);
    
    Fi = eval(Fii);
    diff_argu = (-inv(eval(J))*Fi)';
    err = diff_argu*diff_argu';
    num = num+1
end
num;
argu = (argu + diff_argu);
alpha = argu(1);
beta = argu(2);
gama = argu(3);
X0 = argu(4);
Y0 = argu(5);
Z0 = argu(6);
%argu

%解 C1 坐标

Xc = Xc1; Yc = Yc1; Zc = Zc1; XB = B1(1); YB = B1(2); ZB = B1(3); XA = A1(1); YA = A1(2); ZA = A1(3);
Xi = A2(1)-A3(1);  Yi = A2(2)-A3(2);  Zi = A2(3)-A3(3);

M1 = Xi*XB+Yi*YB+Zi*ZB;
K1 = 2*(YA-YB)-2*(XA-XB)*Yi/Xi;
K2 = 2*(ZA-ZB)-2*(XA-XB)*Zi/Xi;
M2 = 2*(XA-XB)*M1/Xi+(XB^2+YB^2+ZB^2-(XA^2+YA^2+ZA^2)+q1^2-2*s^2);
K4 = -K1/K2;
M4 = -M2/K2;
K3 = -(Zi*K4+Yi)/Xi;
M3 = (M1-Zi*M4)/Xi;

ca = eval(K3^2+K4^2+1);
cb = eval(-(2*XB*K3+2*YB+2*ZB*K4)+2*K3*M3+2*K4*M4);
cc = eval(M3^2+M4^2 -(2*XB*M3+2*ZB*M4)+XB^2+YB^2+ZB^2-s^2); 

if real((((cb.^2-4*ca*cc)))) > 0
    Yc = real(((-cb+sqrt(cb.^2-4*ca*cc))/2/ca));
    Xc = eval(K3*Yc+M3);
    Zc = eval(Yc*K4+M4);
% Xc = eval((M1-Yc*Yi-Zc*Zi)/Xi)
    Xc1 = Xc;
    Yc1 = Yc;
    Zc1 = Zc;

else
    sprintf('不满足方程求解条件')
end

fc1_1 = eval(norm(vB1C1)-s);
fc1_2 = eval(vB1C1 * vAxi_A1');
fc1_3 = eval(norm(vA1C1)-sqrt(q1^2-s^2));
fc1 = [fc1_1 fc1_2 fc1_3];
 
% 解C2坐标
%代入验证  解C2 C3
%解一元二次方程组 

XB = B2(1); YB = B2(2); ZB = B2(3); XA = A2(1); YA = A2(2); ZA = A2(3);
Xi = A1(1)-A3(1);  Yi = A1(2)-A3(2); Zi;
Zi = A1(3)-A3(3);

M1 = Xi*XB+Yi*YB+Zi*ZB;
K1 = 2*(YA-YB)-2*(XA-XB)*Yi/Xi;
K2 = 2*(ZA-ZB)-2*(XA-XB)*Zi/Xi;
M2 = 2*(XA-XB)*M1/Xi+(XB^2+YB^2+ZB^2-(XA^2+YA^2+ZA^2)+q2^2-2*s^2);
K4 = -K1/K2;
M4 = -M2/K2;
K3 = -(Zi*K4+Yi)/Xi;
M3 = (M1-Zi*M4)/Xi;

ca = eval(K3^2+K4^2+1);
cb = eval(-(2*XB*K3+2*YB+2*ZB*K4)+2*K3*M3+2*K4*M4);
cc = eval(M3^2+M4^2 -(2*XB*M3+2*ZB*M4)+XB^2+YB^2+ZB^2-s^2); 

if real((((cb.^2-4*ca*cc)))) > 0
    Yc = real(((-cb-sqrt(cb.^2-4*ca*cc))/2/ca));
Xc = eval(K3*Yc+M3);
Zc = eval(Yc*K4+M4);
% Xc = eval((M1-Yc*Yi-Zc*Zi)/Xi)
    Xc2 = Xc;
    Yc2 = Yc;
    Zc2 = Zc;
else
    sprintf('不满足方程求解条件')
end

fc2_1 = eval(norm(vB2C2)-s);
fc2_2 = eval(vB2C2 * vAxi_A2');
fc2_3 = eval(norm(vA2C2)-sqrt(q2^2-s^2));
fc2 = [fc2_1 fc2_2 fc2_3];


XB = B3(1); YB = B3(2); ZB = B3(3); XA = A3(1); YA = A3(2); ZA = A3(3);
Xi = A1(1)-A2(1);  Yi = A1(2)-A2(2); 
Zi = A1(3)-A2(3);

M1 = Xi*XB+Yi*YB+Zi*ZB;
K1 = 2*(YA-YB)-2*(XA-XB)*Yi/Xi;
K2 = 2*(ZA-ZB)-2*(XA-XB)*Zi/Xi;
M2 = 2*(XA-XB)*M1/Xi+(XB^2+YB^2+ZB^2-(XA^2+YA^2+ZA^2)+q3^2-2*s^2);
K4 = -K1/K2;
M4 = -M2/K2;
K3 = -(Zi*K4+Yi)/Xi;
M3 = (M1-Zi*M4)/Xi;

ca = eval(K3^2+K4^2+1);
cb = eval(-(2*XB*K3+2*YB+2*ZB*K4)+2*K3*M3+2*K4*M4);
cc = eval(M3^2+M4^2 -(2*XB*M3+2*ZB*M4)+XB^2+YB^2+ZB^2-s^2); 

if real((((cb.^2-4*ca*cc)))) > 0
    Yc = real(((-cb-sqrt(cb.^2-4*ca*cc))/2/ca));
%     Xc = real(eval((-K2*Yc-M2)/K1));
%     Zc = real(eval((M1-(Xi*Xc+Yi*Yc))/Zi));
Xc = eval(K3*Yc+M3);
Zc = eval(Yc*K4+M4);
% Xc = eval((M1-Yc*Yi-Zc*Zi)/Xi)
    Xc3 = Xc;
    Yc3 = Yc;
    Zc3 = Zc;
else
    sprintf('不满足方程求解条件')
end

fc3_1 = eval(norm(vB3C3)-s);
fc3_2 = eval(vB3C3 * vAxi_A3');
fc3_3 = eval(norm(vA3C3)-sqrt(q3^2-s^2));
fc3 = [fc3_1 fc3_2 fc3_3];

% 5.解质点的坐标
%表示各个质点 修正球铰受力  计算电机加的力  共10个质点  C语言 需要把计算力的部分 放到一块
%力的大小为常数  力的作用点为三个坐标
%10个点的坐标  凑齐 C点坐标 暂作为电机、连接块部分质心 需要修正 

syms FGy Gx Gy Gz real
syms Fx_b1 Fy_b1 Fz_b1 Fx_b2 Fy_b2 Fz_b2 Fx_b3 Fy_b3 Fz_b3 real    %球铰受力 底座给球铰的力
syms TG1_1 TG1_2 TG1_3 TG2_1 TG2_2 TG2_3 TG3_1 TG3_2 TG3_3 TGm real 
%TG1_1 扭矩-第一条支链 配重块 _2 导轨  _3 电机 滚子  Gm  末端

%定义各点力的大小  方向
FB1 = [Fx_b1,Fy_b1,Fz_b1];
FB2 = [Fx_b2,Fy_b2,Fz_b2];
FB3 = [Fx_b3,Fy_b3,Fz_b3];
FG = [0,FGy,0];
G = [Gx,Gy,Gz];
%  求力矩
B0 = [0;0;20];
TB1 = cross(B1-B0,FB1);
TB2 = cross(B2-B0,FB2);
TB3 = cross(B3-B0,FB3);
%解Fz 列绕x轴的力矩方程 TG(1) Fz_b1*b Fz_b2*b/2 Fz_b3*b/2 
%因为B1 B2 B3与原点在同一平面内 所以y方向没有对于绕x的扭矩
%定义常数
FG_1 = [0,-438*10,0];%单位g   力的单位 mN  配重块重量  clump weight
FG_2 = [0,-285.75*10,0]; %导轨重量  guide
FG_3 = [0,-342.5*10,0]; %BC间的链接架 Bracket 总质量为 329.1（连接架+球铰+导轨座） + 79（电机） = 355.3g
FGm = [0,-(86+30+30)*10,0];  %操作端  Moving platform
Q = 310;%导轨总长度 310mm
P1 = clc_P(eval(vA1C1),eval(C1),Q-sqrt(q1^2-s^2));  %P为配重块坐标
P2 = clc_P(eval(vA2C2),eval(C2),Q-sqrt(q2^2-s^2));
P3 = clc_P(eval(vA3C3),eval(C3),Q-sqrt(q3^2-s^2));

D1 = eval((P1+A1')/2);      %D为导轨重心
D2 = eval((P2+A2')/2);
D3 = eval((P3+A3')/2);
Gm = [X0 Y0 Z0];
%求所有物块重力  扭矩TG
TG1_1 = cross(P1-B0',FG_1);  %对B0点取距
TG2_1 = cross(P2-B0',FG_1);
TG3_1 = cross(P3-B0',FG_1);

TG1_2 = cross(D1-B0',FG_2);
TG2_2 = cross(D2-B0',FG_2);
TG3_2 = cross(D3-B0',FG_2);

TG1_3 = cross(C1-B0',FG_3); %求连接架对B0点的力矩   连接架的重心等效到C点
TG2_3 = cross(C2-B0',FG_3);
TG3_3 = cross(C3-B0',FG_3);
TGm = cross(Gm-B0',FGm);

TG = eval(TG1_1+TG2_1+TG3_1  +TG1_2+TG2_2+TG3_2  +TG1_3+TG2_3+TG3_3  +TGm);
FG = FG_1*3+FG_2*3+FG_3*3+FGm;
%要满足力平衡   力矩平衡  对x轴 力矩平衡得到Fz_b2  之后由对称关系  得到Fz_b3  由合力为0 得出Fz_b1
%求解B1 B2 B3处的受力情况，相同方向间，不同球铰B的受力比例关系是已知的  未知的是随末端位置变化的TG(1),即所有支链对x轴的矩
%会随着位置变化而变化  关键在这一步力矩平衡
Fz_b2 = solve(2*b*Fz_b2 + b/2*Fz_b2 + b/2*Fz_b2 - TG(1),Fz_b2);
Fz_b3 = Fz_b2;
Fz_b1 = -2*Fz_b2;
%解Fy 列重力方向力平衡方程  Fy_b2的大小与末端位置无关 这个地方？对不对啊啊？？
Fy_b2 = -(FG(2))/3;
Fy_b3 = Fy_b2;
Fy_b1 = Fy_b2;
%解Fx 列绕z轴的转矩方程 TG(3)+Fx_b1*3*b=0
Fx_b2 = -TG(3)/2/b;
Fx_b1 = -2*Fx_b2; %Fx_b2与Fx_b1方向相反
Fx_b3 = Fx_b2;
%F 的单位为 mN
FB1 = vpa(eval(FB1),8);
FB2 = vpa(eval(FB2),8);
FB3 = vpa(eval(FB3),8);

% 求解P副 主动关节需要补偿的力  列力平衡方程
% 原始思路 P副推力 与 部件所受重力 球铰受力做力平衡  不对！：FB1 FB2 FB3在叠加后只剩重力项

%在单个P副处做力平衡  球铰受力 连接块、电机重力 P副处的三个方向力平衡{但作用点不同 不知能否直接叠加}
n1 = eval((vA1C1)/norm((vA1C1)));  %AC方向 P副的力
if n1(3)>0
    n1 = -n1;
end
vn2 = A2'-A3';   %平行于R副转轴方向
n2 = eval(vn2/norm(vn2));
vn3 = cross(vA1C1,vn2);    %垂直于  AC与R副转轴所构成的平面
n3 = eval(vn3/norm(vn3));

N = ([n1;n2;n3])';
bn = (-eval(FB1+FG_3))';
F1 = N\bn;

clear n1 n2 n3 vn2 vn3 N bn
n1 = eval((vA2C2)/norm((vA2C2)));  %AC方向 P副的力
if n1(3)>0
    n1 = -n1;
end
vn2 = A1'-A3';   %平行于R副转轴方向
n2 = eval(vn2/norm(vn2));
vn3 = cross(vA2C2,vn2);    %垂直于  AC与R副转轴所构成的平面
n3 = eval(vn3/norm(vn3));
if n2(2)>0
    n2 = -n2;
end
if n3(1)>0
    n3 = -n3;
end

N = ([n1;n2;n3])';
bn = (-eval(FB2+FG_3))';
F2 = N\bn;

clear n1 n2 n3 vn2 vn3 N bn
n1 = eval((vA3C3)/norm((vA3C3)));  %AC方向 P副的力
if n1(3)>0
    n1 = -n1;
end
vn2 = A1'-A2';   %平行于R副转轴方向
n2 = eval(vn2/norm(vn2));
vn3 = cross(vA3C3,vn2);    %垂直于  AC与R副转轴所构成的平面
n3 = eval(vn3/norm(vn3));

if n2(2)>0
    n2 = -n2;
end
if n3(1)<0
    n3 = -n3;
end
N = ([n1;n2;n3])';
bn = (-eval(FB3+FG_3))';
F3 = N\bn;

FMotor(i,:) = [i F1(1)  F2(1)  F3(1) ];
FMotor(i,:)
end

%% 画图
figure(9)
hold on
plot(result(:,1),result(:,8),'r')
plot(result(:,1),result(:,9),'g')

figure(10)
hold on
plot(FMotor(:,1),FMotor(:,2))
plot(FMotor(:,1),FMotor(:,3))
plot(FMotor(:,1),FMotor(:,4))

