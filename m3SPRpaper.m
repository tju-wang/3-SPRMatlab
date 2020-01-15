%解算机构重力补偿 电机力
%条件  q1-q3   求解Fm1 - Fm3
%步骤1.求解3-SPR机构正解  求末端X Y Z坐标点  2.求解三个球铰受力(重力等效) 
%3.力平衡方程 求解P副受力  4.比较受力相对大小  乘系数 确定机构补偿重力的PWM上下限

%相较于motor_force.m文件  修正了B点坐标  修正求B点力矩 
%开始整理论文后创建的文件
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
RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)];


RZ_alpha = [cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
RX_beta = [1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
RZ_gama = [cos(gama),-sin(gama),0;sin(gama),cos(gama),0;0,0,1];
%欧拉角
%R = RZ_alpha*RX_beta*RZ_gama  %相对于动坐标系  XYZ欧拉角
R = RX_alpha*RY_beta*RX_gama

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3))
%R  =[ux,vx,wx;uy,vy,wy;uz,vz,wz]

%修改坐标系  x轴正方向向  右   y正方向向上   z 正方向由静平台指向动平台
A1_o = [0;a;0];
A2_o = [-(3.^(1/2)/2)*a;-1/2*a;0];
A3_o = [(3^(1/2)/2)*a;-1/2*a;0];

%B1 = [0;b;20];
% B2 = [-(3^(1/2)/2)*b;-1/2*b;20];
% B3 = [(3^(1/2)/2)*b;-1/2*b;20];
B1 = [0;b;0];
B2 = [-(3^(1/2)/2)*b;-1/2*b;0];
B3 = [(3^(1/2)/2)*b;-1/2*b;0];

% X0 = (b*uy*(ux-3*vy)+2*Z0*wx)/(2*wz);
% Y0 = (b*ux*(vy-ux)+2*b*vx*uy+2*Z0*wy)/(2*wz);
X0 =  (b*uy*(3*vy-ux)+2*Z0*wx)/(2*wz);    %论文结果 
Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*wy)/(2*wz);
gama = alpha
% X0 = 2*Z0*wx/(2*wz)
% Y0 = 2*Z0*wy/(2*wz)

Ao = [X0;Y0;Z0];
A1 = R*A1_o + Ao;
A2 = R*A2_o + Ao;
A3 = R*A3_o + Ao;

%列向量  q1  q2  q3
BA = [ A1-B1,A2-B2,A3-B3]

C1 = R*(A2_o-A3_o)
C2 = R*(A1_o-A3_o)
C3 = R*(A1_o-A2_o)

C = [C1,C2,C3]

%求逆解的时候用  下边三个式子值为0  得出alpha = gama X0  Y0 的表达式
% BA(:,1)'*(C1)
% BA(:,2)'*(C2)
% BA(:,3)'*(C3)

%按照论文结果  表达q1  q2  q3
syms q1 q2 q3  qVar
%Fi处的BA 应该表达为关于 alpha beta  Z0的式子
%q1 = (BA(:,1)'*BA(:,1))
qVAr = []
Fi = [BA(:,1)'*BA(:,1)-q1^2;BA(:,1)'*BA(:,1)-q2^2; BA(:,1)'*BA(:,1)-q3^2 ]
Fii = [BA(:,1)'*BA(:,1)-q1^2;BA(:,1)'*BA(:,1)-q2^2; BA(:,1)'*BA(:,1)-q3^2 ]

argu = [alpha beta Z0];
J = []
J = jacobian(Fi,argu)

% Fi = [ BA(1,:)*BA(1,:)'-q11^2;BA(2,:)*BA(2,:)'-q12^2; BA(3,:)*BA(3,:)'-q21^2;BA(4,:)*BA(4,:)'-q22^2; BA(5,:)*BA(5,:)'-q31^2;BA(6,:)*BA(6,:)'-q32^2 ];
% Fii =[ BA(1,:)*BA(1,:)'-q11^2;BA(2,:)*BA(2,:)'-q12^2; BA(3,:)*BA(3,:)'-q21^2;BA(4,:)*BA(4,:)'-q22^2; BA(5,:)*BA(5,:)'-q31^2;BA(6,:)*BA(6,:)'-q32^2 ];


% Fi = [ BA(1,:)*BA(1,:)'-q11^2;BA(2,:)*BA(2,:)'-q12^2; BA(3,:)*BA(3,:)'-q21^2;BA(4,:)*BA(4,:)'-q22^2; BA(5,:)*BA(5,:)'-q31^2;BA(6,:)*BA(6,:)'-q32^2 ];
% Fii =[ BA(1,:)*BA(1,:)'-q11^2;BA(2,:)*BA(2,:)'-q12^2; BA(3,:)*BA(3,:)'-q21^2;BA(4,:)*BA(4,:)'-q22^2; BA(5,:)*BA(5,:)'-q31^2;BA(6,:)*BA(6,:)'-q32^2 ];
% 
% argu = [alpha beta gama X0 Y0 Z0];
% J = []
% J = jacobian(Fi,argu)

%% 正解 根据逆解结果计算
clear all
clc
syms alpha beta  gama X0 Y0 Z0 real
syms a b ux uy uz vx vy vz wx wy wz s real%a b 分别为静动平台三角形外接圆半径  s为机构偏距
%C点坐标 P副长度q
syms Xc1 Xc2 Xc3 Yc1 Yc2 Yc3 Zc1 Zc2 Zc3 q1 q2 q3 real

RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)];
% R = RX_alpha*RY_beta*RX_gama

R = RX_alpha*RY_beta*RX_alpha

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3))

%R = [ux,vx,wx;uy,vy,wy;uz,vz,wz]

%修改坐标系  x轴正方向向  右   y正方向向上   z 正方向由静平台指向动平台
A1_o = [0;a;0];
A2_o = [-(3.^(1/2)/2)*a;-1/2*a;0];
A3_o = [(3^(1/2)/2)*a;-1/2*a;0];

Ao = [X0;Y0;Z0];
A1 = R*A1_o + Ao;
A2 = R*A2_o + Ao;
A3 = R*A3_o + Ao;

%B1 = [0;b;20];
% B2 = [-(3^(1/2)/2)*b;-1/2*b;20];
% B3 = [(3^(1/2)/2)*b;-1/2*b;20];
B1 = [0;b;0];
B2 = [-(3^(1/2)/2)*b;-1/2*b;0];
B3 = [(3^(1/2)/2)*b;-1/2*b;0];

%逆解得出的结果   此时需要对欧拉角做处理  使得R=R = RX_alpha*RY_beta*RX_alpha   即直接将gama替换为alpha  这样才能做到X0  Y0 到最后的变量q1 q2 q3表达为 alpha beta 与Z0的式子
gama = alpha
X0 =  (b*uy*(3*vy-ux)+2*Z0*(vz*uy-vy*uz))/(2*(ux*vy-vx*uy))   %论文结果 
Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*(uz*vx-vz*ux))/(2*(ux*vy-vx*uy))

BA = [ A1-B1,A2-B2,A3-B3]
%表达变量qVar[]  1*3
qVar = [BA(:,1)'*BA(:,1),BA(:,2)'*BA(:,2),BA(:,3)'*BA(:,3)]
Fi = [qVar(1) - q1^2;qVar(2) - q2^2;qVar(3) - q3^2]
Fii = [qVar(1) - q1^2;qVar(2) - q2^2;qVar(3) - q3^2]

argu = [alpha beta Z0];
J = []
J = jacobian(Fi,argu)

% 开始数值运算  求正解
s = 62;  %测量得 58.5mm + 6.5/2
alpha = 0.1;
beta = 0.15;
% gama = 0.12;
Z0 = 280;

% X0 = 10;
% Y0 = 10;
%q1 =234.93;   q2 =275.6786;   q3 =289.5438;
% q1 =236.4224;   q2 =278.7179;   q3 =291.7681;
q1 =261;   q2 =261;   q3 = 261;

a = 41.56; %动平台外接圆半径
b = 80;
argu = eval(argu);

err = 1;
Fi = eval(Fi);
% Fi = eval([-q1^2+(B1-A1)'*(B1-A1)-s.^2;-q2^2+(B2-A2)'*(B2-A2)-s.^2;-q3^2+(B3-A3)'*(B3-A3)-s.^2])
%  diff_argu' = - inv(eval(J))*Fi;
diff_argu = [0 0 0];
num = 0;
while(err>1.0e-6)
    point = 1
    argu = (argu + diff_argu)
    alpha = (argu(1))
    beta = (argu(2));
    Z0 = (argu(3));
%     for i = 1:1:2        %限制角度范围  在+-pi之间
%        while(argu(i)<-pi)
%           argu(i) = argu(i)+pi; 
%        end
%        while(argu(i)>pi)
%           argu(i) = argu(i)-pi; 
%        end
%     end
    num
    Fi = eval(Fii);
    diff_argu = (-inv(eval(J))*Fi)';
    err = eval(diff_argu*diff_argu')
    num = num+1
end
num;
argu = (argu + diff_argu);
argu
X0
Y0


%% 正解 转化为6杆机构
clear all
clc
syms alpha beta  gama X0 Y0 Z0 real
syms a b ux uy uz vx vy vz wx wy wz s real%a b 分别为静动平台三角形外接圆半径  s为机构偏距
%C点坐标 P副长度q
syms Xc1 Xc2 Xc3 Yc1 Yc2 Yc3 Zc1 Zc2 Zc3 q1 q2 q3 real

RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
% RX_gama = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)];

RZ_alpha = [cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
RX_beta = [1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
RZ_gama = [cos(gama),-sin(gama),0;sin(gama),cos(gama),0;0,0,1];
%欧拉角
%R = RZ_alpha*RX_beta*RZ_gama  %相对于动坐标系  XYZ欧拉角
R = RX_alpha*RY_beta*RX_gama

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3));
%R  =[ux,vx,wx;uy,vy,wy;uz,vz,wz]

%修改坐标系  x轴正方向向  右   y正方向向上   z 正方向由静平台指向动平台
A1_o = [0;a;0];
A2_o = [-(3.^(1/2)/2)*a;-1/2*a;0];
A3_o = [(3^(1/2)/2)*a;-1/2*a;0];

%B1 = [0;b;20];
% B2 = [-(3^(1/2)/2)*b;-1/2*b;20];
% B3 = [(3^(1/2)/2)*b;-1/2*b;20];
B1 = [0;b;0];
B2 = [-(3^(1/2)/2)*b;-1/2*b;0];
B3 = [(3^(1/2)/2)*b;-1/2*b;0];

Ao = [X0;Y0;Z0];
A1 = R*A1_o + Ao;
A2 = R*A2_o + Ao;
A3 = R*A3_o + Ao;


C1 = R*(A2_o-A3_o);
C2 = R*(A1_o-A3_o);
C3 = R*(A1_o-A2_o);
C = [C1,C2,C3];

%定义6根杆
r1 = A2' - A3';
r1 = r1/norm(r1);
r2 = A1' - A3';
r2 = r2/norm(r2);
r3 = A1' - A2';
r3 = r3/norm(r3);
L = 100;     %假设的R副长度的一半是20mm
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

% 3.数值运算 机构正解
% 开始数值运算  结构常数等进行赋值
s = 62;  %测量得 58.5mm + 6.5/2
alpha = 0.3;
beta = -0.3;
gama = 0.3;
Z0 = 250;
X0 = 10;
Y0 = 10;
%q1 =234.93;   q2 =275.6786;   q3 =289.5438;
% q1 =236.4224;   q2 =278.7179;   q3 =291.7681;
%q1 =261;   q2 =261;   q3 = 261;
%q1 = 287;  q2 = 309;  q3 = 322;

pq = [290.0163  269.2276  255.1586];
q1 = pq(1);  q2 = pq(2);  q3 = pq(3);

a = 41.56; %动平台外接圆半径
b = 80;
argu = eval(argu);

err = 1;err_2 = 10;
Fi = eval(Fi);
% Fi = eval([-q1^2+(B1-A1)'*(B1-A1)-s.^2;-q2^2+(B2-A2)'*(B2-A2)-s.^2;-q3^2+(B3-A3)'*(B3-A3)-s.^2])
%  diff_argu' = - inv(eval(J))*Fi;
diff_argu = [0 0 0 0 0 0];
num = 0;
while(err>1.0e-6 || err_2>1.0e-6)
    argu = (argu + diff_argu);
    alpha = argu(1);
    beta = argu(2);
    gama = argu(3);
    X0 = argu(4);
    Y0 = argu(5);
    Z0 = argu(6);
%     for i = 1:1:3        %限制角度范围  在+-pi之间
%        while(argu(i)<-pi)
%           argu(i) = argu(i)+pi; 
%        end
%        while(argu(i)>pi)
%           argu(i) = argu(i)-pi; 
%        end
%     end
    Fi = eval(Fii);
    diff_argu = (-inv(eval(J))*Fi)';
    err = diff_argu(1)*diff_argu(1)+diff_argu(2)*diff_argu(2)+diff_argu(3)*diff_argu(3);
    err_2 = diff_argu(4)*diff_argu(4)+diff_argu(5)*diff_argu(5)+diff_argu(6)*diff_argu(6);
    num = num+1
end
num;
alpha = argu(1);
beta = argu(2);
gama = argu(3);
X0 = argu(4);
Y0 = argu(5);
Z0 = argu(6);
argu