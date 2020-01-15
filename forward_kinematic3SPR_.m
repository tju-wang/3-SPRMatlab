%% 将3-SPR拓展为Stewart机构 在R转轴处  将一个点拓展为两个点  构成六条支链 求解得到六个变量
clear all
clc
syms alpha beta  gama X0 Y0 Z0 real
syms a b ux uy uz vx vy vz wx wy wz s real %a b 分别为静动平台三角形外接圆半径  s为机构偏距
RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
% RX_gama = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)];
RZ_gama = [cos(gama),-sin(gama),0;sin(gama),cos(gama),0;0,0,1];

%欧拉角
%R = RX_alpha*RY_beta*RZ_gama;
% R = RZ_gama*RY_beta*RX_alpha;
R = RX_gama*RY_beta*RX_alpha;

% A2_o = [-(3.^(1/2)/2)*a;1/2*a;0]
% A3_o = [(3^(1/2)/2)*a;1/2*a;0]
% A1_o = [0;-a;0]
% 
% B2 = [-(3^(1/2)/2)*b;1/2*b;0]
% B3 = [(3^(1/2)/2)*b;1/2*b;0]
% B1 = [0;-b;0]
A1_o = [0;-a;0];
A2_o = [(3.^(1/2)/2)*a;1/2*a;0];
A3_o = [-(3^(1/2)/2)*a;1/2*a;0];

B1 = [0;-b;20];
B2 = [(3^(1/2)/2)*b;1/2*b;20];
B3 = [-(3^(1/2)/2)*b;1/2*b;20];

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
Fii = [ BA(1,:)*BA(1,:)'-q11^2;BA(2,:)*BA(2,:)'-q12^2; BA(3,:)*BA(3,:)'-q21^2;BA(4,:)*BA(4,:)'-q22^2; BA(5,:)*BA(5,:)'-q31^2;BA(6,:)*BA(6,:)'-q32^2 ];

argu = [alpha beta gama X0 Y0 Z0];
J = []
J = jacobian(Fi,argu);

% 开始数值运算  结构常数等进行赋值
s = 62  %测量得 58.5mm + 6.5/2
alpha = 0.5;
beta = 0.6;
gama = 0.8;
Z0 = 100;
X0 = 10;
Y0 = 10;
q1 =271;   q2 =268;   q3 =261;
q1 =281;   q2 =271;   q3 = 261;    
a = 41.56; %动平台外接圆半径
b = 80;
argu = eval(argu);

err = 1
Fi = eval(Fi);


%%  迭代
% Fi = eval([-q1^2+(B1-A1)'*(B1-A1)-s.^2;-q2^2+(B2-A2)'*(B2-A2)-s.^2;-q3^2+(B3-A3)'*(B3-A3)-s.^2])
%  diff_argu' = - inv(eval(J))*Fi;
diff_argu = [0 0 0 0 0 0]
num = 0
while(err>1.0e-15)
   argu = (argu + diff_argu);
   for i = 1:1:3        %限制角度范围  在+-pi之间
       while(argu(i)<-pi)
          argu(i) = argu(i)+pi; 
       end
       while(argu(i)>pi)
          argu(i) = argu(i)-pi; 
       end

   end
   
   alpha = argu(1);
   beta = argu(2);
   gama = argu(3);
   X0 = argu(4);
   Y0 = argu(5);
   Z0 = argu(6);
%    Z0 = argu(3);
%    gama = alpha;
%    Fi = eval([-q1^2+(B1-A1)'*(B1-A1)-s.^2;-q2^2+(B2-A2)'*(B2-A2)-s.^2;-q3^2+(B3-A3)'*(B3-A3)-s.^2]);
    Fi = eval(Fii);
    diff_argu = (-inv(eval(J))*Fi)';
   err = diff_argu*diff_argu';
   num = num+1
end
num
argu
err
% X0 = eval(X0)
% Y0 = eval(Y0)


% Q = [q1,q2,q3];
% q1 =271.0712;   q2 = 256;   q3 =260;
% while(1)
%     
%     
%     q2 = q2+1;
%     err = 1;
%     % Fi = eval([-q1^2+(B1-A1)'*(B1-A1)-s.^2;-q2^2+(B2-A2)'*(B2-A2)-s.^2;-q3^2+(B3-A3)'*(B3-A3)-s.^2])
%     %  diff_argu' = - inv(eval(J))*Fi;
%     diff_argu = [0 0 0 0 0 0];
%     num = 0;
%     while(err>1.0e-25)
%        argu = (argu + diff_argu);
%        alpha = argu(1);
%        beta = argu(2);
%        gama = argu(3);
%        X0 = argu(4);
%        Y0 = argu(5);
%        Z0 = argu(6);
% %    Z0 = argu(3);
% %    gama = alpha;
% %    Fi = eval([-q1^2+(B1-A1)'*(B1-A1)-s.^2;-q2^2+(B2-A2)'*(B2-A2)-s.^2;-q3^2+(B3-A3)'*(B3-A3)-s.^2]);
%        Fi = eval(Fii);
%        diff_argu = (-inv(eval(J))*Fi)'
%        err = diff_argu*diff_argu';
%        num = num+1;
%     end
%     num
%     argu
%     Q = [q1,q2,q3]
%     
% end













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
% A2_o = [-(3.^(1/2)/2)*a;1/2*a;0]
% A3_o = [(3^(1/2)/2)*a;1/2*a;0]
% A1_o = [0;-a;0]
% 
% B2 = [-(3^(1/2)/2)*b;1/2*b;0]
% B3 = [(3^(1/2)/2)*b;1/2*b;0]
% B1 = [0;-b;0]
A1_o = [0;-a;0]
A2_o = [(3.^(1/2)/2)*a;1/2*a;0]
A3_o = [-(3^(1/2)/2)*a;1/2*a;0]

B1 = [0;-b;0]
B2 = [(3^(1/2)/2)*b;1/2*b;0]
B3 = [-(3^(1/2)/2)*b;1/2*b;0]

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3))
% X0 = (b*uy*(ux-3*vy)+2*Z0*wx)/(2*wz)
% Y0 = (b*ux*(vy-ux)+2*b*vx*uy+2*Z0*wy)/(2*wz)
X0 =  (b*uy*(3*vy-ux)+2*Z0*wx)/(2*wz)    %论文结果 
Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*wy)/(2*wz)
% X0 = 2*Z0*wx/(2*wz)
% Y0 = 2*Z0*wy/(2*wz)

alpha = gama

Ao = [X0;Y0;Z0]
A1 = R*A1_o + Ao
A2 = R*A2_o + Ao
A3 = R*A3_o + Ao


f = [(B1-A1)'*(B1-A1)-s.^2;(B2-A2)'*(B2-A2)-s.^2;(B3-A3)'*(B3-A3)-s.^2]
argu = [alpha beta Z0]
% JJ=[diff(f(1),alpha),diff(f(1),beta),diff(f(1),Z0);
%     diff(f(2),alpha),diff(f(2),beta),diff(f(2),Z0);
%     diff(f(3),alpha),diff(f(3),beta),diff(f(3),Z0)];
J = []
J = jacobian(f,argu)

% J_1 = inv(J)
% 开始数值运算  结构常数等进行赋值
s = 62  %测量得 58.5mm + 6.5/2
alpha = 0.2;
beta = 0.3;
gama = alpha;
Z0 = 250;
q1 =287;   q2 = 309;   q3 = 322;

a = 41.56 %动平台外接圆半径
b = 80
argu = eval(argu)

err = 1
Fi = eval([-q1^2+(B1-A1)'*(B1-A1)-s.^2;-q2^2+(B2-A2)'*(B2-A2)-s.^2;-q3^2+(B3-A3)'*(B3-A3)-s.^2])
%  diff_argu' = - inv(eval(J))*Fi;
diff_argu = [0 0 0]
num = 0
while(err>0.000001)
   argu = (argu + diff_argu);
   alpha = argu(1);
   beta = argu(2);
   Z0 = argu(3);
   gama = alpha;
   Fi = eval([-q1^2+(B1-A1)'*(B1-A1)-s.^2;-q2^2+(B2-A2)'*(B2-A2)-s.^2;-q3^2+(B3-A3)'*(B3-A3)-s.^2]);
   diff_argu = (-inv(eval(J))*Fi)';
   err = diff_argu(1)*diff_argu(1)+diff_argu(2)*diff_argu(2);
   num = num+1
end
num
argu
err
X0 = eval(X0)
Y0 = eval(Y0)
