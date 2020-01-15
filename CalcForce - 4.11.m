%��������������� �����
%����  q1-q3   ���Fm1 - Fm3
%����1.���3-SPR��������  ��ĩ��X Y Z�����  2.��������������(������Ч) 
%3.��ƽ�ⷽ�� ���P������  4.�Ƚ�������Դ�С  ��ϵ�� ȷ����������������PWM������

%�����motor_force.m�ļ�  ������B������  ������B������ 
%% 1.
%��������
clear all
clc
syms alpha beta  gama X0 Y0 Z0 real
syms a b ux uy uz vx vy vz wx wy wz s real%a b �ֱ�Ϊ����ƽ̨���������Բ�뾶  sΪ����ƫ��
%C������ P������q
syms Xc1 Xc2 Xc3 Yc1 Yc2 Yc3 Zc1 Zc2 Zc3 q1 q2 q3 real
RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
% RX_gama = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
% RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)];
RZ_gama = [cos(gama),-sin(gama),0;sin(gama),cos(gama),0;0,0,1];
%ŷ����
R = RX_alpha*RY_beta*RZ_gama;
A1_o = [0;-a;0]
A2_o = [(3.^(1/2)/2)*a;1/2*a;0]
A3_o = [-(3^(1/2)/2)*a;1/2*a;0]

B1 = [0;-b;20]
B2 = [(3^(1/2)/2)*b;1/2*b;20]
B3 = [-(3^(1/2)/2)*b;1/2*b;20]

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3));
% X0 = (b*uy*(ux-3*vy)+2*Z0*wx)/(2*wz);
% Y0 = (b*ux*(vy-ux)+2*b*vx*uy+2*Z0*wy)/(2*wz);
% X0 =  (b*uy*(3*vy-ux)+2*Z0*wx)/(2*wz);    %���Ľ�� 
% Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*wy)/(2*wz);
% X0 = 2*Z0*wx/(2*wz)
% Y0 = 2*Z0*wy/(2*wz)

Ao = [X0;Y0;Z0];
A1 = R*A1_o + Ao;
A2 = R*A2_o + Ao;
A3 = R*A3_o + Ao;
%�����õ���6��֧��
r1 = A2' - A3';
r1 = r1/norm(r1);
r2 = A1' - A3';
r2 = r2/norm(r2);
r3 = A1' - A2';
r3 = r3/norm(r3);
L = 40;     %�����R�����ȵ�һ����20mm
A11 = A1' + r1*L;
A12 = A1' - r1*L;
A21 = A2' + r2*L;
A22 = A2' - r2*L;
A31 = A3' + r3*L;
A32 = A3' - r3*L;

%���6��֧��
syms q11 q12 q21 q22 q31 q32 q1 q2 q3 real

q11 = sqrt(q1^2 + L^2);   q12 = q11;
q21 = sqrt(q2^2 + L^2);   q22 = q21;
q31 = sqrt(q3^2 + L^2);   q32 = q31;

BA = [ B1'-A11;B1'-A12; B2'-A21;B2'-A22; B3'-A31;B3'-A32];
Fi = [ BA(1,:)*BA(1,:)'-q11^2;BA(2,:)*BA(2,:)'-q12^2; BA(3,:)*BA(3,:)'-q21^2;BA(4,:)*BA(4,:)'-q22^2; BA(5,:)*BA(5,:)'-q31^2;BA(6,:)*BA(6,:)'-q32^2 ];
Fii =[ BA(1,:)*BA(1,:)'-q11^2;BA(2,:)*BA(2,:)'-q12^2; BA(3,:)*BA(3,:)'-q21^2;BA(4,:)*BA(4,:)'-q22^2; BA(5,:)*BA(5,:)'-q31^2;BA(6,:)*BA(6,:)'-q32^2 ];

argu = [alpha beta gama X0 Y0 Z0];
J = []
J = jacobian(Fi,argu);

%% 2.��C������ ������ı�� ��������
% 3.��ƽ�ⷽ�� ���P������   ����ƫ����P������C1 C2 C3
%��C�������Ϊ VBC vAC��ֱ vACģ����֪ vBC��R��������vAxi ��ֱ
C1 = [Xc1,Yc1,Zc1]
C2 = [Xc2,Yc2,Zc2]
C3 = [Xc3,Yc3,Zc3]

vA1C1 = vector(A1',C1) 
vA2C2 = vector(A2',C2)
vA3C3 = vector(A3',C3)

vB1C1 = vector(B1,C1) 
vB2C2 = vector(B2,C2)
vB3C3 = vector(B3,C3)

vAxi_A1 = A2' - A3';
vAxi_A2 = A1' - A3';
vAxi_A3 = A1' - A2';
%% 3.��ֵ���� ��������
% ��ʼ��ֵ����  �ṹ�����Ƚ��и�ֵ
s = 62  %������ 58.5mm + 6.5/2
alpha = 0.1;
beta = 0.15;
gama = 0.12;
Z0 = 280;
X0 = 10;
Y0 = 10;
q1 =261;   q2 =261;   q3 =261;
    
a = 41.56; %��ƽ̨���Բ�뾶
b = 80;
argu = eval(argu);

err = 1
Fi = eval(Fi);
% Fi = eval([-q1^2+(B1-A1)'*(B1-A1)-s.^2;-q2^2+(B2-A2)'*(B2-A2)-s.^2;-q3^2+(B3-A3)'*(B3-A3)-s.^2])
%  diff_argu' = - inv(eval(J))*Fi;
diff_argu = [0 0 0 0 0 0]
num = 0
while(err>1.0e-15)
    argu = (argu + diff_argu);
    alpha = argu(1);
    beta = argu(2);
    gama = argu(3);
    X0 = argu(4);
    Y0 = argu(5);
    Z0 = argu(6);
    for i = 1:1:3        %���ƽǶȷ�Χ  ��+-pi֮��
       while(argu(i)<-pi)
          argu(i) = argu(i)+pi; 
       end
       while(argu(i)>pi)
          argu(i) = argu(i)-pi; 
       end
    end
    Fi = eval(Fii);
    diff_argu = (-inv(eval(J))*Fi)'
    err = diff_argu*diff_argu';
    num = num+1
end
num
argu = (argu + diff_argu);
alpha = argu(1);
beta = argu(2);
gama = argu(3);
X0 = argu(4);
Y0 = argu(5);
Z0 = argu(6);
argu
%% 4.��ֵ���� ��C������
%��C������
% fc1_1 = eval(norm(vB1C1)-s);
% fc1_2 = eval(vB1C1 * vAxi_A1);
% fc1_3 = eval(norm(vA1C1)-q1);
% [Xc1,Yc1,Zc1] = solve(fc1_1,fc1_2,fc1_3,Xc1,Yc1,Zc1)

% [Xc2,Yc2,Zc2] = solve(fc2_1,fc2_2,fc2_3,Xc2,Yc2,Zc2)

% fc3_1 = eval(norm(vB3C3)-s);
% fc3_2 = eval(vB3C3 * vAxi_A3);
% fc3_3 = eval(norm(vA3C3)-q3);
% [Xc3,Yc3,Zc3] = solve(fc3_1,fc3_2,fc3_3,Xc3,Yc3,Zc3)

%���C���� ��ֵ����C������ʽ
syms Xc Yc Zc XB YB ZB XA YA ZA Xi Yi Zi real
Xc = Xc1; Yc = Yc1; Zc = Zc1; XB = B1(1); YB = B1(2); ZB = B1(3); XA = A1(1); YA = A1(2); ZA = A1(3);
Xi = A2(1)-A3(1);  Yi = A2(2)-A3(2);  

 Zi = A2(3)-A3(3);
M1 = Xi*XB+Yi*YB+Zi*ZB;
K1 = (2*(XA-XB)-2*(ZA-ZB)*Xi/Zi);
K2 = 2*(YA-YB)-2*(ZA-ZB)*Yi/Zi;
M2 = 2*(ZA-ZB)*M1/Zi+(XB.^2+YB.^2+ZB.^2)-(XA.^2+YA.^2+ZA.^2)+q1.^2-2*s.^2;
K3 = -K2/K1;
M3 = -M2/K1;
K4 = -K1/K2;
M4 = -M2/K2;
K5 = -(K3*Xi+Yi)/Zi;
M5 = (M1-Xi*M3)/Zi;
K6 = -(Xi+K4*Yi)/Zi;
M6 = (M1-Yi*M4)/Zi;
%��һԪ���η�����
ca = eval(K3.^2+K5.^2+1);
cb = eval(2*K3*M3+2*K5*M5-2*YB-2*XB*K3-2*ZB*K5);
cc = eval(M3.^2+M5.^2-2*XB*M3-2*ZB*M5+XB.^2+YB.^2+ZB.^2-s.^2);
if real((((cb.^2-4*ca*cc)))) > 0
    Yc = real(((-cb-sqrt(cb.^2-4*ca*cc))/2/ca));
    Xc = real(eval((-K2*Yc-M2)/K1));
    Zc = real(eval((M1-(Xi*Xc+Yi*Yc))/Zi));
    Xc1 = Xc
    Yc1 = Yc
    Zc1 = Zc
% K3 = eval(K3)
% K5 = eval(K5)
% K4 = eval(K4)
% K6 = eval(K6)
else
    sprintf('�����㷽���������')
end

fc1_1 = eval(norm(vB1C1)-s);
fc1_2 = eval(vB1C1 * vAxi_A1');
fc1_3 = eval(norm(vA1C1)-sqrt(q1^2-s^2));
fc1 = [fc1_1 fc1_2 fc1_3]
 
%������֤  ��C2 C3
%��һԪ���η����� 

XB = B2(1); YB = B2(2); ZB = B2(3); XA = A2(1); YA = A2(2); ZA = A2(3);
Xi = A1(1)-A3(1);  Yi = A1(2)-A3(2); Zi;
Zi = A1(3)-A3(3);

M1 = Xi*XB+Yi*YB+Zi*ZB;
K1 = (2*(XA-XB)-2*(ZA-ZB)*Xi/Zi);
K2 = 2*(YA-YB)-2*(ZA-ZB)*Yi/Zi;
M2 = 2*(ZA-ZB)*M1/Zi+(XB.^2+YB.^2+ZB.^2)-(XA.^2+YA.^2+ZA.^2)+q2.^2-2*s.^2;
K3 = -K2/K1;
M3 = -M2/K1;
K4 = -K1/K2;
M4 = -M2/K2;
K5 = -(K3*Xi+Yi)/Zi;
M5 = (M1-Xi*M3)/Zi;
K6 = -(Xi+K4*Yi)/Zi;
M6 = (M1-Yi*M4)/Zi;

ca = eval(K4.^2+K6.^2+1);
cb = eval(2*K4*M4+2*K6*M6-2*XB-2*K4*YB-2*ZB*K6);
cc = eval(M4.^2+M6.^2-2*YB*M4-2*ZB*M6+XB.^2+YB.^2+ZB.^2-s.^2);
% K3 = eval(K3)
% K5 = eval(K5)
% K4 = eval(K4)
% K6 = eval(K6)
if real((((cb.^2-4*ca*cc)))) > 0
    Xc = real(((-cb+sqrt(cb.^2-4*ca*cc))/2/ca));
    Yc = real(eval(K4*Xc+M4));
    Zc = real(eval(K6*Xc+M6));
    Xc2 = Xc
    Yc2 = Yc
    Zc2 = Zc
else
    sprintf('�����㷽���������')
end
fc2_1 = eval(norm(vB2C2)-s);
fc2_2 = eval(vB2C2 * vAxi_A2');
fc2_3 = eval(norm(vA2C2)-sqrt(q2^2-s^2));
fc2 = [fc2_1 fc2_2 fc2_3]

XB = B3(1); YB = B3(2); ZB = B3(3); XA = A3(1); YA = A3(2); ZA = A3(3);
Xi = A1(1)-A2(1);  Yi = A1(2)-A2(2); 
Zi = A1(3)-A2(3);
M1 = Xi*XB+Yi*YB+Zi*ZB;
K1 = (2*(XA-XB)-2*(ZA-ZB)*Xi/Zi);
K2 = 2*(YA-YB)-2*(ZA-ZB)*Yi/Zi;
M2 = 2*(ZA-ZB)*M1/Zi+(XB.^2+YB.^2+ZB.^2)-(XA.^2+YA.^2+ZA.^2)+q3.^2-2*s.^2;
K3 = -K2/K1;
M3 = -M2/K1;
K4 = -K1/K2;
M4 = -M2/K2;
K5 = -(K3*Xi+Yi)/Zi;
M5 = (M1-Xi*M3)/Zi;
K6 = -(Xi+K4*Yi)/Zi;
M6 = (M1-Yi*M4)/Zi;

ca = eval(K4.^2+K6.^2+1);
cb = eval(2*K4*M4+2*K6*M6-2*XB-2*K4*YB-2*ZB*K6);
cc = eval(M4.^2+M6.^2-2*YB*M4-2*ZB*M6+XB.^2+YB.^2+ZB.^2-s.^2);
% K3 = eval(K3)
% K5 = eval(K5)
% K4 = eval(K4)
% K6 = eval(K6)
if real((((cb.^2-4*ca*cc)))) > 0
    Xc = real(((-cb-sqrt(cb.^2-4*ca*cc))/2/ca));
    Yc = real(eval(K4*Xc+M4));
    Zc = real(eval(K6*Xc+M6));
    Xc3 = Xc
    Yc3 = Yc
    Zc3 = Zc
else
    sprintf('�����㷽���������')
end
fc3_1 = eval(norm(vB3C3)-s);
fc3_2 = eval(vB3C3 * vAxi_A3');
fc3_3 = eval(norm(vA3C3)-sqrt(q3^2-s^2));
fc3 = [fc3_1 fc3_2 fc3_3]

%% 5.���ʵ������
%��ʾ�����ʵ� �����������  �������ӵ���  ��10���ʵ�  C���� ��Ҫ�Ѽ������Ĳ��� �ŵ�һ��
%���Ĵ�СΪ����  �������õ�Ϊ��������
%10���������  ���� C������ ����Ϊ��������ӿ鲿������ ��Ҫ���� 

syms FGy Gx Gy Gz real
syms Fx_b1 Fy_b1 Fz_b1 Fx_b2 Fy_b2 Fz_b2 Fx_b3 Fy_b3 Fz_b3 real    %������� ��������µ���
syms TG1_1 TG1_2 TG1_3 TG2_1 TG2_2 TG2_3 TG3_1 TG3_2 TG3_3 TGm real 
%TG1_1 Ť��-��һ��֧�� ���ؿ� _2 ����  _3 ��� ����

%����������Ĵ�С  ����
FB1 = [Fx_b1,Fy_b1,Fz_b1];
FB2 = [Fx_b2,Fy_b2,Fz_b2];
FB3 = [Fx_b3,Fy_b3,Fz_b3];
FG = [0,FGy,0];
G = [Gx,Gy,Gz];
%  ������
B0 = [0;0;20]
TB1 = cross(B1-B0,FB1);
TB2 = cross(B2-B0,FB2);
TB3 = cross(B3-B0,FB3);
%��Fz ����x������ط��� TG(1) Fz_b3*b Fz_b2*b/2 Fz_b3*b/2 
%��ΪB1 B2 B3��ԭ����ͬһƽ���� ����y����û�ж�����x��Ť��
%���峣��
FG_1 = [0,438*10,0];%��λg   ���ĵ�λ mN
FG_2 = [0,285.75*10,0];
FG_3 = [0,342.5*10,0];
FGm = [0,(86+30)*10,0];
Q = 310;%�����ܳ��� 310mm
P1 = clc_P(eval(vA1C1),eval(C1),Q-sqrt(q1^2-s^2));  %PΪ���ؿ�����
P2 = clc_P(eval(vA2C2),eval(C2),Q-sqrt(q2^2-s^2));
P3 = clc_P(eval(vA3C3),eval(C3),Q-sqrt(q3^2-s^2));

D1 = eval((P1+A1')/2);      %DΪ��������
D2 = eval((P2+A2')/2);
D3 = eval((P3+A3')/2);
Gm = [X0 Y0 Z0];
%�������������  Ť��TG
TG1_1 = cross(P1-B0',FG_1);  %��B0��ȡ��
TG2_1 = cross(P2-B0',FG_1);
TG3_1 = cross(P3-B0',FG_1);

TG1_2 = cross(D1-B0',FG_2);
TG2_2 = cross(D2-B0',FG_2);
TG3_2 = cross(D3-B0',FG_2);

TG1_3 = cross(C1-B0',FG_3);
TG2_3 = cross(C2-B0',FG_3);
TG3_3 = cross(C3-B0',FG_3);
TGm = cross(Gm-B0',FGm);

TG = eval(TG1_1+TG2_1+TG3_1  +TG1_2+TG2_2+TG3_2  +TG1_3+TG2_3+TG3_3  +TGm);
FG = FG_1*3+FG_2*3+FG_3*3+FGm;
%Ҫ������ƽ��   ����ƽ��
Fz_b2 = solve(2*b*Fz_b2 + b/2*Fz_b2 + b/2*Fz_b2 + TG(1),Fz_b2);
Fz_b3 = Fz_b2;
Fz_b1 = -2*Fz_b2;
%��Fy ������������ƽ�ⷽ��
Fy_b2 = -(FG(2))/3;
Fy_b3 = Fy_b2;
Fy_b1 = Fy_b2;
%��Fx ����z���ת�ط��� TG(3)+Fx_b1*3*b=0
Fx_b2 = -TG(3)/3/b;
Fx_b1 = -2*Fx_b2; %Fx_b2��Fx_b1�����෴
Fx_b3 = Fx_b2;
%F �ĵ�λΪ mN
FB1 = vpa(eval(FB1),8)
FB2 = vpa(eval(FB2),8)
FB3 = vpa(eval(FB3),8)

%% ���P�� �����ؽ���Ҫ��������  ����ƽ�ⷽ��
% ԭʼ˼· P������ �� ������������ �����������ƽ��  ���ԣ���FB1 FB2 FB3�ڵ��Ӻ�ֻʣ������
% n1 = eval((vA1C1)/norm((vA1C1)))
% n2 = eval(vA2C2/norm(vA2C2))
% n3 = eval(vA3C3/norm(vA3C3))
% syms f1 f2 f3
% N = ([n1;n2;n3])
% b = (-eval(FB1+FB2+FB3+FG))'
% X = N\b

%�ڵ���P��������ƽ��  ������� ���ӿ顢������� P����������������ƽ��{�����õ㲻ͬ ��֪�ܷ�ֱ�ӵ���}
n1 = eval((vA1C1)/norm((vA1C1)));  %AC���� P������
if n1(3)>0
    n1 = -n1;
end
vn2 = A2'-A3';   %ƽ����R��ת�᷽��
n2 = eval(vn2/norm(vn2));
vn3 = cross(vA1C1,vn2);    %��ֱ��  AC��R��ת�������ɵ�ƽ��
n3 = eval(vn3/norm(vn3));

N = ([n1;n2;n3])'
bn = (-eval(FB1+FG_3))';
F1 = N\bn

clear n1 n2 n3 vn2 vn3 N bn
n1 = eval((vA2C2)/norm((vA2C2)));  %AC���� P������
if n1(3)>0
    n1 = -n1;
end
vn2 = A1'-A3';   %ƽ����R��ת�᷽��
n2 = eval(vn2/norm(vn2));
vn3 = cross(vA2C2,vn2);    %��ֱ��  AC��R��ת�������ɵ�ƽ��
n3 = eval(vn3/norm(vn3));
if n2(2)>0
    n2 = -n2;
end
if n3(1)>0
    n3 = -n3;
end

N = ([n1;n2;n3])';
bn = (-eval(FB2+FG_3))';
F2 = N\bn

clear n1 n2 n3 vn2 vn3 N bn
n1 = eval((vA3C3)/norm((vA3C3)));  %AC���� P������
if n1(3)>0
    n1 = -n1;
end
vn2 = A1'-A2';   %ƽ����R��ת�᷽��
n2 = eval(vn2/norm(vn2));
vn3 = cross(vA3C3,vn2);    %��ֱ��  AC��R��ת�������ɵ�ƽ��
n3 = eval(vn3/norm(vn3));

if n2(2)>0
    n2 = -n2;
end
if n3(1)<0
    n3 = -n3;
end
N = ([n1;n2;n3])';
bn = (-eval(FB3+FG_3))';
F3 = N\bn

 