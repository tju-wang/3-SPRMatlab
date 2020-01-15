clear all
clc
syms alpha beta  gama X0 Y0 Z0
syms a b ux uy uz vx vy vz wx wy wz s  %a b �ֱ�Ϊ����ƽ̨���������Բ�뾶  sΪ����ƫ��
RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
% RX_gama = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)];

%ŷ����
R = RX_alpha*RY_beta*RX_gama;
% A1_o = [0;-a;0]
% A2_o = [-(3.^(1/2)/2)*a;1/2*a;0]
% A3_o = [(3^(1/2)/2)*a;1/2*a;0]
% 
% B1 = [0;-b;0]
% B2 = [-(3^(1/2)/2)*b;1/2*b;0]
% B3 = [(3^(1/2)/2)*b;1/2*b;0]

A1_o = [0;-a;0];
A2_o = [(3.^(1/2)/2)*a;1/2*a;0];
A3_o = [-(3^(1/2)/2)*a;1/2*a;0];

% B1 = [0;-b;0];
% B2 = [(3^(1/2)/2)*b;1/2*b;0];
% B3 = [-(3^(1/2)/2)*b;1/2*b;0];
B1 = [0;-b;20];
B2 = [(3^(1/2)/2)*b;1/2*b;20];
B3 = [-(3^(1/2)/2)*b;1/2*b;20];

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3))
X0 = (b*uy*(ux-3*vy)+2*Z0*wx)/(2*wz)
Y0 = (b*ux*(vy-ux)+2*b*vx*uy+2*Z0*wy)/(2*wz)

% X0 =  (b*uy*(3*vy-ux)+2*Z0*wx)/(2*wz)    %���Ľ��  ����  
% Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*wy)/(2*wz)
%��֤������ķ���  �Ŵ����  1. alpha = 0.12  beta= 0.1  Z0=50  ����q1-3  q1=1903 q2=1941
%X0=-18 Y0=-36  ���ĵĹ�ʽ X0=28 Y0=11.4
%q3 = 1972    2.��alpha=0 beta=0  ����q1=q2=q3=1938   3.��Ӧ ����q1  q2<q3 
%���� X>0   Y<0  ���

alpha = gama

Ao = [X0;Y0;Z0];
A1 = R*A1_o + Ao;
A2 = R*A2_o + Ao;
A3 = R*A3_o + Ao;

% q1 = sqrt((B1-A1)'*(B1-A1)-s.^2);    %q��ʾ P������ + ������������ȣ�����ʼ��ƽ�У�
% q2 = sqrt((B2-A2)'*(B2-A2)-s.^2);
% q3 = sqrt((B3-A3)'*(B3-A3)-s.^2);

%��ֵ����
%������άģ�� ��������  ���D=16mm С��d=10mm
s = 62;  %������ 58.5mm + 6.5/2
alpha = 0.2;
beta = 0.1;
gama = alpha;
Z0 = 260;

a = 41.56 %��ƽ̨���Բ�뾶 
b = 80; %ʵ��ֵΪ80
q1 = eval(sqrt((B1-A1)'*(B1-A1)));    %q��ʾ P������ ��ʱ B1��λ��������� q��ֵ�� 
q2 = eval(sqrt((B2-A2)'*(B2-A2)));    %�������ͶӰ���ƶ����ϵĵ���ת����֮��ĳ���  B1��������Ҫ�޸�
q3 = eval(sqrt((B3-A3)'*(B3-A3)));

Pq = [q1 q2 q3]
X0 = eval(X0);
Y0 = eval(Y0);
A0 = [X0 Y0 Z0]


% X00 = eval(X00)
% Y00 = eval(Y00)



%����  �� X0 Y0 Z0 ��֪����µ�ת�� alpha beta gama
% f1 = X0 - (b*uy*(3*vy-ux)+2*Z0*wx)/(2*wz)
% f2 = Y0 - (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*wy)/(2*wz)

% f1 = X0 - (b*(3/2*uy*vy-uy*ux)+Z0*(uz*vy-vz*uy))/(vx*uy-ux*vy)
% 
% f2 = Y0 - (1/2*b*ux*(vy-ux)+b*uy*vx-Z0*(vz*ux-uz*vx))/(vy*ux-uy*vx)
% 
% df1_alpha = diff(f1,alpha)
% df1_beta = diff(f2,beta)
% 
% df2_alpha = diff(f2,alpha)
% df2_beta = diff(f2,beta)
% 
% X0 = (b*(3/2*uy*vy-uy*ux)+Z0*(uz*vy-vz*uy))/(vx*uy-ux*vy)
% Y0 = (1/2*b*ux*(vy-ux)+b*uy*vx-Z0*(vz*ux-uz*vx))/(vy*ux-uy*vx)
% 
% alpha = 0
% beta = 0
% err = 10
% X0 = 5
% Y0 = 5
% Z0 = 260
% a = 41.56 %��ƽ̨���Բ�뾶
% b = 80
% s = 62  %������ 58.5mm + 6.5/2
% while(err>0.00001)
%    gama = alpha
%    
%    alpha1 = (alpha-f1/df1_alpha);
%    beta1 = (beta - f1/df1_beta);
%    
%    err = (alpha-alpha1)
%    alpha= alpha1
%    beta = beta1
%     
% end
% alpha
% beta

%% ���
clear all
clc
syms alpha beta  gama X0 Y0 Z0
syms a b ux uy uz vx vy vz wx wy wz s  %a b �ֱ�Ϊ����ƽ̨���������Բ�뾶  sΪ����ƫ��
RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
% RX_gama = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
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

% X0 = (b*uy*(ux-3*vy)+2*Z0*wx)/(2*wz)
% Y0 = (b*ux*(vy-ux)+2*b*vx*uy+2*Z0*wy)/(2*wz)

% X0 =  (b*uy*(3*vy-ux)+2*Z0*wx)/(2*wz)    %���Ľ��  ����  
% Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*wy)/(2*wz)
%��֤������ķ���  �Ŵ����  1. alpha = 0.12  beta= 0.1  Z0=50  ����q1-3  q1=1903 q2=1941
%X0=-18 Y0=-36  ���ĵĹ�ʽ X0=28 Y0=11.4
%q3 = 1972    2.��alpha=0 beta=0  ����q1=q2=q3=1938   3.��Ӧ ����q1  q2<q3 
%���� X>0   Y<0  ���

alpha = gama

Ao = [X0;Y0;Z0];
A1 = R*A1_o + Ao;
A2 = R*A2_o + Ao;
A3 = R*A3_o + Ao;

% q1 = sqrt((B1-A1)'*(B1-A1)-s.^2);    %q��ʾ P������ + ������������ȣ�����ʼ��ƽ�У�
% q2 = sqrt((B2-A2)'*(B2-A2)-s.^2);
% q3 = sqrt((B3-A3)'*(B3-A3)-s.^2);

%��ֵ����
%������άģ�� ��������  ���D=16mm С��d=10mm
s = 62;  %������ 58.5mm + 6.5/2
alpha = 0;
beta = 0;
gama = alpha;
Z0 = 240;

a = 41.56 %��ƽ̨���Բ�뾶 
b = 80; %ʵ��ֵΪ80
q1 = eval(sqrt((B1-A1)'*(B1-A1)));    %q��ʾ P������ ��ʱ B1��λ��������� q��ֵ�� 
q2 = eval(sqrt((B2-A2)'*(B2-A2)));    %�������ͶӰ���ƶ����ϵĵ���ת����֮��ĳ���  B1��������Ҫ�޸�
q3 = eval(sqrt((B3-A3)'*(B3-A3)));

Pq = [q1 q2 q3]
X0 = eval(X0);
Y0 = eval(Y0);
A0 = [X0 Y0 Z0]






