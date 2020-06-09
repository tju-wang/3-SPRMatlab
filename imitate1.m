clear all
clc
tic  %��ʱ��ʼ
syms alpha beta  gama X0 Y0 Z0
syms a b ux uy uz vx vy vz wx wy wz s  %a b �ֱ�Ϊ����ƽ̨���������Բ�뾶  sΪ����ƫ��
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

Ao = [X0;Y0;Z0];
A1 = R*A1_o + Ao;
A2 = R*A2_o + Ao;
A3 = R*A3_o + Ao;

%��ֵ����
%������άģ�� ��������  ���D=16mm С��d=10mm
s = 62;  %������ 58.5mm + 6.5/2
alpha = 0.12;
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
numm = 1
result = []     %��������Ľ��

for i=1:10
    alpha = -0.2 + 0.04*i;
    for j=1:5
        beta = -0.1 + 0.04*j;    
        gama = alpha;
        Z0 = 260;
        X0 = (b*uy*(3*vy-ux)+2*Z0*(vz*uy-vy*uz))/(2*(ux*vy-vx*uy));   %���Ľ�� 
        Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*(uz*vx-vz*ux))/(2*(ux*vy-vx*uy));
        alpha = gama;

        q1 = eval(sqrt((B1-A1)'*(B1-A1)));    %q��ʾ P������ ��ʱ B1��λ��������� q��ֵ�� 
        q2 = eval(sqrt((B2-A2)'*(B2-A2)));    %�������ͶӰ���ƶ����ϵĵ���ת����֮��ĳ���  B1��������Ҫ�޸�
        q3 = eval(sqrt((B3-A3)'*(B3-A3)));

        Pq = [q1 q2 q3];
        X0 = eval(X0);
        Y0 = eval(Y0);
        A0 = [X0 Y0 Z0];
        result(numm,:) = [numm Pq alpha beta gama A0];
        numm = numm+1
    end
    
end
time1 = toc
result

% figure(1)
% plot(result(:,5),result(:,6))
% figure(2)
% plot(result(:,1),result(:,2))

%% ʹ���������  ����

clear alpha beta gama X0 Y0 Z0 q1 q2 q3
syms alpha beta  gama X0 Y0 Z0 q1 q2 q3
%�޸�����ϵ  x����������  ��   y����������   z �������ɾ�ƽָ̨��ƽ̨
X0 =  (b*uy*(3*vy-ux)+2*Z0*(vz*uy-vy*uz))/(2*(ux*vy-vx*uy))   %���Ľ�� 
Y0 = (b*ux*(ux-vy)-2*b*vx*uy+2*Z0*(uz*vx-vz*ux))/(2*(ux*vy-vx*uy))

alpha = gama

Ao = [X0;Y0;Z0]
A1 = R*A1_o + Ao
A2 = R*A2_o + Ao
A3 = R*A3_o + Ao


f = [(B1-A1)'*(B1-A1);(B2-A2)'*(B2-A2);(B3-A3)'*(B3-A3)];
argu = [alpha beta Z0];
J = []
J = jacobian(f,argu);

% J_1 = inv(J)
% ��ʼ��ֵ����  �ṹ�����Ƚ��и�ֵ
s = 62  %������ 58.5mm + 6.5/2
alpha = -0.1;
beta = 0.1;
gama = alpha;
Z0 = 200;
argu = eval(argu)
numm=1
kineResult3 = []
err = [];  err_2 = [];
tic
for i=1:10
   
    for j=1:5
        q1 = result(numm,2); q2 = result(numm,3);  q3 = result(numm,4);
        Fi = eval([-q1^2+(B1-A1)'*(B1-A1);-q2^2+(B2-A2)'*(B2-A2);-q3^2+(B3-A3)'*(B3-A3)]);
        diff_argu = [0 0 0];
        num = 1;
        err(numm,num) = 1;  err_2(numm,num) = 1;
        while(err(numm,num)>1.0e-4 || err_2(numm,num)>1.0e-4 )
           argu = (argu + diff_argu);
           alpha = argu(1);
           beta = argu(2);
           Z0 = argu(3);
           gama = alpha;
           Fi = eval([-q1^2+(B1-A1)'*(B1-A1);-q2^2+(B2-A2)'*(B2-A2);-q3^2+(B3-A3)'*(B3-A3)]);
           diff_argu = (-inv(eval(J))*Fi)';
          
           num = num+1
           err(numm,num) = diff_argu(1)*diff_argu(1)+diff_argu(2)*diff_argu(2);
           err_2(numm,num) = (diff_argu(3)*diff_argu(3));
           
        end
        kineResult3(numm,:) = [numm q1 q2 q3 argu(1) argu(2) argu(1) eval(X0) eval(Y0) Z0 num];
        numm = numm+1
    end
end
time2 = toc
kineResult3

%% 6��  ����
clear alpha beta  gama X0 Y0 Z0 a b ux uy uz vx vy vz wx wy wz s Xc1 Xc2 Xc3 Yc1 Yc2 Yc3 Zc1 Zc2 Zc3 q1 q2 q3
clc
syms alpha beta  gama X0 Y0 Z0 real
syms a b ux uy uz vx vy vz wx wy wz s real%a b �ֱ�Ϊ����ƽ̨���������Բ�뾶  sΪ����ƫ��
%C������ P������q
syms Xc1 Xc2 Xc3 Yc1 Yc2 Yc3 Zc1 Zc2 Zc3 q1 q2 q3 real

RX_alpha = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
RY_beta = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
% RX_gama = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)]
RX_gama = [1,0,0;0,cos(gama),-sin(gama);0,sin(gama),cos(gama)];

RZ_alpha = [cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
RX_beta = [1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
RZ_gama = [cos(gama),-sin(gama),0;sin(gama),cos(gama),0;0,0,1];
%ŷ����
%R = RZ_alpha*RX_beta*RZ_gama  %����ڶ�����ϵ  XYZŷ����
R = RX_alpha*RY_beta*RZ_gama;

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3));
%R  =[ux,vx,wx;uy,vy,wy;uz,vz,wz]

%�޸�����ϵ  x����������  ��   y����������   z �������ɾ�ƽָ̨��ƽ̨
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

%����6����
r1 = A2' - A3';
r1 = r1/norm(r1);
r2 = A1' - A3';
r2 = r2/norm(r2);
r3 = A1' - A2';
r3 = r3/norm(r3);
L = 100;     %�����R�����ȵ�һ����20mm
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
J = [];
J = jacobian(Fi,argu);

% 3.��ֵ���� ��������
% ��ʼ��ֵ����  �ṹ�����Ƚ��и�ֵ
s = 62;  %������ 58.5mm + 6.5/2
a = 41.56; %��ƽ̨���Բ�뾶
b = 80;
alpha = 0.3;
beta = -0.3;
gama = 0.3;
Z0 = 250;
X0 = 10;
Y0 = 10;
argu = eval(argu);
numm=1
kineResult6 = []
  err = [];err_2 = [];
tic
for i=1:10
   for j=1:5
    q1 = result(numm,2); q2 = result(numm,3);  q3 = result(numm,4);
    Fi = eval(Fii);
    diff_argu = [0 0 0 0 0 0];
    num = 1;
    err(numm,num)=1; err_2(numm,num)=1;
    while(err(numm,num)>1.0e-4 || err_2(numm,num)>1.0e-4)
        argu = (argu + diff_argu);
        alpha = argu(1);
        beta = argu(2);
        gama = argu(3);
        X0 = argu(4);
        Y0 = argu(5);
        Z0 = argu(6);
        Fi = eval(Fii);
        diff_argu = (-inv(eval(J))*Fi)';
        num = num+1
        err(numm,num) = diff_argu(1)*diff_argu(1)+diff_argu(2)*diff_argu(2);
        err_2(numm,num) = diff_argu(6)*diff_argu(6);
        
    end
    kineResult6(numm,:) = [numm q1 q2 q3 argu num];
    numm = numm+1
   end
end
time3 = toc
kineResult6

%% ��ͼ
% %x  y �Ĺ켣��״
% figure(9)
% plot(result(:,8),result(:,9),'r -.')

figure(1)
plot(result(:,1),log(kineResult6(:,13)),'r -.')
hold on
plot(result(:,1),log(kineResult3(:,13)),'g -.')
title('���������Ա�')
hold off

figure(2)
plot(result(:,1),(result(:,8)),'b -.')
hold on
plot(result(:,1),(kineResult6(:,8)),'r *')
hold on
plot(result(:,1),(kineResult3(:,8)),'g s')
title('x���������֤')
hold off

figure(3)
plot(result(:,1),real(result(:,9)),'b -.')
hold on
plot(result(:,1),real(kineResult6(:,9)),'r *')
hold on
plot(result(:,1),real(kineResult3(:,9)),'g s')
title('y���������֤')
hold off

figure(4)
plot(result(:,1),real(result(:,5)),'b -.')
hold on
plot(result(:,1),real(kineResult6(:,5)),'r *')
hold on
plot(result(:,1),real(kineResult3(:,5)),'g s')
title('alpha�������֤')
hold off

figure(5)
plot(result(:,1),real(result(:,6)),'b -.')
hold on
plot(result(:,1),real(kineResult6(:,6)),'r *')
hold on
plot(result(:,1),real(kineResult3(:,6)),'g s')
title('beta�������֤')
hold off

figure(6)
plot(result(:,1),power((result(:,8)-kineResult6(:,8)),2),'r -.')
hold on
plot(result(:,1),power((result(:,8)-kineResult3(:,8)),2),'g -.')
title('x�ᾫ��')
hold off

figure(7)
plot(result(:,1),power((result(:,9)-kineResult6(:,9)),2),'r -.')
hold on
plot(result(:,1),power((result(:,9)-kineResult3(:,9)),2),'g -.')
title('y�ᾫ��')
hold off

figure(8)
var = 10
plot(result(:,1),power((result(:,var)-kineResult6(:,var)),2),'r -.')
hold on
plot(result(:,1),power((result(:,var)-kineResult3(:,var)),2),'g -.')
title('y�ᾫ��')
hold off


