%% 输出Jacobian

f = fopen('C:\Users\Mrwang\Desktop\Haptic Project\3-SPR Kinematic\MATLAB\J.txt','wt')
for kk = 1:1:6
    for i = 1:1:6
         fprintf(f,'\n\n %d,%d  \n',kk,i)
         fprintf(f,'%s\n',char(J(kk,i)))
    end
end

fclose(f)
%% 输出Fi
f = fopen('C:\Users\Mrwang\Desktop\Haptic Project\3-SPR Kinematic\MATLAB\Fi.txt','wt')
for kk = 1:1:6
    fprintf(f,'\n\n %d  \n',kk)
    fprintf(f,'%s\n',char(Fii(kk)))
end

fclose(f)



%%
alpha = 0.2;
beta = 0.32;
gama = 0.2;
X0 = 20;
Y0 = 40;
Z0 = 200;
a = 41.56;
b = 80;
s = 62.0;
KK = vpa(eval(J),10)
% 
% Matrix = [ 7,2,4,5,6,8;
%            8,9,6,8,7,1;
%            8,7,7,7,9,1;
%            2,5,8,9,7,6;
%            15,85,44,12,3,6;
%            4,8,7,5,6,2]
q1 = 200;  q2 = 200; q3 = 205;
FF = eval(Fi)

%% 解C1点坐标
% clear q1 q2 q3 M1 M2 M3 M4 K1 K2 K3 K4 ca cb cc alpha beta gama X0 Y0 Z0
syms q1 q2 q3
syms Xc Yc Zc XB YB ZB XA YA ZA Xi Yi Zi Xc1 Yc1 Zc1 ca cb cc real

Xc = Xc1; Yc = Yc1; Zc = Zc1; XB = B1(1); YB = B1(2); ZB = B1(3); XA = A1(1); YA = A1(2); ZA = A1(3);
Xi = A2(1)-A3(1);  Yi = A2(2)-A3(2);  
Zi = A2(3)-A3(3);

% XB = B2(1); YB = B2(2); ZB = B2(3); XA = A2(1); YA = A2(2); ZA = A2(3);
% Xi = A1(1)-A3(1);  Yi = A1(2)-A3(2); Zi;
% Zi = A1(3)-A3(3);
% XB = B3(1); YB = B3(2); ZB = B3(3); XA = A3(1); YA = A3(2); ZA = A3(3);
% Xi = A1(1)-A2(1);  Yi = A1(2)-A2(2); 
% Zi = A1(3)-A2(3);


M1 = Xi*XB+Yi*YB+Zi*ZB;
K1 = 2*(YA-YB)-2*(XA-XB)*Yi/Xi;
K2 = 2*(ZA-ZB)-2*(XA-XB)*Zi/Xi;
M2 = 2*(XA-XB)*M1/Xi+(XB^2+YB^2+ZB^2-(XA^2+YA^2+ZA^2)+q3^2-2*s^2);
K4 = -K1/K2;
M4 = -M2/K2;
K3 = -(Zi*K4+Yi)/Xi;
M3 = (M1-Zi*M4)/Xi;


ca = (K3^2+K4^2+1);
cb = (-(2*XB*K3+2*YB+2*ZB*K4)+2*K3*M3+2*K4*M4);
cc = (M3^2+M4^2 -(2*XB*M3+2*ZB*M4)+XB^2+YB^2+ZB^2-s^2); 
X = K3*Yc+M3
Z = Yc*K4+M4


f = fopen('C:\Users\Mrwang\Desktop\Haptic Project\3-SPR Kinematic\MATLAB\solve.txt','wt');
fprintf(f,'ca\n%s\n\n',char(ca));
fprintf(f,'cb\n%s\n\n',char(cb));
fprintf(f,'cc\n%s\n\n',char(cc));
fprintf(f,'CY\n%s\n\n',char(X));
fprintf(f,'CZ\n%s\n\n',char(Z));
fprintf(f,'vB1C1\n%s\n\n',char(vB1C1));
fprintf(f,'vAxi_A1\n%s\n\n',char(vAxi_A1));
fprintf(f,'vA1C1\n%s\n\n',char(vA1C1));



fprintf(f,'先替换3^(1/2)  再找 ^ 然后替换 ^2');

% ca = eval(K3^2+K4^2+1);
% cb = eval(-(2*XB*K3+2*YB+2*ZB*K4)+2*K3*M3+2*K4*M4);
% cc = eval(M3^2+M4^2 -(2*XB*M3+2*ZB*M4)+XB^2+YB^2+ZB^2-s^2); 

% 
% if real((((cb.^2-4*ca*cc)))) > 0
%     Yc = real(((-cb+sqrt(cb.^2-4*ca*cc))/2/ca));
%     Xc = eval(K3*Yc+M3);
%     Zc = eval(Yc*K4+M4);
% % Xc = eval((M1-Yc*Yi-Zc*Zi)/Xi)
%     Xc1 = Xc
%     Yc1 = Yc
%     Zc1 = Zc
% % K3 = eval(K3)
% % K5 = eval(K5)
% % K4 = eval(K4)
% % K6 = eval(K6)
% else
%     sprintf('不满足方程求解条件')
% end

% fc1_1 = eval(norm(vB1C1)-s);
% fc1_2 = eval(vB1C1 * vAxi_A1');
% fc1_3 = eval(norm(vA1C1)-sqrt(q1^2-s^2));
% fc1 = [fc1_1 fc1_2 fc1_3]

% fc2_1 = eval(norm(vB2C2)-s);
% fc2_2 = eval(vB2C2 * vAxi_A2');
% fc2_3 = eval(norm(vA2C2)-sqrt(q2^2-s^2));
% fc2 = [fc2_1 fc2_2 fc2_3]

%% C3点坐标
syms Zi;

ca = (K4.^2+K6.^2+1)
cb = (2*K4*M4+2*K6*M6-2*XB-2*K4*YB-2*ZB*K6)
cc = (M4.^2+M6.^2-2*YB*M4-2*ZB*M6+XB.^2+YB.^2+ZB.^2-s.^2)
Yc = ((K4*Xc+M4))
Zc = ((K6*Xc+M6))

vB3C3
vAxi_A3
vA3C3

 f = fopen('C:\Users\Mrwang\Desktop\Haptic Project\3-SPR Kinematic\MATLAB\solve.txt','wt');
fprintf(f,'ca\n%s\n\n',char(ca));
fprintf(f,'cb\n%s\n\n',char(cb));
fprintf(f,'cc\n%s\n\n',char(cc));
fprintf(f,'CY\n%s\n\n',char(Yc));
fprintf(f,'CZ\n%s\n\n',char(Zc));
fprintf(f,'vB1C1\n%s\n\n',char(vB3C3));
fprintf(f,'vAxi_A1\n%s\n\n',char(vAxi_A3));
fprintf(f,'vA1C1\n%s\n\n',char(vA3C3));
Zi = A1(3)-A2(3);
fprintf(f,'Zi\n%s\n\n',char(Zi));

fprintf(f,'先替换3^(1/2)  再找 ^ 然后替换 ^2');

ca = eval(K4.^2+K6.^2+1);
cb = eval(2*K4*M4+2*K6*M6-2*XB-2*K4*YB-2*ZB*K6);
cc = eval(M4.^2+M6.^2-2*YB*M4-2*ZB*M6+XB.^2+YB.^2+ZB.^2-s.^2);
% K3 = eval(K3)
% K5 = eval(K5)
% K4 = eval(K4)
% K6 = eval(K6)
if real((((cb.^2-4*ca*cc)))) > 0
    Xc = real(((-cb-sqrt(cb.^2-4*ca*cc))/2/ca))
    Yc = real(eval(K4*Xc+M4))
    Zc = real(eval(K6*Xc+M6))
    Xc3 = Xc; Yc3 = Yc; Zc3 = Zc;
else
    sprintf('不满足方程求解条件')
end
fc3_1 = eval(norm(vB3C3)-s);
fc3_2 = eval(vB3C3 * vAxi_A3');
fc3_3 = eval(norm(vA3C3)-sqrt(q3^2-s^2));
fc3 = [fc3_1 fc3_2 fc3_3]

%% 解方程
A = [6,-2,3;3,-1,5;2,1,5];
b = [1;2;3]
B = [A b]

X = A\b
%% 求C2
syms Zi

ca = (K4.^2+K6.^2+1);
cb = (2*K4*M4+2*K6*M6-2*XB-2*K4*YB-2*ZB*K6);
cc = (M4.^2+M6.^2-2*YB*M4-2*ZB*M6+XB.^2+YB.^2+ZB.^2-s.^2);


f = fopen('C:\Users\Mrwang\Desktop\Haptic Project\3-SPR Kinematic\MATLAB\solve.txt','wt');
fprintf(f,'ca\n%s\n\n',char(ca));
fprintf(f,'cb\n%s\n\n',char(cb));
fprintf(f,'cc\n%s\n\n',char(cc));
fprintf(f,'CY\n%s\n\n',char(K4*Xc+M4));
fprintf(f,'CZ\n%s\n\n',char(K6*Xc+M6));
fprintf(f,'vB1C1\n%s\n\n',char(vB2C2));
fprintf(f,'vAxi_A1\n%s\n\n',char(vAxi_A2));
fprintf(f,'vA1C1\n%s\n\n',char(vA2C2));
Zi = A1(3)-A3(3);
fprintf(f,'Zi\n%s\n\n',char(Zi));


fprintf(f,'先替换3^(1/2)  再找 ^ 然后替换 ^2');


ca = eval(K4.^2+K6.^2+1);
cb = eval(2*K4*M4+2*K6*M6-2*XB-2*K4*YB-2*ZB*K6);
cc = eval(M4.^2+M6.^2-2*YB*M4-2*ZB*M6+XB.^2+YB.^2+ZB.^2-s.^2)


if (((cb.^2-4*ca*cc))) > 0
    Xc = real(((-cb+sqrt(cb.^2-4*ca*cc))/2/ca));
    Yc = eval(K4*Xc+M4);
    Zc = eval(K6*Xc+M6);
    Xc2 = Xc
    Yc2 = Yc
    Zc2 = Zc
else
    sprintf('不满足方程求解条件')
end
fc2_1 = eval(norm(vB2C2)-s);
fc2_2 = eval(vB2C2 * vAxi_A2');
fc2_3 = eval(norm(vA2C2)-sqrt(q2^2-s^2));
fc2 = [fc2_1 fc2_2 fc2_3]



















