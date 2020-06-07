%导入数据后 分析
%和计算力的部分  1.修改了y轴方向  2.修改了求C点时解方程的符号
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
% A1_o = [0;-a;0];
% A2_o = [(3.^(1/2)/2)*a;1/2*a;0];
% A3_o = [-(3^(1/2)/2)*a;1/2*a;0];
% 
% B1 = [0;-b;20];
% B2 = [(3^(1/2)/2)*b;1/2*b;20];
% B3 = [-(3^(1/2)/2)*b;1/2*b;20];
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

vA1C1 = vector(A1',C1); 
vA2C2 = vector(A2',C2);
vA3C3 = vector(A3',C3);

vB1C1 = vector(B1,C1);
vB2C2 = vector(B2,C2);
vB3C3 = vector(B3,C3);

vAxi_A1 = A2' - A3';
vAxi_A2 = A1' - A3';
vAxi_A3 = A1' - A2';

s = 62;  %测量得 58.5mm + 6.5/2
a = 41.56; %动平台外接圆半径
b = 80;
rq_min = 90
rq_max = 275
err_rq = 7  %diff为5
err_C = 2
spaceData = []  %分析过后的数据集合
spaceNdata = []
%sprspace(numm,:) = [numm q1 q2 q3 argu color err err_2 num];
%% 求解C1  C2  C3的坐标
syms Xc Yc Zc XB YB ZB XA YA ZA Xi Yi Zi real
load ./spaceData/space6.31.mat
rq = []
spNum = 1;
for num = 1:size(sprspace,1)  %从1到sprspace的列数
    
%A1需要更新
alpha = sprspace(num,5);
beta = sprspace(num,6);
gama = sprspace(num,7);
X0 = sprspace(num,8);
Y0 = sprspace(num,9);
Z0 = sprspace(num,10);
q1 = sprspace(num,2);
q2 = sprspace(num,3);
q3 = sprspace(num,4);

%求C1  C2  C3
Xc = Xc1; Yc = Yc1; Zc = Zc1; XB = B1(1); YB = B1(2); ZB = B1(3); XA = A1(1); YA = A1(2); ZA = A1(3);
Xi = A2(1)-A3(1);   Yi = A2(2)-A3(2);    Zi = A2(3)-A3(3);

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
    
    Xc1 = Xc;
    Yc1 = Yc;
    Zc1 = Zc;
else
    sprintf('C1 不满足方程求解条件')
    sprspace(num,:)
    break;
     continue;
end

fc1_1 = eval(norm(vB1C1)-s);
fc1_2 = eval(vB1C1 * vAxi_A1');
fc1_3 = eval(norm(vA1C1)-sqrt(q1^2-s^2));
fc1 = [fc1_1 fc1_2 fc1_3];
if(fc1_1+fc1_2+fc1_3 > 0.0001)
    sprintf('C1 求解有误')
end

%求解C2
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
    sprintf('C2 不满足方程求解条件')
    continue;
end

fc2_1 = eval(norm(vB2C2)-s);
fc2_2 = eval(vB2C2 * vAxi_A2');
fc2_3 = eval(norm(vA2C2)-sqrt(q2^2-s^2));
fc2 = [fc2_1 fc2_2 fc2_3];

if(fc2_1+fc2_2+fc2_3 > 0.0001)
    sprintf('C2 求解有误')
end
%求解C3
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
    sprintf('C3 不满足方程求解条件')
     continue;
end

fc3_1 = eval(norm(vB3C3)-s);
fc3_2 = eval(vB3C3 * vAxi_A3');
fc3_3 = eval(norm(vA3C3)-sqrt(q3^2-s^2));
fc3 = [fc3_1 fc3_2 fc3_3];
if(fc3_1+fc3_2+fc3_3 > 0.0001)
    sprintf('C3 求解有误')
end
num

color = 0;
rq(num,:) = [norm(eval(vA1C1)) norm(eval(vA2C2)) norm(eval(vA3C3))];
rq(num,:);

if abs(rq(num,1) - rq_min)<err_rq || abs(rq(num,2) - rq_min)<err_rq || abs(rq(num,3) - rq_min)<err_rq
    color = bitor(color,2^1);
end
if abs(rq(num,1) - rq_max)<err_rq || abs(rq(num,2) - rq_max)<err_rq || abs(rq(num,3) - rq_max)<err_rq
    color = bitor(color,2^2);
end
if abs(eval(B1(3)) - eval(C1(3)))<err_C || abs(eval(B2(3)) - eval(C2(3)))<err_C || abs(eval(B3(3)) - eval(C3(3)))<err_C
    color = bitor(color,2^3);
end
sprspace(num,11) = color;
spaceNdata(num,:) = [sprspace(num,1:11) rq(num,:) eval(C1) eval(C2) eval(C3)];
if (rq(num,1) > rq_min) && (rq(num,2) > rq_min) && (rq(num,3) > rq_min)
    if (rq(num,1) < rq_max+err_C) && (rq(num,2) < rq_max+err_C) && (rq(num,3) < rq_max+err_C)
        if (eval(B1(3)) <= eval(C1(3))) && (eval(B2(3)) <= eval(C2(3))) && (eval(B3(3)) <= eval(C3(3)))
            spaceData(spNum,:) = [spNum sprspace(num,2:11) rq(num,:) eval(C1) eval(C2) eval(C3)];
            spNum = spNum+1;
        end
       
    end
end
end  %end..数据总循环
%% 画图
%scatter3(sprspace(numm,8),sprspace(numm,9),sprspace(numm,10),'b o')



figure(20)
hold on  
scatter3(spaceData(:,8),spaceData(:,9),spaceData(:,10),400,'o','filled','MarkerFaceColor',[162/255 205/255 90/255])


%ones(length(msdata(i,1)))*[0.1,0.3,0]
%scatter3(sprspace(numm,8),sprspace(numm,9),sprspace(numm,10),400,'o','filled','c',[1 1 0])
% rqNum = 1;
% for i=1:size(spaceNdata)
%    if (spaceNdata(i,12) - rq_max)>0 || (spaceNdata(i,13) - rq_max)>0 || (spaceNdata(i,14) - rq_max)>0
%       rqLength(rqNum,:)  = spaceNdata(i,:);
%       rqNum = rqNum+1;
%    end
%    i
% end
% figure(25)
% hold on
% scatter3(rqLength(:,8),rqLength(:,9),rqLength(:,10),'b o')
% for i=1,size(rqLength,1)
%     if spaceNdata(rqLength(i,1),11) == 4 || spaceNdata(rqLength(i,1),11) == 12
%     i;
%     else
%         i
%         sprintf('error')
%         break;
%     end
% 
% end


figure(21)
hold on 
get = 0;
get = 2^1;
rq_min_num = 0
for i = 1:size(spaceData,1)
    numm = i
 
   if bitand(spaceData(numm,11),get) == get
       scatter3(spaceData(numm,8),spaceData(numm,9),spaceData(numm,10),'b o')
       rq_min_num = rq_min_num + 1;
   end
end

figure(22)
hold on 
get = 0;
get = 2^2;
rq_max_num = 0
for i = 1:size(spaceData,1)
    numm = i
   if bitand(spaceData(numm,11),get) == get
       scatter3(spaceData(numm,8),spaceData(numm,9),spaceData(numm,10),'g o')
       rq_max_num = rq_max_num + 1;
   end
end

figure(23)
hold on 
get = 0;
get = 2^3;
err_C_num = 0;
for i = 1:size(spaceData,1)
    numm = i

   if bitand(spaceData(numm,11),get) == get
       scatter3(spaceData(numm,8),spaceData(numm,9),spaceData(numm,10),'r o')
       err_C_num = err_C_num+1;
   end
end

figure(24)
hold on 
get = 0;
get = 2^3;
err_C_num = 0;
for i = 1:size(spaceData,1)
    numm = i

   if spaceData(numm,11) ~= 0
       if bitand(spaceData(numm,11),2) == 2
           scatter3(spaceData(numm,8),spaceData(numm,9),spaceData(numm,10),'g o')
       end
       if bitand(spaceData(numm,11),4) == 4
           scatter3(spaceData(numm,8),spaceData(numm,9),spaceData(numm,10),'r o')
       end
      if bitand(spaceData(numm,11),8) == 8
           scatter3(spaceData(numm,8),spaceData(numm,9),spaceData(numm,10),'y o')
       end
       err_C_num = err_C_num+1;
   end
end



