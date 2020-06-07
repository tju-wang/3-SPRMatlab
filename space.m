%% 6杆  正解 
%求解工作空间  一定要求解设备实际的工作空间才有意义
%求得模型的正解之后，需要求真实模型的各个点的坐标  C点  P副  P副长度使用C点和配重块之间的长度来限制
%干涉条件  1.P副长度 最长最短  2.支架干涉  怎么算？？  3.动平台与桌面  配重块与桌面之间的干涉
%标记干涉种类  位操作？  在一定干涉范围内修改标记位  最后画对应点的颜色
%计算和分析分开整

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
R = RX_alpha*RY_beta*RZ_gama;

[ux,vx,wx,uy,vy,wy,uz,vz,wz] = deal(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3));
%R  =[ux,vx,wx;uy,vy,wy;uz,vz,wz]

%修改坐标系  x轴正方向向  右   y正方向向上   z 正方向由静平台指向动平台
A1_o = [0;a;0];
A2_o = [-(3.^(1/2)/2)*a;-1/2*a;0];
A3_o = [(3^(1/2)/2)*a;-1/2*a;0];

B1 = [0;b;20];
B2 = [-(3^(1/2)/2)*b;-1/2*b;20];
B3 = [(3^(1/2)/2)*b;-1/2*b;20];
% B1 = [0;b;0];
% B2 = [-(3^(1/2)/2)*b;-1/2*b;0];
% B3 = [(3^(1/2)/2)*b;-1/2*b;0];

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
a = 41.56; %动平台外接圆半径
b = 80;
alpha = 0.3;
beta = -0.3;
gama = 0.3;
Z0 = 250;
X0 = 10;
Y0 = 10;
argu = eval(argu);
numm=1
sprspace = []
errflag = 0;
color = 0;
qmin = 100; qmax = 295; qdiff = 2;
tic
rember_err = []
for q1 = qmin:qdiff:qmax
for q2 = qmin:qdiff:qmax
   for q3 = qmin:qdiff:qmax
       %先对杆长情况进行判断
       color = 0;
       errq(1) = abs(q1-q2);
       errq(2) = abs(q1-q3);
       errq(3) = abs(q2-q3);
       if max(errq(:))>=sqrt(3)*(b-a)    %三角形 两边之差小于第三边   大于or等于时，不计算工作空间
            rember_err(numm,:) = [numm q1 q2 q3];
           continue;
       end
    err = 1;err_2 = 10;
    Fi = eval(Fii);
    diff_argu = [0 0 0 0 0 0];
    num = 0;
    while(err>1.0e-4 || err_2>1.0e-4)
        argu = (argu + diff_argu);
        alpha = argu(1);
        beta = argu(2);
        gama = argu(3);
        X0 = argu(4);
        Y0 = argu(5);
        Z0 = argu(6);
        Fi = eval(Fii);
        diff_argu = (-inv(eval(J))*Fi)';
        err = diff_argu(1)*diff_argu(1)+diff_argu(2)*diff_argu(2);
        err_2 = diff_argu(6)*diff_argu(6);
        num = num+1;
        if num > 50
            errflag = 1;
            break;
        end
    end
    if errflag ~= 1
        %求解出结果  对结果是否是机构的工作空间进行验算
        sprspace(numm,:) = [numm q1 q2 q3 argu color err err_2 num];
        numm = numm+1
    else
        errflag = 0;
        errq = [q1 q2 q3]
    end
    
   end    
end
end



time = toc
sprspace

%% 分析数据  改变颜色
num = max(sprspace(:))
numm = 1;
spaceRed = []
redNum = 1;
spaceBlue = []
blueNum = 1;
for i = 1:num
    if (sprspace(numm,2) == qmin) || sprspace(numm,3) == qmin || sprspace(numm,4) == qmin
        sprspace(numm,11) = 1;
        spaceRed(redNum,:) = [redNum sprspace(numm,:)];
        redNum = redNum+1;
    end
    if (sprspace(numm,2) == qmax) || sprspace(numm,3) == qmax || sprspace(numm,4) == qmax
       sprspace(numm,11) = bitor(sprspace(numm,11),2); 
        spaceBlue(blueNum,:) =[blueNum sprspace(numm,:)];
        blueNum = blueNum+1;
    end
    numm = numm+1;
    
end


%区分杆长限制边界的颜色
num = max(sprspace(:))
numm = 1;
figure(1)
hold on
for i = 1:num    
    if sprspace(numm,11)==1
        scatter3(sprspace(numm,8),sprspace(numm,9),sprspace(numm,10),'b o')
    elseif sprspace(numm,11)>=2
        scatter3(sprspace(numm,8),sprspace(numm,9),sprspace(numm,10),'g o')
    else
        scatter3(sprspace(numm,8),sprspace(numm,9),sprspace(numm,10),'r o')
    end
    numm = numm+1
end

%scatter3(sprspace(:,8),sprspace(:,9),sprspace(:,10),'r o')
datasize = size(spaceRed(:,9))
grid_diff = 300/datasize(1)
[Mx My] = meshgrid(-150:grid_diff:150);
Mz = griddata(spaceRed(:,9),spaceRed(:,10),spaceRed(:,11),Mx,My)

figure(2)
mesh(Mx,My,Mz)
hold on
Hz = griddata(spaceBlue(:,9),spaceBlue(:,10),spaceBlue(:,11),Mx,My)
mesh(Mx,My,Hz)


