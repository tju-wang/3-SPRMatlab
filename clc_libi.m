
function [ D ] = clc_libi(R,L)
%传入力的作用点坐标 与 要对点取距的点坐标  返回力臂
dx = sqrt((L(2)-R(2)).^2+(L(3)-R(3)).^2);
dy = sqrt((L(1)-R(1)).^2+(L(3)-R(3)).^2);
dz = sqrt((L(1)-R(1)).^2+(L(2)-R(2)).^2);
D = [dx,dy,dz];

end

