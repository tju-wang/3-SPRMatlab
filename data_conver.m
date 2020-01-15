function [ D ] = data_conver(M)
%传入力的作用点坐标 返回
dx = sqrt(M(2).^2+M(3).^2);
dy = sqrt(M(1).^2+M(3).^2);
dz = sqrt(M(1).^2+M(2).^2);
D = [dx,dy,dz];

end