function [P] = clc_P(n,C,L) 
%传入参数 方向向量 坐标起始点（点C1 C2 or C3） 长度L  求出P坐标
    n = n/norm(n);
    if n(3) > 0
        n = -n;
    end
    P = C + n*L;
end