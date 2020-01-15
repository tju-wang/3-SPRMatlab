function [ T ] = clc_liju(F,L)
%传入力的大小与对应点的力臂  求力矩
T = cross(L,F);
end