function [P] = clc_P(n,C,L) 
%������� �������� ������ʼ�㣨��C1 C2 or C3�� ����L  ���P����
    n = n/norm(n);
    if n(3) > 0
        n = -n;
    end
    P = C + n*L;
end