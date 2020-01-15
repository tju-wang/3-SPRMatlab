function  [V] = vector(M1,M2)
    V(1) = M2(1)-M1(1);
    V(2) = M2(2)-M1(2);
    V(3) = M2(3)-M1(3);
end