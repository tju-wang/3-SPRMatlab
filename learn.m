
[X,Y] = meshgrid(-3:1:3)

Z = peaks(X,Y)
%Z = []
subplot(3,1,1);
mesh(X,Y,Z);
subplot(3,1,2);
meshc(X,Y,Z);
subplot(3,1,3);meshz(X,Y,Z);
axis([-3 3 -3 3 -10 5]);
%%
% ����1
T = [1:100];D = [1:100]; K = rand(1,100);
% ��ֵ
[X,Y,Z]=griddata(T,D,K,linspace(min(T),max(T))',linspace(min(D),max(D)),'v4');
figure,surf(X,Y,Z);


