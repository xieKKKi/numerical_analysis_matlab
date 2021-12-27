
clear;
syms x0 x1 A0 A1 a b;
% syms x0 x1 A0 A1;
% a = 0;
% b = 1;
[x, x1, A0, A1]=solve(A0+A1==b-a,...
    A0*x0+A1*x1==1/2*(b^2-a^2),...
    A0*(x0^2)+A1*(x1^2)==1/2*(b^3-a^3),...
    A0*(x0^3)+A1*(x1^3)==1/2*(b^4-a^4))
