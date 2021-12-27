clear;

a = -1;
b = 1;
e = 0.5 * 10^(-7);
k = 1;

while(1)
    t1 = T3(2^k, a,b);
    t2 = T3(2^(k+1), a,b);
    k
    if(abs(t2-t1) <= e)
        fprintf('计算结果: %.9f \n', t2);
        break;
    end
    k=k+1;
end

function [z] = f(x)  %被积函数
    z = 1 / (1 + 100 * x * x);
end

function [t] = T(n, a, b)
% 将x对应的区间[a,b]作n等分
    t = 0;
    h = (b-a)/n;
    for k = 0:n-1
        xi = a+k*h;
        xi1 = a+(k+1)*h;
        t = t + f(xi) + f(xi1);
    end
    t = h /2 *t;
end


function [t] = T1(n, a, b)  % Simpson
    t = 4/3*T(2*n, a, b) - 1/3*T(n, a, b);
end

function [t] = T2(n, a, b)  % Cotes
    t = 16/15*T1(2*n, a, b) - 1/15*T1(n, a, b);
end

function [t] = T3(n, a, b)  % Romberg
    t = 64/63*T2(2*n, a, b) - 1/63*T2(n, a, b);
end