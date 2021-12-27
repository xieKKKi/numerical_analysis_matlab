clear;

a = 0;
b = pi/3;
c = 0;
d = pi/6;
e = 0.000005;
k = 1;

while(1)
    t1 = T3(2^k, 2^k, a,b,c,d);
    t2 = T3(2^(k+1), 2^(k+1), a,b,c,d);
    if(abs(t2-t1) <= e)
        fprintf('计算结果: %f \n 等分数: %d \n', t2, 2^(k+1+2));
        break;
    end
    k=k+1;
end

function [z] = f(x, y)
%被积函数
    z = tan(x^2+y^2);
end

function [t] = T(m, n, a, b, c, d)
% 课本公式5.7.7 
% 将x对应的区间[a,b]作m等分，y对应的[c,d]作n等分
    t = 0;
    h = (b-a)/m;
    k = (d-c)/n;
        for i = 0:m-1
            xi = a+i*h;
            xi1 = a+(i+1)*h;
            for j = 0:n-1
                yj = c+j*k;
                yj1 = c+(j+1)*k;
                t = t + (h*k)/4 * (f(xi,yj) + f(xi, yj1) + f(xi1,yj) + f(xi1, yj1));
            end
        end
end

% T1、T2、T3为使用Richardson外推思想得到的求积公式
function [t] = T1(m, n, a, b, c, d)
    t = 4/3*T(2*m, 2*n, a, b, c, d) - 1/3*T(m, n, a, b, c, d);
end

function [t] = T2(m, n, a, b, c, d)
    t = 16/15*T1(2*m, 2*n, a, b, c, d) - 1/15*T1(m, n, a, b, c, d);
end

function [t] = T3(m, n, a, b, c, d)
    t = 64/63*T2(2*m, 2*n, a, b, c, d) - 1/63*T2(m, n, a, b, c, d);
end