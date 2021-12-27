clear;

x0 = 0;  % x初值
xn = 1.5; % x终值
y0 = 3; % y初值
h = 0.1; % 步长

[x_precise,y_precise] = f_precise(x0, xn, y0, h);   % 精确解

[x,y] = RK4(x0, xn, y0, h); % 4阶龙格库塔法
error = cal_error(y, y_precise);
disp('RK4: ');
disp('---------------------------------------------------------------------------');
disp('i                xi                        yi                    y(xi)                 y(xi)-yi');
disp('---------------------------------------------------------------------------');
for i = 1:size(x)
    fprintf('%d       %.8f       %.8f       %.8f       %.8f\n', i, x(i), y(i), y_precise(i), error(i));
end
disp('---------------------------------------------------------------------------');

[x,y] = AB4(x, y ,h); % AB4
error = cal_error(y, y_precise);
disp('AB4: ');
disp('---------------------------------------------------------------------------');
disp('i                xi                        yi                    y(xi)                 y(xi)-yi');
disp('---------------------------------------------------------------------------');
for i = 1:size(x)
    fprintf('%d       %.8f       %.8f       %.8f       %.8f\n', i, x(i), y(i), y_precise(i), error(i));
end
disp('---------------------------------------------------------------------------');

[x,y] = AB4_AM4(x, y ,h); % AB4_AM4
error = cal_error(y, y_precise);
disp('AB4_AM4: ');
disp('---------------------------------------------------------------------------');
disp('i                xi                        yi                    y(xi)                 y(xi)-yi');
disp('---------------------------------------------------------------------------');
for i = 1:size(x)
    fprintf('%d       %.8f       %.8f       %.8f       %.8f\n', i, x(i), y(i), y_precise(i), error(i));
end
disp('---------------------------------------------------------------------------');

[x,y] = AB4_AM4_Richardson(x, y ,h); % AB4_AM4_Richardson
error = cal_error(y, y_precise);
disp('AB4_AM4_Richardson: ');
disp('---------------------------------------------------------------------------');
disp('i                xi                        yi                    y(xi)                 y(xi)-yi');
disp('---------------------------------------------------------------------------');
for i = 1:size(x)
    fprintf('%d       %.8f       %.8f       %.8f       %.8f\n', i, x(i), y(i), y_precise(i), error(i));
end
disp('---------------------------------------------------------------------------');

function[z] = f(x, y)
    z = -x*x*y*y;
end

function[x, y] = f_precise(x0, xn, y0, h)
    n = (xn-x0)/h;
    x = zeros(n+1, 1);
    y = zeros(n+1, 1);
    x(1) = x0;
    y(1) = y0;
    for i = 1:n
        x(i+1) = x(i)+h;
        y(i+1) = 3/(1+x(i+1)^3);
    end
end

function[error] = cal_error(y, y_precise)
    error = zeros(size(y));
    for i = 1:size(y)
        error(i) = y_precise(i) - y(i);
    end
end

function[x, y] = RK4(x0, xn, y0, h)
    n = (xn-x0)/h;
    x = zeros(n+1, 1);
    y = zeros(n+1, 1);
    x(1) = x0;
    y(1) = y0;
    for i = 1:n
        x(i+1) = x(i)+h;
        k1 = f(x(i), y(i));
        k2 = f(x(i)+h/2, y(i)+h*k1/2);
        k3 = f(x(i)+h/2, y(i)+h*k2/2);
        k4 = f(x(i)+h, y(i)+h*k3);
        y(i+1) = y(i) + (k1+2*k2+2*k3+k4)*h/6;
    end
end

function[x, y] = AB4(x ,y ,h)
    for i = 4:size(x)-1
        y(i+1) = y(i) + h/24*(55 * f(x(i),y(i)) - 59 * f(x(i-1),y(i-1)) + 37 * f(x(i-2),y(i-2)) - 9 * f(x(i-3),y(i-3)));
    end
end

function[x, y] = AB4_AM4(x ,y ,h)
yp = y;
    for i = 4:size(x)-1
        yp(i+1) = y(i) + h/24*(55 * f(x(i),y(i)) - 59 * f(x(i-1),y(i-1)) + 37 * f(x(i-2),y(i-2)) - 9 * f(x(i-3),y(i-3)));
        y(i+1) = y(i) + h/24*(9 * f(x(i+1),yp(i+1)) + 19 * f(x(i),y(i)) - 5 * f(x(i-1),y(i-1)) + f(x(i-2),y(i-2)));
    end
end

function[x, y] = AB4_AM4_Richardson(x ,y ,h)
yp = y;
yc = y;
    for i = 4:size(x)-1
        yp(i+1) = y(i) + h/24*(55 * f(x(i),y(i)) - 59 * f(x(i-1),y(i-1)) + 37 * f(x(i-2),y(i-2)) - 9 * f(x(i-3),y(i-3)));
        yc(i+1) = y(i) + h/24*(9 * f(x(i+1),yp(i+1)) + 19 * f(x(i),y(i)) - 5 * f(x(i-1),y(i-1)) + f(x(i-2),y(i-2)));
        y(i+1) = 251/270*yc(i+1) + 19/270*yp(i+1);
    end
end