x=[0 1 2 3 4 5 6 7 8 9 10];
y=[2.51 3.30 4.04 4.70 5.22 5.54 5.78 5.40 5.57 5.70 5.80];
y0=0.8 ;          %  S'(x0)=f'(x0)=y0   
yn=0.2;          %  S'(xn)=f'(xn)=yn

x1=0.5:1:9.5;
s1=resultCubicSplineInterpolation(x,y,x1,y0,yn);
Answer = sprintf(' \n S(i+0.5): \n i=0:  %.6f \n i=1:  %.6f \n i=2:  %.6f \n i=3:  %.6f \n i=4:  %.6f \n i=5:  %.6f \n i=6:  %.6f \n i=7:  %.6f \n i=8:  %.6f \n i=9:  %.6f \n',s1(1), s1(2),s1(3),s1(4),s1(5),s1(6),s1(7),s1(8),s1(9),s1(10));
fprintf(Answer);

x0=0:0.01:10;
s=resultCubicSplineInterpolation(x,y,x0,y0,yn);
plot(x0,s)        %绘制第一边界条件插值函数图像
hold on
grid on
plot(x,y,'o')
%axis([0.2 0.55 0.4 0.75])
xlabel('X'), ylabel('Y')
title('插值点与三次样条函数') 
legend('三次样条插值点坐标','插值点')

function [D,h,A,g,M]=cubicSplineInterpolation1(X,Y,y0,yn)
%        自然边界条件的三次样条函数(第一种边界条件)
%        此函数为M值求值函数
%        D,h,A,g,M输出量分别为系数矩阵D，插值宽度h，差商表A，g值,M值 
         n=length(X); 
         A=zeros(n,n);A(:,1)=Y';D=zeros(n,n);g=zeros(n,1);
         for  j=2:n
            for i=j:n
                A(i,j)=(A(i,j-1)- A(i-1,j-1))/(X(i)-X(i-j+1));
            end
         end
         
         for i=1:n-1
             h(i)=X(i+1)-X(i);
         end
         for i=1:n
             D(i,i)=2; 
             D(1,2)=1;
             D(n,n-1)=1;
             if (i==1)
                 g(i,1)=6/h(i)*(A(2,2)-y0); 
             elseif (i==n) 
                     g(i,1)=6/h(i-1)*(yn-A(i,2));
             else 
                 g(i,1)=(6/(h(i-1)+h(i)))*(A(i+1,2)-A(i,2));
             end
           
         end  
         for i=1:n-2
             u(i)=h(i)/(h(i)+h(i+1));
             n(i)=1-u(i);  
             D(i+1,i+2)=n(i);
             D(i+1,i)=u(i);             
         end
         M=D\g;
end

function s=resultCubicSplineInterpolation(X,Y,x,y0,yn)
%        三次样条插值函数第一类型代码 
%        s函数表示三次样条插值函数插值点对应的函数值
%        根据三次样条参数函数求出的D,h,A,g,M
%        x表示求解插值点函数点，X为已知插值点        
         [D,h,A,g,M]=cubicSplineInterpolation1(X,Y,y0,yn);
         n=length(X); m=length(x);    
         for i=1:n-1
             a0 = Y(i);
             a1 = A(i+1,1)-A(i,1)-(M(i,1)/3+M(i+1,1)/6)*h(i);
             a2 = (M(i,1))/2;
             a3 = (M(i+1,1)-M(i,1))/(6*h(i));
             fprintf('S(x) = %.6f + %.6f * (x-%d) + %.6f * (x-%d)^2 + %.6f * (x-%d)^3', a0,a1, X(i),a2, X(i),a3, X(i));
             fprintf('  ( %d<x<%d )\n', X(i),X(i+1));
        end
         
         for t=1:m
            for i=1:n-1
               if (x(t)<=X(i+1))&&(x(t)>=X(i))
                  a0 = Y(i);
                  a1 = A(i+1,1)-A(i,1)-(M(i,1)/3+M(i+1,1)/6)*h(i);
                  a2 = (M(i,1))/2;
                  a3 = (M(i+1,1)-M(i,1))/(6*h(i));
                  s(t)=a0+a1*(x(t)-X(i))+a2*(x(t)-X(i))^2+a3*(x(t)-X(i))^3;
                  break;
               else
                   s(t)=0; 
               end
            end
         end
end
