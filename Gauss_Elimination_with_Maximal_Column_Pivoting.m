clear;

R = [31, -13,0,0,0,-10,0,0,0;
    -13,35,-9,0,-11,0,0,0,0;
    0,-9,31,-10,0,0,0,0,0;
    0,0,-10,79,-30,0,0,0,-9;
    0,0,0,-30,57,-7,0,-5,0;
    0,0,0,0,-7,47,-30,0,0;
    0,0,0,0,0,-30,41,0,0;
    0,0,0,0,-5,0,0,27,-2;
    0,0,0,-9,0,0,0,-2,29];

V = [-15;27;-23;0;-20;12;-7;7;10];

I_1 = gauss_elimination(R,V);
I = vpa(I_1, 5)

function x=gauss_elimination(A,b)  	%输入矩阵A和列向量b，返回解向量x
[ni,~]=size(b);
if rank(A)~=rank([A,b])				%若系数矩阵秩和增广矩阵秩不相等，则无解
    fprintf('无解\n')
    return
else if rank(A)<ni							%若系数矩阵秩和增广矩阵秩相等，但是其秩小于未知量个数，则无穷解
        fprintf('无穷解\n')
        return
    else        
        for j=1:ni
            [~,ti]=max(A(j:ni,j));			%找出该列中按模最大的元素
            A([ti+j-1,j],:)=A([j,ti+j-1],:);%交换行
            for i=j+1:ni						%消去过程
                d=-A(i,j)/A(j,j);
                A(i,:)=A(i,:)+d*A(j,:);
                b(i)=b(i)+d*b(j);
            end
        end
        x=zeros(size(b));					%初始化解向量
        x(ni)=b(ni)/A(ni,ni);
        for i=ni-1:-1:1						%回代过程
            x(i)=(b(i)-sum(x.*A(i,:)'))/A(i,i);
        end
    end
end
end
