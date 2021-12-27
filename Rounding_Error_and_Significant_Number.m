clear;

N1 = single(10^2);
true_result1 = single((3/2-1/N1-1/(N1+1))/2);
result1_1 = from_large_to_small(N1);
result1_2 = from_small_to_large(N1);
Answer1 = sprintf(' \n true_result1: %.6f \n result1_1:  %.6f \n result1_2:  %.6f \n',true_result1, result1_1,result1_2);
fprintf(Answer1);

N2 = single(10^4);
true_result2 = single((3/2-1/N2-1/(N2+1))/2);
result2_1 = from_large_to_small(N2);
result2_2 = from_small_to_large(N2);
Answer2 = sprintf(' \n true_result2: %.6f \n result2_1:  %.6f \n result2_2:  %.6f \n',true_result2, result2_1,result2_2);
fprintf(Answer2);

N3 = single(10^6);
true_result3 = single((3/2-1/N3-1/(N3+1))/2);
result3_1 = from_large_to_small(N3);
result3_2 = from_small_to_large(N3);
Answer3 = sprintf(' \n true_result3: %.6f \n result3_1:  %.6f \n result3_2:  %.6f \n',true_result3, result3_1,result3_2);
fprintf(Answer3);

function [result] = from_large_to_small(N)
%从大到小顺序
result = single(0);
for i=2:N
    result = result + 1/(i*i-1);
end
end

function [result] = from_small_to_large(N)
%从小到大顺序
result = single(0);
for i=N:-1:2
    result = result + 1/(i*i-1);
end
end