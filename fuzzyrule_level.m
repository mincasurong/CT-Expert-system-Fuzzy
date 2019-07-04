function [num,sum_num,mu,mu_num] = fuzzyrule_specific(a,b,Level_mbs0)

num = zeros(5,1);
mu_num = zeros(5,1);
mu = zeros(1,1);

while(1), k=1;  num(k) = a(1) * b(1); 
    mu_num(k,1) = num(k) * Level_mbs0(1); 
    break;
end

while(1), k=2;  num(k) = a(2) * b(2); 
    mu_num(k,1) = num(k) * Level_mbs0(2); 
    break;
end

while(1), k=3;  num(k) = a(3) * b(3); 
    mu_num(k,1) = num(k) * Level_mbs0(3); 
    break;
end

while(1), k=4;  num(k) = a(4) * b(4); 
    mu_num(k,1) = num(k) * Level_mbs0(4); 
    break;
end

while(1), k=5;  num(k) = a(5) * b(5); 
    mu_num(k,1) = num(k) * Level_mbs0(5); 
    break;
end


% a, b, c ,d
sum_num=sum(num);

if sum_num == 0,                                           % Undefined membership function -> Maintain
    mu = [1];
    fprintf('\n<Error Mode> Decision = Maintain \n')
    
else, mu = sum(mu_num)/sum_num; end
end
