function [num,sum_num,mu,mu_num] = fuzzyrule_specific(a,b,c,d,e,KH_mbs0,SPD_mbs0)

num = zeros(15,1);
mu_num = zeros(15,2);
mu = zeros(1,2);

while(1), k=1;  num(k) = a(1) * b(1) * c(1) * d(1); % TTOO Easy
    mu_num(k,1) = num(k) * KH_mbs0(5); % PH
    mu_num(k,2) = num(k) * SPD_mbs0(5); % PH
    break;
end
while(1), k=2;  num(k) = a(1) * b(1) * c(2) * d(2); % Too Easy
    mu_num(k,1) = num(k) * KH_mbs0(4); % PL
    mu_num(k,2) = num(k) * SPD_mbs0(5); % PH
    break;
end
while(1), k=3;  num(k) = a(1) * b(1) * c(3) * d(3); % Easy
    mu_num(k,1) = num(k) * KH_mbs0(3); % ZE
    mu_num(k,2) = num(k) * SPD_mbs0(4); % PL
    break;
end

while(1), k=4;  num(k) = a(2) * b(2) * c(1) * d(1); % Too Easy
    mu_num(k,1) = num(k) * KH_mbs0(4); % PL
    mu_num(k,2) = num(k) * SPD_mbs0(5); % PH    
    break;
end
while(1), k=5;  num(k) = a(2) * b(2) * c(2) * d(2); % Easy
    mu_num(k,1) = num(k) * KH_mbs0(3); % ZE
    mu_num(k,2) = num(k) * SPD_mbs0(4); % PL    
    break;
end
while(1), k=6;  num(k) = a(2) * b(2) * c(3) * d(3); % Timely easy
    mu_num(k,1) = num(k) * KH_mbs0(4);  % PL
    mu_num(k,2) = num(k) * SPD_mbs0(3); % ZE    
    break;
end

while(1), k=7;  num(k) = a(3) * b(3) * c(1) * d(1); % Easy
    mu_num(k,1) = num(k) * KH_mbs0(3); % ZE
    mu_num(k,2) = num(k) * SPD_mbs0(4); % PL    
    break;
end
while(1), k=8;  num(k) = a(3) * b(3) * c(2) * d(2); % Maintain
    mu_num(k,1) = num(k) * KH_mbs0(3); % ZE
    mu_num(k,2) = num(k) * SPD_mbs0(3); % ZE    
    break;
end
while(1), k=9;  num(k) = a(3) * b(3) * c(3) * d(3); % Hard
    mu_num(k,1) = num(k) * KH_mbs0(3); % ZE
    mu_num(k,2) = num(k) * SPD_mbs0(2); % NL    
    break;
end

while(1), k=10;  num(k) = a(4) * b(4) * c(1) * d(1); % Timely hard
    mu_num(k,1) = num(k) * KH_mbs0(2); % NL
    mu_num(k,2) = num(k) * SPD_mbs0(3); % ZE    
    break;
end
while(1), k=11;  num(k) = a(4) * b(4) * c(2) * d(2); % Hard
    mu_num(k,1) = num(k) * KH_mbs0(3); % ZE
    mu_num(k,2) = num(k) * SPD_mbs0(2); % NL    
    break;
end
while(1), k=12;  num(k) = a(4) * b(4) * c(3) * d(3); % Too Hard
    mu_num(k,1) = num(k) * KH_mbs0(2); % NL
    mu_num(k,2) = num(k) * SPD_mbs0(2); % NL    
    break;
end

while(1), k=13;  num(k) = a(5) * b(5) * c(1) * d(1); % Hard
    mu_num(k,1) = num(k) * KH_mbs0(3); % ZE
    mu_num(k,2) = num(k) * SPD_mbs0(2); % NL    
    break;
end
while(1), k=14;  num(k) = a(5) * b(5) * c(2) * d(2); % Too Hard
    mu_num(k,1) = num(k) * KH_mbs0(2); % NL
    mu_num(k,2) = num(k) * SPD_mbs0(2); % NL    
    break;
end
while(1), k=15;  num(k) = a(5) * b(5) * c(3) * d(3); % Hard
    mu_num(k,1) = num(k) * KH_mbs0(1); % ZE
    mu_num(k,2) = num(k) * SPD_mbs0(1); % NL    
    break;
end

% a, b, c ,d
sum_num=sum(num);

if sum_num == 0,                                           % Undefined membership function -> Maintain
    mu = [KH_mbs0(3) SPD_mbs0(3)];
    fprintf('\n<Error Mode> Decision = Maintain \n')
    
elseif sum_num > 0 && e(2) == max(e) %&& e(1) > 0.5     % Foot drop situation -> KH up
    mu = [e(2)*KH_mbs0(2) e(2)*SPD_mbs0(4)];
    if e(2)*KH_mbs0(2) < e(2)*KH_mbs0(3)                        % Avoid KH down in Foot drop
        mu(2) = KH_mbs0(3)+(KH_mbs0(3)-e(2)*KH_mbs0(2));
    end
    fprintf('\n<Foot Drop> Decision = KH small up \n')
elseif sum_num > 0 && e(1) == max(e) %&& e(1) > 0.5     % Foot drop situation -> KH up
    mu = [e(1)*KH_mbs0(1) e(1)*SPD_mbs0(5)];
    if e(1)*KH_mbs0(1) < KH_mbs0(2)                        % Avoid KH down in Foot drop
        mu(2) = KH_mbs0(2)+(KH_mbs0(2)-e(1)*KH_mbs0(1));
    end
    fprintf('\n<Foot Drop> Decision = KH large up \n')
    
else, for k=1:2, mu(1,k) = sum(mu_num(:,k))/sum_num; end  % Normal situation
end
