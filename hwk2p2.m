function [x] = hwk2p2(x0,n)

syms x1 x2;
fenton(x1,x2) = (12 + x1^2 + ((1 + x2^2)/(x1^2)) + (((x1^2)*(x2^2) + 100)/((x1*x2)^4)))*(1/10);

path = [];
x = x0;

for i = 1:n
    path = [path;x];
    B = double(subs(hessian(fenton),[x1,x2],x));
    [~,PD_flag] = chol(B);
    count = 0;
    
    while PD_flag ~= 0
        B = B + (count + rand(1))*eye(2);
        [~,PD_flag] = chol(B);
        count = count + 1;
    end
    
    G = double(subs(gradient(fenton),[x1,x2],x));
    p = -B\G;
    
    alpha = 1;
    c = 0.5;
    rho = 0.75;
    
    step_val = subs(fenton,[x1,x2],x + alpha*p');
    requirement = subs(fenton,[x1,x2],x) + c*alpha*G'*p;
    select_alpha = requirement - step_val;
    
    while select_alpha < 0
        alpha = rho*alpha;
        step_val = subs(fenton,[x1,x2],x + alpha*p');
        requirement = subs(fenton,[x1,x2],x) + c*alpha*G'*p;
        select_alpha = requirement - step_val;
    end
    
    x = x + alpha * p';
    
end
      
    

    
    

