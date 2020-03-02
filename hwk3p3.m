function [x] = hwk3p3(x0,niter,eps)

[rows,~] = size(x0);
if rows == 1
    x0 = x0';
end

path = [];
norm_grad = Inf;
x = x0;

for i = 1:niter
    if norm_grad > eps
        path = [path;x'];
        [f,G,B] = rosenbrocknfgH(x);
        norm_grad = norm(G);
        [~,PD_flag] = chol(B);
        count = 0;
    
        while PD_flag ~= 0
            B = B + (count + rand(1))*eye(length(x0));
            [~,PD_flag] = chol(B);
            count = count + 1;
        end
        
        p = -B\G;
    
        alpha = 1;
        c = 0.5;
        rho = 0.75;
    
        [step_val,~,~] = rosenbrocknfgH(x + alpha*p);
        requirement = f + c*alpha*G'*p;
        select_alpha = requirement - step_val;
    
        while select_alpha < 0
            alpha = rho*alpha;
            [step_val,~,~] = rosenbrocknfgH(x + alpha*p);
            requirement = f + c*alpha*G'*p;
            select_alpha = requirement - step_val;
        end
    
        x = x + alpha * p;
    end
end



