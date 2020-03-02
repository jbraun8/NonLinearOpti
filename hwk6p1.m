function [x] = hwk6p1(x_in,tol,iter)

n = length(x_in);
norm_grad = Inf;
x = x_in;

[f,g,~] = rosenbrocknfgH(x);
H = eye(n)/norm(g);
nfunc = 1;

c = 0.5;


for i = 1:iter
    if norm_grad > tol
        p = -H*g;
        
        alpha = 1;
        rho = 0.75;
        
        [step_val,~,~] = rosenbrocknfgH(x + alpha*p);
        nfunc = nfunc + 1;
        requirement = f + c*alpha*g'*p;
        select_alpha = requirement - step_val;
    
        while select_alpha < 0
            alpha = rho*alpha;
            [step_val,~,~] = rosenbrocknfgH(x + alpha*p);
            nfunc = nfunc + 1;
            requirement = f + c*alpha*g'*p;
            select_alpha = requirement - step_val;
        end  
        
        x = x + alpha*p;
        
        gold = g;
        [f,g,~] = rosenbrocknfgH(x);
        nfunc = nfunc + 1;
        norm_grad = norm(g);
        
        s = alpha*p;
        y = g - gold;
        
        threshold = 0.2*s'*H*s;
        
        if s'*y >= threshold
            theta = 1;
        else
            theta = (0.8*s'*H*s)/(s'*H*s - s'*y);
        end
        
        r = theta*y + (1 - theta)*H*s;
        
        H = H + ((s'*r + r'*H*r)*(s*s')/((s'*r)^2)) - ((H*r*s' + s*r'*H)/(s'*r));
    else
        break
    end
end

end
