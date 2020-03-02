function [x] = hwk6p2(x_in,tol,iter)

n = length(x_in);
norm_grad = Inf;
nfunc = 0;
x = x_in;

c = 0.5;


for i = 0:iter
    if norm_grad > tol
        [f,g,H] = rosenbrocknfgH(x);
        nfunc = nfunc + 1;
        norm_grad = norm(g);
        eps = min(0.5,sqrt(norm_grad))*norm_grad;
        z = 0;
        r = g;
        d = -r;
        
        for j = 0:n
            if d'*H*d <= 0
                if j == 0
                    p = -g
                    break
                else
                    p = z;
                    break
                end
            end
            alpha_hat = (r'*r)/(d'*H*d);
            z = z + alpha_hat*d;
            rold = r;
            r = r + alpha_hat*H*d;
            if norm(r) < eps
                p = z;
                break
            end
            beta = (r'*r)/(rold'*rold);
            d = -r + beta*d;
        end
         
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
        
    else
        break
    end
    
end