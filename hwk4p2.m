function [xf,nfunc,nlin] = hwk4p2(n,iter,epsilon)

x0 = 10*ones(n,1);

delta = 10e-05;
nfunc = 0;
nlin = 0;

path = [];
norm_grad = Inf;
x = x0;

for i = 1:iter
    if norm_grad > epsilon
        path = [path;x'];
        [f,G,B] = rosenbrocknfgH(x);
        nfunc = nfunc + 1;
        norm_grad = norm(G);
        
        % LDL subroutine
        B = full(B);
        [L,D,P] = ldl(B);
        nlin = nlin + 3;
        [Q,V] = eig(D);
        evals = diag(V);
        tau = delta*ones(length(evals),1);
        perturb = tau > evals;
        F = Q * diag(delta*perturb - (perturb.*evals)) * Q';
        E = P' * L * F * L' * P;
        B = B + E;
        %
     
        p = -B\G;
        nlin = nlin + 1;
    
        alpha = 1;
        c = 0.5;
        rho = 0.75;
    
        [step_val,~,~] = rosenbrocknfgH(x + alpha*p);
        nfunc = nfunc + 1;
        requirement = f + c*alpha*G'*p;
        select_alpha = requirement - step_val;
    
        while select_alpha < 0
            alpha = rho*alpha;
            [step_val,~,~] = rosenbrocknfgH(x + alpha*p);
            nfunc = nfunc + 1;
            requirement = f + c*alpha*G'*p;
            select_alpha = requirement - step_val;
        end
    
        x = x + alpha * p;
    end
end

xf = x;