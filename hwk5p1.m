function [x] = hwk5p1(x_in);

n = length(x_in)^2;
nlin = 0;
nfunc = 0;

epsilon = 1e-04;
iter_max = 200;
x = x_in;

amazing = 0.9;
sufficient = 0.1;
wider = 2;
narrower = 0.5;
delta_cap = n^2;

delta = NaN;

for k = 1:iter_max
    if k == 1
        [f,g,H] = rosenbrocknfgH(x);
        nfunc = nfunc + 1;
    end
    if isnan(delta)
        delta = 0.4*norm(g);
    end
    [~,flag] = chol(H);
    
    if flag ~= 0
        p = -(g/norm(g)) * delta;
    else
        
        p_c = (-norm(g)^2/(g'*H*g)) * g;
        p_n = -H\g;
        nlin = nlin + 1;
        p_d = p_c - p_n;
        
        if delta <= norm(p_c)
            p = -delta*p_c/(norm(p_c));
        elseif delta > norm(p_n)
            p = p_n;
        else
            a = p_c'*p_c;
            b = p_n'*p_n;
            c = p_d'*p_d;
            d = (a+b-c)/2;
            tau = (b - delta^2)/(b - d + sqrt(d^2 - a*b + (delta^2)*c));
            p = tau*p_c + (1-tau)*p_n;
        end
       
    end
    
    
    [new_f,new_g,new_H] = rosenbrocknfgH(x + p);
    nfunc = nfunc + 1;
    rho = (f - new_f)/(-g'*p - 0.5*p'*H*p);
    
    if rho >= amazing
        x = x + p;
        delta = wider*delta;
        f = new_f;
        g = new_g;
        H = new_H;
    elseif rho >= sufficient
        x = x + p;
        f = new_f;
        g = new_g;
        H = new_H;
    else
        delta = narrower*delta;
    end
    
    if norm(g) <= epsilon && flag == 0
        break
    end
    
end




