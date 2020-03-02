function [f,g,H] = rosenbrocknfgH(x)
    n = length(x);
    f = 0;
    g = zeros(n,1);
    H = spalloc(n,n,3*n-2);
    for j = 1:n
        if j == 1
            f = f + 100*(x(j+1) - x(j)^2)^2 + (x(j) - 1)^2;
            grad = 400*x(j)^3 - 400*x(j)*x(j+1) + 2*x(j) - 2;
            H(1:2,1) = [1200*x(j)^2 - 400*x(j+1) + 2;-400*x(j)];
        elseif j == n
            grad = 200*x(j) - 200*x(j-1)^2;
            H(n-1:n,n) = [-400*x(j-1);200];
        else
            f = f + 100*(x(j+1) - x(j)^2)^2 + (x(j) - 1)^2;
            grad = 400*x(j)^3 - 400*x(j+1)*x(j) + 202*x(j) - 200*x(j-1)^2 -2;
            H(j-1:j+1,j) = [-400*x(j-1);1200*x(j)^2 - 400*x(j+1) + 202;-400*x(j)];
        end
        g(j) = grad;
    end
end