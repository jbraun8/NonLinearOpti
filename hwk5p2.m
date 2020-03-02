function [time1,err1,time2,err2,time3,err3] = hwk5p2(n)
    
x = 10*ones(n,1);
[~,g,H] = rosenbrock1(x);
[~,~,M] = rosenbrocknfgH(x);

tol = 1e-14;

tic;
R = chol(H);
y = R\(R'\g);
time1 = toc;
err1 = norm(H*y - g);

tic;
[y,~,err2,iter2] = pcg(H,g,tol);
time2 = toc;

tic;
[y,~,err3,iter3] = pcg(H,g,tol,[],M);
time3 = toc;

end