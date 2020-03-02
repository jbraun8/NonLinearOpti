function [ts,td,Ls,Ld] = hwk4p1(n)


ts = zeros(1,5);
td = zeros(1,5);

L_sparse = sparse([1:n,2:n,1:n-1],[1:n,1:n-1,2:n],[2*ones(1,n),-1*ones(1,2*n-2)],n,n);
L_dense = full(L_sparse);

for i = 1:5
    tic
    Ls = chol(L_sparse);
    ts(1,i) = toc;
    
    tic
    Ld = chol(L_dense);
    td(1,i) = toc;
end


    
    


