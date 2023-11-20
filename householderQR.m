function [Q,R] = householderQR(A)
    H = @(u,x) x - u*(u'*x);
    [m,n] = size(A);
    R = zeros(n,n);
    Q = A;
    for j = 1:min(m,n)
        u = house_gen(Q(j:m,j));
        R(j:m,j) = u;
        Q(j:m,j:n) = H(u,Q(j:m,j:n));
        Q(j+1:m,j) = 0;
    end
end