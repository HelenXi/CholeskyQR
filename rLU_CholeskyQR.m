function [Q,R] = rLU_CholeskyQR(A)
    [m,n] = size(A);
    l = n*2;
    A1 = A(randi(m,1,l),:);
    [L,U] = lu(A1);
    X = A*inv(U);
    [Q,R2] = choleskyQR(X);
    R = R2*U;
end