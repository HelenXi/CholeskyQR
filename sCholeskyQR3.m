function [Q,R] = sCholeskyQR3(A)
    [m,n] = size(A);
    G = transpose(A)*A;
    R1 = chol(G + 10^(-15)*eye(n));
    Q1 = A*inv(R1);
    [Q2,R2] = choleskyQR(Q1);
    [Q,R3] = choleskyQR(Q2);
    R = R3*R2*R1;
end