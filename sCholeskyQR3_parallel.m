function [Q,R] = sCholeskyQR3_parallel(A)
    [m,n] = size(A);
    G = transpose(A)*A;
    R1 = chol(G + 10^(-15)*eye(n));
    Q1 = A*multi_process_inverse(R1);
    [Q2,R2] = choleskyQR_parallel(Q1);
    [Q,R3] = choleskyQR_parallel(Q2);
    R = R3*R2*R1;
end