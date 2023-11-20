function [Q,R] = choleskyQR2_parallel(A)
    [Q1,R1] = choleskyQR_parallel(A);
    [Q,R2] = choleskyQR_parallel(Q1);
    R = R2*R1;
end