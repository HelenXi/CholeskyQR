function [Q,R] = choleskyQR2(A)
    [Q1,R1] = choleskyQR(A);
    [Q,R2] = choleskyQR(Q1);
    R = R2*R1;
end