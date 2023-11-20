function [Q,R] = choleskyQR(A)
    G = transpose(A)*A;
    R = chol(G);
    Q = A*inv(R);
end