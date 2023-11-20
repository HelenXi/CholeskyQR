function [Q,R] = choleskyQR_parallel(A)
    G = transpose(A)*A;
    R = chol(G);
    Q = A*multi_process_inverse(R);
end