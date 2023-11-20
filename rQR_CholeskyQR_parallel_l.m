function [Q,R] = rQR_CholeskyQR_parallel_l(A,l)
    [m,n] = size(A);
    A1 = A(randi(m,1,l),:);
    [~,R1] = qr(A1,'econ');
    X = A*multi_process_inverse(R1);
    [Q,R2] = choleskyQR_parallel(X);
    R = R2*R1;
end