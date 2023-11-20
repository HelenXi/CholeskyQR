function [Q,R] = rQR_CholeskyQR(A)
    [m,n] = size(A);
    l = n*2;
    A1 = A(randi(m,1,l),:);
    [~,R1] = qr(A1,'econ');
    X = A*inv(R1);
    [Q,R2] = choleskyQR(X);
    R = R2*R1;
end