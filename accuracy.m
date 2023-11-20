m = 10^5;
n = 100;
l = 2*n;
v = ones(n,1);


tnum = 15;
trial = 20;
householder = zeros(tnum,1);
cholesky = zeros(tnum,1);
cholesky2 = zeros(tnum,1);
cholesky3 = zeros(tnum,1);
lucho = zeros(tnum,1);
qrcho = zeros(tnum,1);

for j = 1:trial
    Y = rand([m,n]);
    [U,S,V] = svd(Y,'econ');
    for i = 1:8
        sigma = 10^(-i);
        for k = 2:n
            v(k) = sigma^((k-1)/(n-1));
        end
        A = U*diag(v)*transpose(V);
        Af = norm(A,'fro');
        [Q,R] = qr(A,"econ");
        householder(i) = householder(i) + norm(Q*R - A,'fro')/Af;
        [Q,R] = choleskyQR(A);
        cholesky(i) = cholesky(i) + norm(Q*R - A,'fro')/Af;
        [Q,R] = choleskyQR2(A);
        cholesky2(i) = cholesky2(i) + norm(Q*R - A,'fro')/Af;
        [Q,R] = sCholeskyQR3(A);
        cholesky3(i) = cholesky3(i) + norm(Q*R - A,'fro')/Af;
        [Q,R] = rLU_CholeskyQR(A);
        lucho(i) = lucho(i) + norm(Q*R - A,'fro')/Af;
        [Q,R] = rQR_CholeskyQR(A);
        qrcho(i) = qrcho(i) + norm(Q*R - A,'fro')/Af;
    end
    for i = 9:tnum
        sigma = 10^(-i);
        for k = 2:n
            v(k) = sigma^((k-1)/(n-1));
        end
        A = U*diag(v)*transpose(V);
        Af = norm(A,'fro');
        [Q,R] = qr(A,"econ");
        householder(i) = householder(i) + norm(Q*R - A,'fro')/Af; 
        [Q,R] = sCholeskyQR3(A);
        cholesky3(i) = cholesky3(i) + norm(Q*R - A,'fro')/Af;
        [Q,R] = rLU_CholeskyQR(A);
        lucho(i) = lucho(i) + norm(Q*R - A,'fro')/Af;
        [Q,R] = rQR_CholeskyQR(A);
        qrcho(i) = qrcho(i) + norm(Q*R - A,'fro')/Af;
    end
end

householder = householder/trial;
cholesky = cholesky/trial;
cholesky2 = cholesky2/trial;
cholesky3 = cholesky3/trial;
lucho = lucho/trial;
qrcho = qrcho/trial;

x = (1:15);
figure;
plot(x, householder, 'o-','Color', '#8B0000', 'LineWidth', 2);
hold on;
plot(x(1:8), cholesky(1:8), '+-','Color', '#006400', 'LineWidth', 2);
plot(x(1:8), cholesky2(1:8), '+-', 'Color', '#00008B', 'LineWidth', 2);
plot(x, cholesky3, '+-', 'Color', '#FF8C00', 'LineWidth', 2);
plot(x, lucho, '*-', 'Color', '#800080', 'LineWidth', 2);
plot(x, qrcho, '*-','Color', '#008B8B', 'LineWidth', 2);

xlabel('conditional number (10^i)');
ylabel('relative error');
title('Accuray of Different Algorithms');
legend('householderQR', 'choleskyQR', 'choleskyQR2', 'scholeskyQR3', 'rLUCholeskyQR', 'rQRCholeskyQR');
grid on;
hold off;


 