n = 100;
l = 2*n;
v = ones(n,1);
sigma = 10^(-5);
for k = 2:n
    v(k) = sigma^((k-1)/(n-1));
end

trial = 6;
householder = zeros(7,1);
cholesky = zeros(7,1);
cholesky2 = zeros(7,1);
cholesky3 = zeros(7,1);
lucho = zeros(7,1);
qrcho = zeros(7,1);

for i = 1:7
    m = 10^4 + 165000*(i-1);
    for j = 1:trial 
        Y = rand([m,n]);
        [U,S,V] = svd(Y,'econ');
        A = U*diag(v)*transpose(V);
        Af = norm(A,'fro');
        [Q,R] = qr(A,"econ");
        householder(i) = householder(i) + norm(Q*R - A,'fro')/Af/trial;
        [Q,R] = choleskyQR(A);
        cholesky(i) = cholesky(i) + norm(Q*R - A,'fro')/Af/trial;
        [Q,R] = choleskyQR2(A);
        cholesky2(i) = cholesky2(i) + norm(Q*R - A,'fro')/Af/trial;
        [Q,R] = sCholeskyQR3(A);
        cholesky3(i) = cholesky3(i) + norm(Q*R - A,'fro')/Af/trial;
        [Q,R] = rLU_CholeskyQR(A);
        lucho(i) = lucho(i) + norm(Q*R - A,'fro')/Af/trial;
        [Q,R] = rQR_CholeskyQR(A);
        qrcho(i) = qrcho(i) + norm(Q*R - A,'fro')/Af/trial;
    end
end

x = 10^4 + 165000.*(0:6);
figure;
plot(x, householder, 'o-','Color', '#8B0000', 'LineWidth', 2);
hold on;
plot(x, cholesky, '+-','Color', '#006400', 'LineWidth', 2);
plot(x, cholesky2, '+-', 'Color', '#00008B', 'LineWidth', 2);
plot(x, cholesky3, '+-', 'Color', '#FF8C00', 'LineWidth', 2);
plot(x, lucho, '*-', 'Color', '#800080', 'LineWidth', 2);
plot(x, qrcho, '*-','Color', '#008B8B', 'LineWidth', 2);

xlabel('number of rows');
ylabel('running time');
title('Stability of Different Algorithms');
legend('householderQR', 'choleskyQR', 'choleskyQR2', 'scholeskyQR3', 'rLUCholeskyQR', 'rQRCholeskyQR');
grid on;
hold off;


 