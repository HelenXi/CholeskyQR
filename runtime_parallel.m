n = 100;
l = 2*n;
v = ones(n,1);
sigma = 10^(-5);
for k = 2:n
    v(k) = sigma^((k-1)/(n-1));
end

num_workers = parpool('Processes');
if num_workers == 0
   parpool('open', 'local', num_workers);
end

trial = 10;
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
        
        num3 = tic();
        [Q,R] = choleskyQR2_parallel(A);
        cholesky2(i) = cholesky2(i) + toc(num3)/trial;
        num4 = tic();
        [Q,R] = sCholeskyQR3_parallel(A);
        cholesky3(i) = cholesky3(i) + + toc(num4)/trial;

        num6 = tic();
        [Q,R] = rQR_CholeskyQR_parallel(A);
        qrcho(i) = qrcho(i) + toc(num6)/trial;
    end
end

delete(gcp("nocreate"));

x = 10^4 + 165000.*(0:6);
figure;
plot(x, cholesky2, '+-', 'Color', '#00008B', 'LineWidth', 2);
hold on;
plot(x, cholesky3, '+-', 'Color', '#FF8C00', 'LineWidth', 2);
plot(x, qrcho, '*-','Color', '#008B8B', 'LineWidth', 2);

xlabel('number of rows');
ylabel('running time');
title('Runtime of Different Algorithms');
legend('choleskyQR2', 'scholeskyQR3','rQRCholeskyQR');
grid on;
hold off;




 