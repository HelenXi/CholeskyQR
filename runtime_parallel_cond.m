m = 10^5;
n = 100;
l = 2*n;
v = ones(n,1);


tnum = 15;
trial = 20;

cholesky2 = zeros(tnum,1);
cholesky3 = zeros(tnum,1);
qrcho = zeros(tnum,1);

num_workers = parpool('Processes');
if num_workers == 0
   parpool('open', 'local', num_workers);
end

for j = 1:trial
    Y = rand([m,n]);
    [U,S,V] = svd(Y,'econ');
    for i = 1:8
        sigma = 10^(-i);
        for k = 2:n
            v(k) = sigma^((k-1)/(n-1));
        end
        A = U*diag(v)*transpose(V);
        num1 = tic();
        [Q,R] = choleskyQR2_parallel(A);
        cholesky2(i) = cholesky2(i) + toc(num1);
        num2 = tic();
        [Q,R] = sCholeskyQR3_parallel(A);
        cholesky3(i) = cholesky3(i) + toc(num2);
        num3 = tic()
        [Q,R] = rQR_CholeskyQR_parallel(A);
        qrcho(i) = qrcho(i) + toc(num3);
    end
    for i = 9:tnum
        sigma = 10^(-i);
        for k = 2:n
            v(k) = sigma^((k-1)/(n-1));
        end
        A = U*diag(v)*transpose(V);
        num2 = tic();
        [Q,R] = sCholeskyQR3_parallel(A);
        cholesky3(i) = cholesky3(i) + toc(num2);
        num3 = tic()
        [Q,R] = rQR_CholeskyQR_parallel(A);
        qrcho(i) = qrcho(i) + toc(num3);
    end 
end

delete(gcp("nocreate"));

cholesky2 = cholesky2/trial;
cholesky3 = cholesky3/trial;
qrcho = qrcho/trial;

x = (1:15);
figure;
plot(x(1:8), cholesky2(1:8), '+-', 'Color', '#00008B', 'LineWidth', 2);
hold on;
plot(x, cholesky3, '+-', 'Color', '#FF8C00', 'LineWidth', 2);
plot(x, qrcho, '*-','Color', '#008B8B', 'LineWidth', 2);

xlabel('conditional number (10^i)');
ylabel('running time');
title('Runtime of Different Algorithms');
legend('choleskyQR2', 'scholeskyQR3','rQRCholeskyQR');
grid on;
hold off;


 