trial = 20;

acr1 = zeros(21,1);
rt1 = zeros(21,1);
acr2 = zeros(21,1);
rt2 = zeros(21,1);
acr3 = zeros(21,1);
rt3 = zeros(21,1);
acr4 = zeros(21,1);
rt4 = zeros(21,1);

num_workers = parpool('Processes');
if num_workers == 0
   parpool('open', 'local', num_workers);
end

m = 10^5;
n = 100;
v = ones(n,1);
sigma = 10^(-5);
for k = 2:n
    v(k) = sigma^((k-1)/(n-1));
end

for i = 1:21
    l = 100 + 10*(i-1);
    for j = 1:trial 
        Y = rand([m,n]);
        [U,S,V] = svd(Y,'econ');
        A = U*diag(v)*transpose(V);
        Af = norm(A,'fro');
        num = tic();
        [Q,R] = rQR_CholeskyQR_parallel_l(A,l);
        rt1(i) = rt1(i) + toc(num);
        acr1(i) = acr1(i) + norm(Q*R - A,'fro')/Af/trial;
    end
end

m = 10^4;
n = 100;
v = ones(n,1);
sigma = 10^(-5);
for k = 2:n
    v(k) = sigma^((k-1)/(n-1));
end

for i = 1:21
    l = 100 + 10*(i-1);
    for j = 1:trial 
        Y = rand([m,n]);
        [U,S,V] = svd(Y,'econ');
        A = U*diag(v)*transpose(V);
        Af = norm(A,'fro');
        num = tic();
        [Q,R] = rQR_CholeskyQR_parallel_l(A,l);
        rt2(i) = rt2(i) + toc(num);
        acr2(i) = acr2(i) + norm(Q*R - A,'fro')/Af/trial;
    end
end

m = 10^4;
n = 300;
v = ones(n,1);
sigma = 10^(-5);
for k = 2:n
    v(k) = sigma^((k-1)/(n-1));
end

for i = 1:21
    l = 300 + 20*(i-1);
    for j = 1:trial 
        Y = rand([m,n]);
        [U,S,V] = svd(Y,'econ');
        A = U*diag(v)*transpose(V);
        Af = norm(A,'fro');
        num = tic();
        [Q,R] = rQR_CholeskyQR_parallel_l(A,l);
        rt3(i) = rt3(i) + toc(num);
        acr3(i) = acr3(i) + norm(Q*R - A,'fro')/Af/trial;
    end
end

m = 10^4;
n = 100;
v = ones(n,1);
sigma = 10^(-13);
for k = 2:n
    v(k) = sigma^((k-1)/(n-1));
end

for i = 1:21
    l = 100 + 10*(i-1);
    for j = 1:trial 
        Y = rand([m,n]);
        [U,S,V] = svd(Y,'econ');
        A = U*diag(v)*transpose(V);
        Af = norm(A,'fro');
        num = tic();
        [Q,R] = rQR_CholeskyQR_parallel_l(A,l);
        rt4(i) = rt4(i) + toc(num);
        acr4(i) = acr4(i) + norm(Q*R - A,'fro')/Af/trial;
    end
end

delete(gcp("nocreate"));

x = 1 + (0.1).*(0:20);
figure;
plot(x, acr1, 'o-','Color', '#8B0000', 'LineWidth', 2);
hold on;
plot(x, acr2, '+-','Color', '#006400', 'LineWidth', 2);
plot(x, acr3, '+-', 'Color', '#00008B', 'LineWidth', 2);
plot(x, acr4, '+-', 'Color', '#FF8C00', 'LineWidth', 2);

xlabel('sampling rate');
ylabel('relative error');
title('accuracy of different sampling rates');
legend('m = 10^5,n = 100,sigma = 10^5$', 'm = 10^4,n = 100,sigma = 10^5', 'm = 10^4,n = 300,sigma = 10^5', 'm = 10^4,n = 100,sigma = 10^13');
grid on;
hold off;

figure;
plot(x, rt1, 'o-','Color', '#8B0000', 'LineWidth', 2);
hold on;
plot(x, rt2, '+-','Color', '#006400', 'LineWidth', 2);
plot(x, rt3, '+-', 'Color', '#00008B', 'LineWidth', 2);
plot(x, rt4, '+-', 'Color', '#FF8C00', 'LineWidth', 2);

xlabel('sampling rate');
ylabel('running time');
title('runtime of different sampling rates');
legend('m = 10^5,n = 100,sigma = 10^5$', 'm = 10^4,n = 100,sigma = 10^5', 'm = 10^4,n = 300,sigma = 10^5', 'm = 10^4,n = 100,sigma = 10^13');
grid on;
hold off;


 