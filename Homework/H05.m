%% problem 1
a = 0;
b = 1;
c = 0;
d = 1;
iter = 0;
lambda = 1.9;
for n = [1,3]
    for N = [15, 151]
        iter = iter +1;
        deltaX = (b - a)/(N + 2);
        deltaY = (d - c)/(N + 2);
        x = linspace(a+deltaX, b-deltaX, N);
        y = linspace(c+deltaY, d-deltaY, N);
        u0 = zeros(N);
        [u, k, res] = sor(lambda, u0, x, y, n, 1e-4, 2500);
        contour(x, y, u, 20);
        xlabel('x');
        ylabel('y');
        title(strcat('n=',num2str(n),' and N=',num2str(N)));
        saveas(gcf, strcat('Figures/05_',num2str(iter),'.png'), 'png')
    end
end

n = 1;
N = 151;
deltaX = (b - a)/(N + 2);
deltaY = (d - c)/(N + 2);
x = linspace(a+deltaX, b-deltaX, N);
y = linspace(c+deltaY, d-deltaY, N);
u0 = zeros(N);
figure;
hold on
for lambda = [1, 0.5, 1.5, 1.95]
    [u, k, res] = sor(lambda, u0, x, y, n, 1e-4, 2500);
    plot(2:k,log(res(2:k)./res(1)))
end
legend('lambda = 1', 'lambda = 0.5', 'lambda = 1.5', 'lambda = 1.95');
title('Residuals for relaxations');
xlabel('# of iterations');
ylabel('log(residual)');
saveas(gcf, 'Figures/05_5.png', 'png')
hold off

%% Problem 2


%% Problem 3
a = 0;
b = 1;
c = 0;
d = 1;
N = 151;
deltaX = (b - a)/(N + 2);
deltaY = (d - c)/(N + 2);
x = linspace(a+deltaX, b-deltaX, N);
y = linspace(c+deltaY, d-deltaY, N);

m = @(u) multByA(u, x, y, deltaX, deltaY);
rhsFunc = @(x, y) 2000*(sin(4*pi*x)*sin(4*pi*y))^5;

u0 = zeros(N);
rhs = zeros(N);
for i = 1:N
    for j = 1:N
        rhs(i, j) = rhsFunc(x(i), y(j));
    end
end
[u, iter, res] = conjugateGradient(m, u0, rhs, 1e-5, N^2);
contour(x, y, u');
xlabel('x');
ylabel('y');
title('Poisson Equation Solution by Conjugate Gradient');
saveas(gcf, 'Figures/05_8.png', 'png');

plot(1:iter,log(res/res(1)))
xlabel('# of Iterations');
ylabel('log(residual)');
title('Residual vs Iteration');
saveas(gcf, 'Figures/05_9.png', 'png');
