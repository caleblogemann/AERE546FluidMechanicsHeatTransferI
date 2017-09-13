%% Problem 2a
n = 100;
deltaX = 1/n;
nM = 3;
sol = zeros(nM, n+1);
sol(:,1) = ones(nM,1)';
for iter = [1, 5, 9; 1:nM]
    M = iter(1);
    i = iter(2);
    mainDiagonal  = (-M*deltaX^2 - 2)*ones(n-1,1);
    lowerDiagonal = ones(n-2,1);
    upperDiagonal = ones(n-2,1);
    RHS = zeros(n-1,1);
    % boundary conditions
    RHS(1) = -1;
    sol(i,2:end-1) = tridiag(n-1, mainDiagonal, lowerDiagonal, upperDiagonal, RHS);
end

x = linspace(0, 1, n+1);
plot(x, sol(1,:), 'k+', x, sol(2,:), 'k--', x, sol(3,:), 'ko');
legend('M = 1', 'M = 5', 'M = 9');
title('Dirichlet Boundary Conditions');
xlabel('x');
ylabel('T');
saveas(gcf, 'Figures/02_01.png', 'png');

% -dT/dx at x = 1 can be approximated by -dT/dx = (T_(n-1) - T_n)/deltaX
% T_n = sol(i, end) = sol(i, n+1)
heatFlux = zeros(nM, 1);
for i = 1:nM
    heatFlux(i) = (sol(i, end-1) - sol(i,end))/deltaX;
end
disp(heatFlux);

%% Problem 2b
n = 100;
deltaX = 1/n;
nM = 3;
sol2 = zeros(nM, n+1);
sol2(:,1) = ones(nM,1)';
for iter = [1, 5, 9; 1:nM]
    M = iter(1);
    i = iter(2);
    mainDiagonal  = (-M*deltaX^2 - 2)*ones(n,1);
    mainDiagonal(end) = (-M*deltaX^2 - 1);
    lowerDiagonal = ones(n-1,1);
    upperDiagonal = ones(n-1,1);
    RHS = zeros(n,1);
    % boundary conditions
    RHS(1) = -1;
    sol2(i,2:end) = tridiag(n, mainDiagonal, lowerDiagonal, upperDiagonal, RHS);
end

exactSol1 = @(x) (exp(-3)/(exp(-3) - exp(3))) * exp(3*x) + (-exp(3)/(exp(-3) - exp(3))) * exp(-3*x);
exactSol2 = @(x) (exp(-3)/(exp(-3) + exp(3))) * exp(3*x) + (exp(3)/(exp(-3) + exp(3))) * exp(-3*x);
x = linspace(0, 1, n+1);
plot(x, sol(3,:), 'k--', x, exactSol1(x), 'k-');
legend('Numerical Solution', 'Exact Solution ');
title('Cold Reservoir');
xlabel('x');
ylabel('T');
saveas(gcf, 'Figures/02_02.png', 'png');

plot(x, sol2(3,:), 'k--', x, exactSol2(x), 'k-');
legend('Numerical Solution', 'Exact Solution');
title('Insulated');
xlabel('x');
ylabel('T');
saveas(gcf, 'Figures/02_03.png', 'png');

disp(sol2(:,end))

%% Problem 2c
n = 100;
deltaX = 1/n;
sol3 = zeros(1, n+1);
sol3(1,1) = 1;
M = 9;
mainDiagonal  = (-M*deltaX^2 - 2)*ones(n,1);
mainDiagonal(end) = (-M*deltaX^2 - 1);
lowerDiagonal = ones(n-1,1);
upperDiagonal = ones(n-1,1);
x = linspace(deltaX, 1, n);
RHS = (-100*deltaX^2)*((x.^2).*((1-x).^2));
% boundary conditions
RHS(1) = RHS(1) - 1;
sol3(1,2:end) = tridiag(n, mainDiagonal, lowerDiagonal, upperDiagonal, RHS);

x = linspace(0,1,n+1);
plot(x, sol3(1,:), 'k-');
title('Heat Source');
xlabel('x');
ylabel('T');
saveas(gcf, 'Figures/02_04.png', 'png');
