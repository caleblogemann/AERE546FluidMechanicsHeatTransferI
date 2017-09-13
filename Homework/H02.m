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
plot(x, sol(1,:), x, sol(2,:), x, sol(3,:));
legend('M = 1', 'M = 5', 'M = 9');

%% Problem 2b



%% Problem 2c

