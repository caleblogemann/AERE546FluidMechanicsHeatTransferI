%% Problem 1
n = 121;
deltaX = 1/n;
x = linspace(deltaX, 1-deltaX, n-1);
tFinal = 0.3;
alpha = 0.5;
deltaT = alpha*deltaX^2;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;

% initial conditions
T0Func = @(x) zeros(size(x));
T0 = T0Func(x);

% Boundary condtions
Tl = 0;
Tr = 1;
RHSFunc = @(t, T) ([T(2:end),Tr]-2*T+[Tl,T(1:end-1)])/deltaX^2;

sol = forwardEuler(RHSFunc, T0, nTimeSteps, deltaT);

plotInterval = 0.04;
nPlots = ceil(tFinal/plotInterval);
timeStepInterval = floor(nTimeSteps/nPlots);
hold on
for n = 1:timeStepInterval:nTimeSteps
    plot([0,x,1], [0,sol(n,:),1], 'k', 'LineWidth', 2);
end
xlim([-0.01,1.01]);
ylim([-0.01,1.01]);
xlabel('x', 'FontSize', 16);
ylabel('T', 'FontSize', 16);
hold off
saveas(gcf, 'Figures/03_01.png', 'png');

tStart = 0.01;
nInitial = ceil(tStart/deltaT);
ht = (ones(nTimeSteps - nInitial+2,1) - sol(nInitial:end,end))/deltaX;
t = linspace(0,tFinal, nTimeSteps - nInitial+2);
plot(t,ht, 'k', 'LineWidth', 2);
xlabel('t', 'FontSize', 16);
ylabel('h_t', 'FontSize', 16);
title('Bulk Heat Transfer Coefficient', 'FontSize', 16);
saveas(gcf, 'Figures/03_02.png', 'png');

%% Problem 2
n = 250;
tFinal = 0.05;
a = 0;
b = 1;
deltaX = (b - a)/n;
x = linspace(a, b, n+1);
deltaT = 0.1*deltaX^2;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;
kappa = @(x) 0.1 + 2.0*exp(-5*x);
dkappa = @(x) -10.0*exp(-5*x);

% Initial conditions
sigma = 0.1;
T0Func = @(x) exp(-(x - 0.5).^2/sigma^2)/(sigma*sqrt(pi));
T0 = T0Func(x);

% Boundary Conditions
% zero flux
RHSFunc = @(t, T) dkappa(x).*(([T(2:end), T(end)] - [T(1), T(1:end-1)])/(2*deltaX)) ...
    + kappa(x).*(([T(2:end),T(end)]-2*T+[T(1),T(1:end-1)])/deltaX^2);

sol = RK2(RHSFunc, T0, nTimeSteps, deltaT);

plotInterval = 0.01;
nPlots = ceil(tFinal/plotInterval);
timeStepInterval = floor(nTimeSteps/nPlots);
hold on
for n = 1:timeStepInterval:nTimeSteps
    plot(x, sol(n,:), 'k', 'LineWidth', 2);
end
xlim([-0.01,1.01]);
%ylim([-0.01,1.01]);
xlabel('x', 'FontSize', 16);
ylabel('T', 'FontSize', 16);
hold off
saveas(gcf, 'Figures/03_03.png', 'png');

intT = sum(deltaX*sol, 2);
t = linspace(0, tFinal, nTimeSteps+1);
plot(t, intT, 'k', 'LineWidth', 2);
xlabel('t', 'FontSize', 16);
ylabel('Integral of T', 'FontSize', 16);
saveas(gcf, 'Figures/03_04.png', 'png');

%% Problem 3
n = 121;
deltaX = 1/n;
x = linspace(0, 1-deltaX, n);
tFinal = 0.4;
alpha = 0.5;
deltaT = alpha*deltaX^2;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;

% initial conditions
T0Func = @(x) zeros(size(x));
T0 = T0Func(x);

% Boundary condtions
% zero flux at x = 0
Tr = 1;
RHSFunc = @(t, T) ([T(2:end),Tr]-2*T+[T(1),T(1:end-1)])/deltaX^2;

sol = crankNicholson(T0, nTimeSteps, deltaT, deltaX);

plotInterval = 0.04;
nPlots = ceil(tFinal/plotInterval);
timeStepInterval = floor(nTimeSteps/nPlots);
hold on
for n = 1:timeStepInterval:nTimeSteps
    plot([x,1], [sol(n,:),1], 'k', 'LineWidth', 2);
end
xlim([-0.01,1.01]);
ylim([-0.01,1.01]);
xlabel('x', 'FontSize', 16);
ylabel('T', 'FontSize', 16);
hold off
saveas(gcf, 'Figures/03_05.png', 'png');

tStart = 0.01;
nInitial = ceil(tStart/deltaT);
ht = (ones(nTimeSteps+1 - nInitial+1,1) - sol(nInitial:end,end))/deltaX;
t = linspace(0,tFinal, nTimeSteps+1 - nInitial+1);
plot(t,ht, 'k', 'LineWidth', 2);
xlabel('t', 'FontSize', 16);
ylabel('h_t', 'FontSize', 16);
title('Bulk Heat Transfer Coefficient', 'FontSize', 16);
saveas(gcf, 'Figures/03_02.png', 'png');
