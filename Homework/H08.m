%% Problem 1
a = -1;
b = 1;

f = @(u) u;
tFinal = 2.0;
deltaX = 0.01;
nCells = (b - a)/deltaX;
x = linspace(a + deltaX/2, b - deltaX/2, nCells);
u0func = @(x) 0.5 + 0.5*(1 - x.^2).^6;
u0 = u0func(x);
exactSol = @(x, t) u0func(mod((x - t) + 1, 2) - 1);

t = [0.0, 0.5, 1.0, 1.5, 2.0];

% 1. 1.a
cfl = 0.3;
deltaT = cfl*deltaX;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;
p = .5;

tIndices = floor(t/deltaT)+1;

sol = upwindDownwind(u0, f, p, deltaX, deltaT, nTimeSteps);
plot(x, sol(tIndices(1:end-1), :), 'k');
saveas(gcf, 'Figures/08_01.png', 'png');
plot(x, sol(tIndices(end), :), 'k');
saveas(gcf, 'Figures/08_02.png', 'png');

% 1.2
p = 1.0;
cfl = 0.5;
deltaT = cfl*deltaX;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;
tIndices = floor(t/deltaT)+1;

sol = upwindDownwind(u0, f, p, deltaX, deltaT, nTimeSteps);
plot(x, sol(tIndices, :), 'k', x, exactSol(x, 2.0), 'o');
saveas(gcf, 'Figures/08_03.png', 'png');

%1.3
p = 0.55;
cfl = 0.1;
deltaT = cfl*deltaX;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;
tIndices = floor(t/deltaT)+1;
sol = upwindDownwind(u0, f, p, deltaX, deltaT, nTimeSteps);
plot(x, sol(tIndices, :), 'k');
saveas(gcf, 'Figures/08_04.png', 'png');

%% 2
a = -1;
b = 1;

f = @(u) 0.5*u^2;
df = @(u) u;
tFinal = 1.2;
nCells = 201;
deltaX = (b - a)/(nCells-1);
x = linspace(a + deltaX/2, b - deltaX/2, nCells);
u0func = @(x) 0.5 + 0.5*(1 - x.^2).^6;
u0 = u0func(x);
t = 0:0.2:tFinal;

% RK2 central
iter = 0;
RHSFunc = @(t, u) central(u, f, deltaX);
for cfl = [0.2, 1.0]
    iter = iter+1;
    deltaT = cfl*deltaX;
    nTimeSteps = ceil(tFinal/deltaT);
    deltaT = tFinal/nTimeSteps;
    tIndices = floor(t/deltaT)+1;
    sol = RK2(RHSFunc, u0, nTimeSteps, deltaT);
    plot(x, sol(tIndices, :), 'k');
    xlabel('x');
    ylabel('u');
    title(strcat('RK2 Central cfl = ', num2str(cfl)));
    saveas(gcf, strcat('Figures/08_0',num2str(4 + iter),'.png'), 'png');
end

% RK2 upwind
iter = 0;
RHSFunc = @(t, u) upwind(u, f, deltaX);
for cfl = [0.2, 1.0]
    iter = iter+1;
    deltaT = cfl*deltaX;
    nTimeSteps = ceil(tFinal/deltaT);
    deltaT = tFinal/nTimeSteps;
    tIndices = floor(t/deltaT)+1;
    sol = RK2(RHSFunc, u0, nTimeSteps, deltaT);
    plot(x, sol(tIndices, :), 'k');
    xlabel('x');
    ylabel('u');
    title(strcat('RK2 Upwind cfl = ', num2str(cfl)));
    saveas(gcf, strcat('Figures/08_0',num2str(6 + iter),'.png'), 'png');
end
% Lax Wendroff
iter = 0;
for cfl = [0.2, 1.0]
    iter = iter+1;
    deltaT = cfl*deltaX;
    nTimeSteps = ceil(tFinal/deltaT);
    deltaT = tFinal/nTimeSteps;
    tIndices = floor(t/deltaT)+1;
    sol = laxWendroff(f, df, u0, deltaT, deltaX, nTimeSteps);
    plot(x, sol(tIndices, :), 'k');
    xlabel('x');
    ylabel('u');
    title(strcat('Lax Wendroff cfl = ', num2str(cfl)));
    saveas(gcf, strcat('Figures/08_0',num2str(8 + iter),'.png'), 'png');
end

% Maccormack
iter = 0;
for cfl = [0.2, 1.0]
    iter = iter+1;
    deltaT = cfl*deltaX;
    nTimeSteps = ceil(tFinal/deltaT);
    deltaT = tFinal/nTimeSteps;
    tIndices = floor(t/deltaT)+1;
    sol = maccormack(f, u0, deltaT, deltaX, nTimeSteps);
    plot(x, sol(tIndices, :), 'k');
    xlabel('x');
    ylabel('u');
    title(strcat('MacCormick cfl = ', num2str(cfl)));
    saveas(gcf, strcat('Figures/08_1',num2str(iter),'.png'), 'png');
end

%% Problem 3
a = -1;
b = 1;

f = @(u) 0.5*u^2;
df = @(u) u;
tFinal = 1.8;
nCells = 201;
deltaX = (b - a)/(nCells-1);
x = linspace(a + deltaX/2, b - deltaX/2, nCells);
u0func = @(x) abs(x) <= 1/3;
u0 = u0func(x);
t = 0:0.3:tFinal;
cfl = 0.2;
deltaT = cfl*deltaX;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;
tIndices = floor(t/deltaT)+1;

RHSFunc = @(t, u) upwind(u, f, deltaX);
sol1 = RK2(RHSFunc, u0, nTimeSteps, deltaT);
RHSFunc = @(t, u) upwindConvection(u, deltaX);
sol2 = RK2(RHSFunc, u0, nTimeSteps, deltaT);
sol3 = laxWendroff(f, df, u0, deltaT, deltaX, nTimeSteps);

plot(x, sol1(tIndices, :), 'k');
xlabel('x');
ylabel('y');
title('RK2 Upwind Conservation Form');
saveas(gcf, 'Figures/08_13.png', 'png');

plot(x, sol2(tIndices, :), 'k');
xlabel('x');
ylabel('y');
title('RK2 Upwind Convection Form');
saveas(gcf, 'Figures/08_14.png', 'png');

plot(x, sol3(tIndices, :), 'k');
xlabel('x');
ylabel('y');
title('Lax-Wendroff');
saveas(gcf, 'Figures/08_15.png', 'png');
