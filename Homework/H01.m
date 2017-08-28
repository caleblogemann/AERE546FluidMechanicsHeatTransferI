% Problem 2a
tFinal = 32;
% initial conditions
x0 = [1, 0];

% damping rate
sigma = 0.0;

RHSFunc = @(t, x) [x(2), -sigma*x(2) - x(1)];

n1 = 21;
sol1 = RK2(RHSFunc, x0, n1, tFinal);
plot(linspace(0, 32, n1+1), sol1(:,1));
saveas(gcf, 'Figures/01_01.png', 'png');
n2 = 101;
sol2 = RK2(RHSFunc, x0, n2, tFinal);
plot(linspace(0, 32, n2+1), sol2(:,1));
saveas(gcf,'Figures/01_02.png', 'png');
n3 = 301;
sol3 = RK2(RHSFunc, x0, n3, tFinal);
plot(linspace(0, 32, n3+1), sol3(:,1));
saveas(gcf,'Figures/01_03.png', 'png');

exactSolFunc = @(t) cos(t);
exactSol = exactSolFunc(linspace(0, 32, 1000));
plot(linspace(0, 32, 1000), exactSol);
saveas(gcf,'Figures/01_04.png', 'png');

plot(linspace(0, 32, n2+1), sol2(:,1), 'ro',...
    linspace(0, 32, n3+1), sol3(:,1), 'b+',...
    linspace(0, 32, 1000), exactSol, 'k');
saveas(gcf,'Figures/01_05.png', 'png');

sigma = 0.5;

RHSFunc = @(t, x) [x(2), -sigma*x(2) - x(1)];

n1 = 21;
sol1 = RK2(RHSFunc, x0, n1, tFinal);
plot(linspace(0, 32, n1+1), sol1(:,1));
saveas(gcf,'Figures/01_06.png', 'png');
n2 = 101;
sol2 = RK2(RHSFunc, x0, n2, tFinal);
plot(linspace(0, 32, n2+1), sol2(:,1));
saveas(gcf,'Figures/01_07.png', 'png');
n3 = 301;
sol3 = RK2(RHSFunc, x0, n3, tFinal);
plot(linspace(0, 32, n3+1), sol3(:,1));
saveas(gcf,'Figures/01_08.png', 'png');

exactSolFunc = @(t) exp(-sigma/2*t).*cos(sqrt(15)/4 * t);
exactSol = exactSolFunc(linspace(0, 32, 1000));
plot(linspace(0, 32, 1000), exactSol);
saveas(gcf,'Figures/01_09.png', 'png');

plot(linspace(0, 32, n2+1), sol2(:,1), 'ro',...
    linspace(0, 32, n3+1), sol3(:,1), 'b+',...
    linspace(0, 32, 1000), exactSol, 'k');
saveas(gcf,'Figures/01_10.png', 'png');

% Problem 2b
tFinal = 32;
x0 = [1,0];
B = 0.2;
RHSFunc = @(t, x) [x(2), B*x(1)^3 - x(1)];
n = 1000;
sol = RK2(RHSFunc, x0, n, tFinal);
plot(linspace(0,32,n+1), sol(:,1));
title('Nonlinear Spring, B = 0.2')
saveas(gcf, 'Figures/01_11.png', 'png');

B = 0.6;
RHSFunc = @(t, x) [x(2), B*x(1)^3 - x(1)];
n = 1000;
sol = RK2(RHSFunc, x0, n, tFinal);
plot(linspace(0,32,n+1), sol(:,1));
title('Nonlinear Spring, B = 0.6')
saveas(gcf, 'Figures/01_12.png', 'png');

B = 0.9;
RHSFunc = @(t, x) [x(2), B*x(1)^3 - x(1)];
n = 1000;
sol = RK2(RHSFunc, x0, n, tFinal);
plot(linspace(0,32,n+1), sol(:,1));
title('Nonlinear Spring, B = 0.9')
saveas(gcf, 'Figures/01_13.png', 'png');

B = 0.999;
RHSFunc = @(t, x) [x(2), B*x(1)^3 - x(1)];
n = 3000;
sol = RK2(RHSFunc, x0, n, tFinal);
plot(linspace(0,32,n+1), sol(:,1));
title('Nonlinear Spring, B = 0.999')
saveas(gcf, 'Figures/01_14.png', 'png');

% Problem 3
tFinal = 32;
% initial conditions
x0 = [1, 0];

% damping rate
sigma = 0.0;

RHSFunc = @(t, x) [x(2), -sigma*x(2) - x(1)];

n1 = 21;
sol1 = AB2(RHSFunc, x0, n1, tFinal);
plot(linspace(0, 32, n1+1), sol1(:,1));
saveas(gcf, 'Figures/01_15.png', 'png');
n2 = 101;
sol2 = AB2(RHSFunc, x0, n2, tFinal);
plot(linspace(0, 32, n2+1), sol2(:,1));
saveas(gcf,'Figures/01_16.png', 'png');
n3 = 301;
sol3 = AB2(RHSFunc, x0, n3, tFinal);
plot(linspace(0, 32, n3+1), sol3(:,1));
saveas(gcf,'Figures/01_17.png', 'png');

exactSolFunc = @(t) cos(t);
exactSol = exactSolFunc(linspace(0, 32, 1000));
plot(linspace(0, 32, 1000), exactSol);
saveas(gcf,'Figures/01_18.png', 'png');

plot(linspace(0, 32, n2+1), sol2(:,1), 'ro',...
    linspace(0, 32, n3+1), sol3(:,1), 'b+',...
    linspace(0, 32, 1000), exactSol, 'k');
saveas(gcf,'Figures/01_19.png', 'png');

sigma = 0.5;

RHSFunc = @(t, x) [x(2), -sigma*x(2) - x(1)];

n1 = 21;
sol1 = AB2(RHSFunc, x0, n1, tFinal);
plot(linspace(0, 32, n1+1), sol1(:,1));
saveas(gcf,'Figures/01_20.png', 'png');
n2 = 101;
sol2 = AB2(RHSFunc, x0, n2, tFinal);
plot(linspace(0, 32, n2+1), sol2(:,1));
saveas(gcf,'Figures/01_21.png', 'png');
n3 = 301;
sol3 = AB2(RHSFunc, x0, n3, tFinal);
plot(linspace(0, 32, n3+1), sol3(:,1));
saveas(gcf,'Figures/01_22.png', 'png');

exactSolFunc = @(t) exp(-sigma/2*t).*cos(sqrt(15)/4 * t);
exactSol = exactSolFunc(linspace(0, 32, 1000));
plot(linspace(0, 32, 1000), exactSol);
saveas(gcf,'Figures/01_23.png', 'png');

plot(linspace(0, 32, n2+1), sol2(:,1), 'ro',...
    linspace(0, 32, n3+1), sol3(:,1), 'b+',...
    linspace(0, 32, 1000), exactSol, 'k');
saveas(gcf,'Figures/01_24.png', 'png');
