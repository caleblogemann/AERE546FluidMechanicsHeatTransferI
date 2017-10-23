%% 1
% x bounds
a = 0;
b = 1;
% y bounds
c = 0;
d = 1;
% number of x and y points
Nx = 201;
Ny = 201;
nGridCells = Nx*Ny;
% spacing
deltaX = (b - a)/(Nx-1);
deltaY = (d - c)/(Ny-1);

deltaT = 0.25*deltaX^2;
tFinal = 0.1;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;

% numbered first left to right, then bottom to top
nFunc = @(i, j) (j - 1)*Nx + i;
iFunc = @(n) mod(n, Nx) + Nx*(mod(n, Nx) == 0);
jFunc = @(n) (n - iFunc(n))/Nx + 1;
reshapeFunc = @(u) reshape(u, [Nx, Ny])';

% intial conditions
sigma = 0.1;
u0func = @(x, y) 1/(pi*sigma^2)*exp(-((x - 0.25)^2 + (y - 0.25)^2)/sigma^2);

RHSFunc = @(t, T) HeatRHS2D(t, T, Nx, Ny, deltaX, deltaY, reshapeFunc, nFunc);

x = linspace(a, b, Nx);
y = linspace(a, b, Ny);

u0 = zeros(nGridCells, 1);
for i = 2:Nx-1
    for j = 2:Ny-1
        u0(nFunc(i, j)) = u0func(x(i), y(j));
    end
end

rk2sol = RK2(RHSFunc, u0, nTimeSteps, deltaT);

initialMax = max(rk2sol(1,:));
for n = 2:nTimeSteps+1
    if (max(rk2sol(n,:)) < 0.01*initialMax)
        disp(n);
        disp((n-1)/nTimeSteps*tFinal);
        break;
    end
end

rk2solMatrix = reshapeFunc(rk2sol(0.02/deltaT,:));
contour(x, y, rk2solMatrix);
title('t = .02');
saveas(gcf,'Figures/04_01.png', 'png');

rk2solMatrix = reshapeFunc(rk2sol(nTimeSteps+1,:));
contour(x, y, rk2solMatrix)
title('t = .1');
saveas(gcf,'Figures/04_02.png', 'png');

%% 2
% x bounds
a = 0;
b = 1;
% y bounds
c = 0;
d = 1;
% number of x and y points
Nx = 201;
Ny = 201;
nGridCells = Nx*Ny;
% spacing
deltaX = (b - a)/(Nx-1);
deltaY = (d - c)/(Ny-1);

deltaT = deltaX^2;
tFinal = 0.25;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;

% numbered first left to right, then bottom to top
nFunc = @(i, j) (j - 1)*Nx + i;
iFunc = @(n) mod(n, Nx) + Nx*(mod(n, Nx) == 0);
jFunc = @(n) (n - iFunc(n))/Nx + 1;
reshapeFunc = @(u) reshape(u, [Nx, Ny])';

x = linspace(a, b, Nx);
y = linspace(a, b, Ny);
u0func = @(x, y) (y == 0)*sin(2*pi*x)^2;

u0 = zeros(nGridCells, 1);
for i = 1:Nx
    for j = 1:Ny
        u0(nFunc(i, j)) = u0func(x(i), y(j));
    end
end

g = @(t, x, y) (y == 0)*sin(2*pi*x).^2;

adiSol = ADI(u0, deltaX, deltaY, Nx, Ny, deltaT, nTimeSteps, g);

iter = 3;
for t = [0.002, 0.01, 0.04, 0.1, 0.25]
    rk2solMatrix = reshapeFunc(adiSol(t/deltaT+1,:));
    contour(x, y, rk2solMatrix);
    title(strcat('t = ',num2str(t)));
    xlabel('x');
    ylabel('y');
    saveas(gcf,strcat('Figures/04_0',num2str(iter),'.png'), 'png');
    iter=iter+1;
end

T = adiSol(:,nFunc(floor(Nx/2),floor(Nx/2)));
plot(0:deltaT:tFinal,T);
ylabel('T(1/2, 1/2)')
xlabel('time');
saveas(gcf, 'Figures/04_08.png', 'png');
