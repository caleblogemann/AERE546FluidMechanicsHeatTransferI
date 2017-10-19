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
deltaX = (b - a)/Nx;
deltaY = (d - c)/Ny;

deltaT = 0.1*deltaX^2;
tFinal = 0.02;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;

% numbered first left to right, then bottom to top
nFunc = @(i, j) (j - 1)*Nx + i;
iFunc = @(n) mod(n, Nx) + Nx*(mod(n, Nx) == 0);
jFunc = @(n) (n - iFunc(n))/Nx + 1;

% intial conditions
sigma = 0.1;
u0func = @(x, y) 1/(pi*sigma^2)*exp(-((x - 0.25)^2 + (y - 0.25)^2)/sigma^2);

RHSFunc = @(t, T) arrayfun(@(i) (T(min(i+1, nGridCells))*(i+1 <= nGridCells || mod(i, Nx)==0) - ...
    2*T(i) + T(max(i-1, 1))*(i-1 >= 1 || mod(i,Nx) == 1))/deltaX^2 + ...
    (T(min(i+Nx, nGridCells))*(i+Nx <= nGridCells) - 2*T(i) + ...
    T(max(i-Nx, 1))*(i - Nx >= 1))/deltaY^2, 1:length(T));

x = linspace(a, b, Nx);
y = linspace(a, b, Ny);

u0 = zeros(nGridCells, 1);
for i = 1:Nx
    for j = 1:Ny
        u0(nFunc(i, j)) = u0func(x(i), y(j));
    end
end

rk2sol = RK2(RHSFunc, u0, nTimeSteps, deltaT);

rk2solMatrix = reshape(rk2sol(nTimeSteps+1,:), [Nx, Ny]);
contour(rk2solMatrix);
