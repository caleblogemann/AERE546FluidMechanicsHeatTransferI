%% 1 a
topBoundary = @(x) (1 - 0.4*exp(-8*(0.6 - x).^2)).*(x <= 0.6) + 0.6*(x > 0.6);
dtopBoundary = @(x) -0.4*exp(-8*(0.6 - x).^2).*(16*(0.6 - x))*(x <= 0.6) + 0*(x > 0.6);
bottomBoundary = @(x) 0;
dbottomBoundary = @(x) 0;

Nx = 101;
Ny = 121;
deltaX = 2/Nx;
xIndex = @(i) -1 + deltaX*i;
deltaY = @(i) topBoundary(xIndex(i))/Ny;
yIndex = @(i, j) deltaY(i)*j;

hold on;
x = xIndex(0:Nx);
plot(x, topBoundary(x),'k');
plot(x, bottomBoundary(x),'k');
for i = 1:Nx
    y = linspace(bottomBoundary(xIndex(i)), topBoundary(xIndex(i)), 1000);
    [x, ~] = hermiteInterpolation(bottomBoundary(xIndex(i)), xIndex(i), -dbottomBoundary(xIndex(i)), topBoundary(xIndex(i)), xIndex(i), -dtopBoundary(xIndex(i)), y);
    plot(x, y, 'k');
end

for j = 1:Ny
    x = zeros(Nx, 1);
    y = zeros(Nx, 1);
    for i = 1:Nx
        y(i) = yIndex(i, j);
        [x(i), ~] = hermiteInterpolation(bottomBoundary(xIndex(i)), xIndex(i), -dbottomBoundary(xIndex(i)), topBoundary(xIndex(i)), xIndex(i), -dtopBoundary(xIndex(i)), y(i));
    end
    plot(x, y, 'k');
end
hold off
xlabel('x');
ylabel('y');
title('H grid');
saveas(gcf, 'Figures/07_01.png', 'png');

%% 1b
Ntheta = 101;
Nr = 81;
outerBoundaryX = @(theta) 3*cos(theta);
outerBoundaryY = @(theta) 3*sin(theta);
outerBoundaryR = @(theta) 3;
innerBoundaryX = @(theta) cos(theta);
innerBoundaryY = @(theta) 0.3*sin(theta);
innerBoundaryR = @(theta) sqrt(innerBoundaryX(theta).^2 + innerBoundaryY(theta).^2);
dydxInnerBoundary = @(theta) -0.3*cos(theta)./sin(theta);
dydxOuterBoundary = @(theta) -cos(theta)./sin(theta);
theta = linspace(0, 2*pi, Ntheta);
figure;
hold on;
plot(outerBoundaryX(theta), outerBoundaryY(theta), 'k', innerBoundaryX(theta), innerBoundaryY(theta), 'k');
for i = 1:Ntheta
    x = linspace(innerBoundaryX(theta(i)), outerBoundaryX(theta(i)), 1000);
    [y, dy] = hermiteInterpolation(innerBoundaryX(theta(i)), innerBoundaryY(theta(i)), ...
        -1/dydxInnerBoundary(theta(i)), outerBoundaryX(theta(i)), ...
        outerBoundaryY(theta(i)), -1/dydxOuterBoundary(theta(i)), x);
    plot(x, y, 'k');
end

deltaR = @(theta) (outerBoundaryR(theta) - innerBoundaryR(theta))/Nr;
for j = 1:Nr
    for i = 1:Ntheta-1
        
        x = linspace(innerBoundaryX(theta(i)), outerBoundaryX(theta(i)), 1000);
        [y, dy] = hermiteInterpolation(innerBoundaryX(theta(i)), innerBoundaryY(theta(i)), ...
            -1/dydxInnerBoundary(theta(i)), outerBoundaryX(theta(i)), ...
            outerBoundaryY(theta(i)), -1/dydxOuterBoundary(theta(i)), x);
        r = sqrt(x.^2 + y.^2);
        dR1 = deltaR(theta(i));
        r1 = innerBoundaryR(theta(i)) + dR1*j;
        [~, index] = min(abs(r - r1));
        
        x1 = x(index);
        y1 = y(index);
        dy1 = dy(index);

        x = linspace(innerBoundaryX(theta(i+1)), outerBoundaryX(theta(i+1)), 1000);
        [y, dy] = hermiteInterpolation(innerBoundaryX(theta(i+1)), innerBoundaryY(theta(i+1)), ...
            -1/dydxInnerBoundary(theta(i+1)), outerBoundaryX(theta(i+1)), ...
            outerBoundaryY(theta(i+1)), -1/dydxOuterBoundary(theta(i+1)), x);
        r = sqrt(x.^2 + y.^2);
        dR2 = deltaR(theta(i+1));
        r2 = innerBoundaryR(theta(i+1)) + dR2*j;
        [~, index] = min(abs(r - r2));
        
        x2 = x(index);
        y2 = y(index);
        dy2 = dy(index);

%         x = linspace(x1, x2, 100);
%         [y, dy] = hermiteInterpolation(x1, y1, -1/dy1, x2, y2, -1/dy2, x);
        plot([x1,x2], [y1,y2], 'k');
    end
end
title('O-grid');
hold off;
saveas(gcf, 'Figures/07_02.png', 'png');

%% problem 1(c)
Ntheta = 101;
Nr = 71;
outerBoundaryX = @(theta) 1.8 - 1.5*cos(theta);
outerBoundaryY = @(theta) 0.75*sin(theta);
outerBoundaryR = @(theta) sqrt(outerBoundaryX(theta).^2 + outerBoundaryY(theta).^2);
%innerBoundaryX = @(theta) cos(theta);
innerBoundaryYLower = @(x) 0*(0.8 <= x && x <= 1.8) - 4*(0.8 - x).^2.*sqrt(x - 0.55).*(0.55 <= x && x <= 0.8);
innerBoundaryYUpper = @(x) 0*(0.8 <= x && x <= 1.8) + 4*(0.8 - x).^2.*sqrt(x - 0.55).*(0.55 <= x && x <= 0.8);
%innerBoundaryR = @(theta) sqrt(innerBoundaryX(theta).^2 + innerBoundaryY(theta).^2);
dydxInnerBoundaryUpper = @(x) 0*(0.8 <= x && x <= 1.8) + (2*(0.8 - x).^2./sqrt(x-0.55) - 8*(0.8 - x).*sqrt(x -0.55)).*(0.55 <= x && x <= 0.8)
dydxInnerBoundaryLower = @(x) 0*(0.8 <= x && x <= 1.8) - (2*(0.8 - x).^2./sqrt(x-0.55) - 8*(0.8 - x).*sqrt(x -0.55)).*(0.55 <= x && x <= 0.8)
dydxOuterBoundary = @(theta) 0.75*cos(theta)./(1.5*sin(theta));
theta = linspace(-pi/2, pi/2, Ntheta);
Nx = 101;
xIndex = @(i) 0.55 + 1.25*((2*i - Nx - 1)/(Nx - 1)).^2;

figure;
hold on;
plot(outerBoundaryX(theta), outerBoundaryY(theta), 'k');

for i = 1:(Nx - 1)/2
    x = linspace(xIndex(i), outerBoundaryX(theta(i)), 1000);
    [y, dy] = hermiteInterpolation(xIndex(i), innerBoundaryYLower(xIndex(i)), ...
        -1/dydxInnerBoundaryLower(xIndex(i)), outerBoundaryX(theta(i)), ...
        outerBoundaryY(theta(i)), -1/dydxOuterBoundary(theta(i)), x);
    plot(x, y, 'k');
end

%for i = (Nx + 1)/2:Nx
%end
