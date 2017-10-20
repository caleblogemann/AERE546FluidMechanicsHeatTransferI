for n = 1:10:nTimeSteps+1
    test = reshape(rk2sol(n,:), [Nx, Ny]);
    contour(x, y, test);
    pause(.01/nTimeSteps);
end