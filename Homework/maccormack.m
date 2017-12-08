function [u] = maccormack(f, u0, deltaT, deltaX, nTimeSteps)
    nGridCells = length(u0);
    u = zeros(nTimeSteps+1, nGridCells);
    u(1, :) = u0;

    a = deltaT/deltaX;
    b = 0.5 * deltaT/deltaX;

    F = zeros(nGridCells,1);
    Ustar = zeros(nGridCells,1);

    for n = 1:nTimeSteps
        % Fluxes
        for j = 1:nGridCells
            F(j) = f(u(n, j));
        end

        % Ustar
        for j = 1:nGridCells-1
            Ustar(j) = u(n, j) - a*(F(j+1) - F(j));
        end
        % periodic boundary conditions
        Ustar(nGridCells) = u(n, 1);

        % Ustar Fluxes
        for j = 1:nGridCells-1;
            F(j) = f(Ustar(j));
        end

        % update solution
        % periodic boundary condition
        u(n+1, 1) = 0.5*(u(n, nGridCells) + Ustar(nGridCells));
        for j = 2:nGridCells
            u(n+1,j) = 0.5*(u(n, j) + Ustar(j)) - b*(F(j) - F(j-1));
        end
    end
end
