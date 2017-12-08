function [u] = upwindDownwind(u0, f, p, deltaX, deltaT, nTimeSteps)

    nCells = length(u0);
    u = zeros(nTimeSteps+1, nCells);
    u(1,:) = u0;

    % vector to store fluxes at cell interfaces
    % F(j) is flux at cell interface x_{j-1/2} on left hand side
    F = zeros(nCells + 1);
    nu = deltaT/deltaX;

    for n = 1:nTimeSteps
        % Compute fluxes
        for j = 1:nCells
            jm1 = j-1;
            % periodic boundary conditions
            if (j == 1)
                jm1 = nCells;
            end
            F(j) = p*f(u(n, jm1)) + (1 - p)*f(u(n, j));
        end

        for j = 1:nCells
            jp1 = j+1;
            if(j == nCells)
                jp1 = 1;
            end
            u(n+1, j) = u(n, j) - nu*(F(jp1) - F(j));
        end
    end
end
