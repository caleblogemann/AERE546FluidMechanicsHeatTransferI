function [rhs] = central(u, f, deltaX)
    nCells = length(u);
    rhs = zeros(1,nCells);

    % vector to store fluxes at cell interfaces
    % F(j) is flux at cell interface x_{j-1/2} on left hand side
    F = zeros(nCells + 1);
    nu = 1/deltaX;

    % Compute fluxes
    for j = 1:nCells
        jm1 = j-1;
        % periodic boundary conditions
        if (j == 1)
            jm1 = nCells;
        end
        F(j) = 0.5*(f(u(jm1)) + f(u(j)));
    end

    for j = 1:nCells
        jp1 = j+1;
        if(j == nCells)
            jp1 = 1;
        end
        rhs(j) = nu*(F(j) - F(jp1));
    end
end
