function [rhs] = upwindConvection(u, deltaX)
    nCells = length(u);
    rhs = zeros(1,nCells);

    nu = 1/deltaX;

    for j = 1:nCells
        jm1 = j-1;
        % periodic boundary conditions
        if (j == 1)
            jm1 = nCells;
        end
        rhs(j) = -u(j)*nu*(u(j) - u(jm1));
    end
end
