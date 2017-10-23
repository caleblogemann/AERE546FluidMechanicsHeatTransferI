function [u] = ADI(u0, deltaX, deltaY, Nx, Ny, deltaT, nTimeSteps, g)
    % to change from natural row wise ordering to natural
    % column wise ordering or vice versa
    % uCol = uRow(rowToCol) or uRow = uCol(colToRow)
    rowToCol = repmat(0:Ny:Nx*Ny-Ny, 1, Ny) + kron(1:Ny,ones(1, Nx));
    colToRow = repmat(0:Nx:Nx*Ny-Nx, 1, Nx) + kron(1:Nx,ones(1, Ny));

    % matrix to represent diffusion in row/x direction
    % Dx*uRow = Dx2*uRow
    e = ones(Nx, 1);
    tridiagonal = 1/(deltaX^2)*spdiags([e, -2*e, e], [-1, 0, 1], Nx, Nx);
    Dx = kron(speye(Ny), tridiagonal);

    % matrix to represent diffusion in col/y direction
    % Dy*uCol = Dx2*uCol
    e = ones(Ny, 1);
    tridiagonal = 1/(deltaY^2)*spdiags([e, -2*e, e], [-1, 0, 1], Ny, Ny);
    Dy = kron(speye(Nx), tridiagonal);

    I = speye(Nx*Ny);

    % solution at all times
    % stored in rowwise numbering
    u = zeros(nTimeSteps, Nx*Ny);
    u(1,:) = u0;

    % create index arrays for boundary
    % indices for left boundary in row-wise ordering or
    zeroIndicesRow = 1:Nx:Nx*Ny;
    % indices for right boundary in row-wise ordering or
    oneIndicesRow = Nx:Nx:Nx*Ny;
    % indices for bottom boundary in column-wise ordering
    zeroIndicesCol = 1:Ny:Nx*Ny;
    % indices for top boundary in column-wise ordering
    oneIndicesCol = Ny:Ny:Nx*Ny;
    for n = 1:nTimeSteps
        % First stage
        % (I + deltaT/2 Dy2)*u
        % rhs in column wise ordering
        rhs = (I + deltaT/2*Dy)*u(n, :)';
        rhs = rhs(rowToCol);
        % add boundary conditions for y-direction at time t = tn
        bottomBoundary = deltaT/(2*deltaY^2)*g(deltaT*(n-1), deltaX*(0:Nx-1), 0);
        topBoundary = deltaT/(2*deltaY^2)*g(deltaT*(n-1), deltaX*(0:Nx-1), 1);
        rhs(zeroIndicesCol) = rhs(zeroIndicesCol) + bottomBoundary';
        rhs(oneIndicesCol) = rhs(oneIndicesCol) + topBoundary';

        % change to row-wise ordering
        rhs = rhs(colToRow);
        % add boundary conditions for x-direction at time t = tn + deltaT/2
        leftBoundary = deltaT/(2*deltaX^2)*g(deltaT*(n-1)+deltaT/2, 0, deltaY*(0:Ny-1));
        rightBoundary = deltaT/(2*deltaX^2)*g(deltaT*(n-1)+deltaT/2, 1, deltaY*(0:Ny-1));
        rhs(zeroIndicesRow) = rhs(zeroIndicesRow) + leftBoundary';
        rhs(oneIndicesRow) = rhs(oneIndicesRow) + rightBoundary';
        % solve for uStar which approximates u at t = tn + k/2
        % uStar is in row-wise ordering
        uStar = (I - deltaT/2*Dx)\rhs;

        % second stage
        % rhs row-wise ordering
        rhs = (I + deltaT/2*Dx)*uStar;
        % add boundary conditions for x-direction at time t = tn + k/2
        % left and right boundaries same as before
        rhs(zeroIndicesRow) = rhs(zeroIndicesRow) + leftBoundary';
        rhs(oneIndicesRow) = rhs(oneIndicesRow) + rightBoundary';
        % change to column-wise ordering
        rhs = rhs(rowToCol);
        % add boundary conditions for y-direction at time t = tn + k
        bottomBoundary = deltaT/(2*deltaY^2)*g(deltaT*n, deltaX*(0:Nx-1), 0);
        topBoundary = deltaT/(2*deltaY^2)*g(deltaT*n, deltaX*(0:Nx-1), 1);
        rhs(zeroIndicesCol) = rhs(zeroIndicesCol) + bottomBoundary';
        rhs(oneIndicesCol) = rhs(oneIndicesCol) + topBoundary';
        % solve for u at time t = tn + k
        % u is in row wise ordering
        temp = ((I - deltaT/2*Dy)\rhs);
        u(n+1,:) = temp(colToRow)';
    end
end

