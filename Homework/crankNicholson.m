function [result] = crankNicholson(x0, nTimeSteps, deltaT, deltaX)
    n = length(x0);
    result = zeros(nTimeSteps+1, n);
    result(1,:) = x0;
    a = 2*deltaX^2/deltaT;

    mainDiagonal = -(a + 2)*ones(n, 1);
    lowerDiagonal = ones(n-1,1);
    upperDiagonal = ones(n-1,1);
    % zero flux on left boundary
    mainDiagonal(1) = -(a + 2) + 1;

    for i = 1:nTimeSteps
        RHS = -[result(i,2:end),1] - (a - 2)*result(i,:) - [result(i,1), result(i,1:end-1)];
        % value at right boundary is one
        RHS(end) = RHS(end) - 1;
        result(i+1, :) = tridiag(n, mainDiagonal, lowerDiagonal, upperDiagonal, RHS);
    end
end
