function [x] = gaussSeidel(A, b, x0, tol, maxIter)
    % A = L + U
    % Ax = b
    % (L + U)x = b
    % Lx + Ux = b
    % Lx = -Ux + b
    % x = -L^{-1} U x + L^{-1}b
    % make iteration
    % x_{k+1} = -L^{-1} U x_k + L^{-1}b

    % For Gauss-Seidel
    % M is lower triangular part of A including diagonal
    % N = A - M

    L = tril(A);
    U = A - L;
    k = 0;
    while(norm(A*x0 - b) > tol && k < maxIter)
        % just backsolving
        x = L\(-U*x0 + b);

        k = k + 1;
        x0 = x;
    end
end
