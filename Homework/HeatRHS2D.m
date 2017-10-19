function [result] = HeatRHS2D(t, T)
    n = length(T);
    result = zeros(n, 1);
    for i = 1:n
        result(i) = 
    end
    
RHSFunc = @(t, T) arrayfun(@(i) (T(i+1) - 2*T(i) + T(i-1))/deltaX^2 + (T(max(i+Nx) - 2*T(i) + T(i - Nx))/deltaY^2, 1:length(T));
end
