function [result] = HeatRHS2D(t, T0, Nx, Ny, deltaX, deltaY, reshapeFunc, nFunc)
    result = zeros(1,Nx*Ny);
    T = reshapeFunc(T0);
    for i = 2:Nx-1
        for j = 2:Ny-1
            result(nFunc(i, j)) = (T(i+1,j) - 2*T(i, j) + T(i-1,j))/deltaX^2 + (T(i, j+1) - 2*T(i, j) + T(i,j-1))/deltaY^2;
        end
    end
end
