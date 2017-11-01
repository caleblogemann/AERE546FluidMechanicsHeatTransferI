function [u, k, res] = sor(lambda, u, x, y, n, tol, maxIter)
    [I, J] = size(u);
    uold = u;
    k = 0;
    mstop = 1;
    res = zeros(1,maxIter);
    residual0 = 0;
    while(mstop && k < maxIter)
        k = k + 1;
        for i = 1:I
            for j = 1:J
                if (i == I)
                    uip1j = 0;
                    uim1j = u(i-1, j);
                elseif (i == 1)
                    uip1j = u(i+1,j);
                    uim1j = sin(3*n*pi*y(j))^2;
                else
                    uip1j = u(i+1,j);
                    uim1j = u(i-1, j);
                end

                if (j == J)
                    uijp1 = 0;
                    uijm1 = u(i,j-1);
                elseif (j == 1)
                    uijp1 = u(i,j+1);
                    uijm1 = sin(3*n*pi*x(i))^2;
                else
                    uijp1 = u(i,j+1);
                    uijm1 = u(i,j-1);
                end

                u(i, j) = (1 - lambda)*u(i, j) + lambda*0.25*(uip1j + uim1j + uijp1 + uijm1);
            end
        end

        deltaU = u - uold;
        residual = sqrt(sum(sum(deltaU.^2))/(I*J));
        res(k) = residual;
        if (k == 1)
            residual0 = residual;
        else
            if (residual/residual0 <= tol)
                mstop = 0;
            else
                uold = u;
            end
        end
    end
    res = res(1:k);
end
