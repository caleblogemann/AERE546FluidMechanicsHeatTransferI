function [u, k, res] = sor3(lambda, u, x, y, b, tol, maxIter)
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
                if (i == I) % right boundary x = 2
                    uip1j = (16/9)*y(j)^2 - (32/9)*y(j) + 7/9;
                    uim1j = u(i-1, j);
                elseif (i == 1) % left boundary x = 0
                    uip1j = u(i+1,j);
                    uim1j = -y(j);
                else
                    uip1j = u(i+1,j);
                    uim1j = u(i-1, j);
                end

                if (j == J) % top boundary, y = 1
                    uijp1 = -1;
                    uijm1 = u(i,j-1);
                elseif (j == 1) % bottom boundary y = 0
                    uijp1 = u(i,j+1);
                    uijm1 = 0;
                else
                    uijp1 = u(i,j+1);
                    uijm1 = u(i,j-1);
                end

                if ( x(i) >= 1 && y(j) <= 1/4)
                    u(i, j) = 0;
                else
                    u(i, j) = (1 - lambda)*u(i, j) + lambda*(uip1j + uim1j + b^2*(uijp1 + uijm1))/(2*(1 + b^2));
                end
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
