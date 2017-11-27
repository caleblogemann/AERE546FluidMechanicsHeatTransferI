function [u, k, res] = sor4(lambda, u, x, y, bt, tol, maxIter)
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
                    uip1j = -y(j) + 0.5;
                    uim1j = u(i-1, j);
                elseif (i == 1) % left boundary x = 0
                    uip1j = u(i+1,j);
                    uim1j = -y(j) + 0.5;
                else
                    uip1j = u(i+1,j);
                    uim1j = u(i-1, j);
                end

                if (j == J) % top boundary, y = 1
                    uijp1 = -0.5;
                    uijm1 = u(i,j-1);
                elseif (j == 1) % bottom boundary y = 0
                    uijp1 = u(i,j+1);
                    uijm1 = 0.5;
                else
                    uijp1 = u(i,j+1);
                    uijm1 = u(i,j-1);
                end

                w = 0;
                if (x(i) >= 1 && x(i) <= 1.3)
                    if (y(j) >= 0.35 && y(j) <= 0.5)
                        w = 50;
                    elseif (y(j) >= 0.5 && y(j) <= 0.65)
                        w = -50;
                    end
                end

                if ( x(i) >= 0.7 && x(i) <= 1 && y(j) >=0.35 && y(j) <= 0.65)
                    u(i, j) = 0;
                else
                    u(i, j) = (1 - lambda)*u(i, j) + lambda*(uip1j + uim1j + bt^2*(uijp1 + uijm1))/(2*(1 + bt^2)) + lambda*w;
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
