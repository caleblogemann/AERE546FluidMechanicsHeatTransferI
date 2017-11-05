function [Au] = multByA2(u, x, y, deltaX, deltaY)
    [I, J] = size(u);
    Au = zeros(I, J);

    for i = 1:I
        for j = 1:J
                % x = 1
                if (i == I)
                    uip1j = -y(j);
                    uim1j = u(i-1, j);
                % x = 0
                elseif (i == 1)
                    uip1j = u(i+1, j);
                    uim1j = -y(j);
                else
                    uip1j = u(i+1, j);
                    uim1j = u(i-1, j);
                end

                % y == 1
                if (j == J)
                    uijp1 = -1;
                    uijm1 = u(i,j-1);
                % y == 0
                elseif (j == 1)
                    uijp1 = u(i,j+1);
                    uijm1 = 0;
                else
                    uijp1 = u(i,j+1);
                    uijm1 = u(i,j-1);
                end

                Au(i, j) = exp(-4*x(i))*((uip1j - uim1j)/(2*deltaX) - (uip1j - 2*u(i,j) + uim1j)/(2*deltaX^2)) ...
                + exp(-4*y(j))*((uijp1 - uijm1)/(2*deltaY) - (uijp1 - 2*u(i, j) + uijm1)/(2*deltaY^2));
        end
    end
end
