function [Au] = multByA(u, x, y, deltaX, deltaY)
    [I, J] = size(u);
    Au = zeros(I, J);

    for i = 1:I
        for j = 1:J
                % x = 1
                if (i == I)
                    uip1j = 0;
                    uim1j = u(i-1, j);
                % x = 0
                elseif (i == 1)
                    uip1j = u(i+1, j);
                    uim1j = sin(4*pi*y(j))^2;
                else
                    uip1j = u(i+1, j);
                    uim1j = u(i-1, j);
                end

                if (j == J)
                    uijp1 = 0;
                    uijm1 = u(i,j-1);
                elseif (j == 1)
                    uijp1 = u(i,j+1);
                    uijm1 = sin(pi*x(i))^2 + 5*sin(4*pi*x(i))^2 + 10*sin(8*pi*x(i))^2;
                else
                    uijp1 = u(i,j+1);
                    uijm1 = u(i,j-1);
                end

                Au(i, j) = (uip1j - 2*u(i, j) + uim1j)/deltaX^2 + (uijp1 - 2*u(i, j) + uijm1)/deltaY^2;
        end
    end
end
