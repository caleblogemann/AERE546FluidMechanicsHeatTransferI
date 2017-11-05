function [u] = solvePoisson(multByA, rhs, I, J)
    A = zeros(I*J);
    iter = 0;
    for j = 1:J
        for i = 1:I
            iter = iter + 1;
            e = zeros(I, J);
            e(i, j) = 1;
            temp = multByA(e);
            for k = 1:J
                A(((k-1)*I+1):k*I,iter) = temp(:,k);
            end
        end
    end
    sol = A\rhs;
    u = zeros(I, J);
    for k = 1:J
        u(:,k) = sol(((k-1)*I+1):k*I);
    end
end

