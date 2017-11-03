function [x, iter, res] = conjugateGradient(multByA, x0, rhs, tol, maxIterations)
    x = x0;
    Ax = multByA(x);
    r = rhs - Ax;
    rsold = sum(dot(r,r));
    
    res = zeros(1,maxIterations);

    p = r;

    iter = 0;
    while(iter < maxIterations)
        iter = iter + 1;

        Ap = multByA(p);
        pAp = sum(dot(p, Ap));
        alpha = rsold/pAp;
        x = x + alpha*p;
        r = r - alpha*Ap;

        rsnew = sum(dot(r,r));
        res(iter) = rsnew;
        if(rsnew/res(1) < tol)
            rsold = rsnew;
            break;
        else
            p = r + rsnew/rsold*p;
            rsold = rsnew;
        end
    end
    if (iter == maxIterations && sqrt(rsold) > tol)
        disp('Conjugate Gradient did not converge');
    end
    res = res(1:iter);
end
