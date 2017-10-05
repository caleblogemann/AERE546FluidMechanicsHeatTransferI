function [result] = forwardEuler(RHSFunc, x0, nTimeSteps, deltaT)
    n = length(x0);
    result = zeros(nTimeSteps+1, n);
    result(1,:) = x0;

    for i = 1:nTimeSteps
        t = (i-1)*deltaT;
        result(i+1, :) = result(i,:) + deltaT*RHSFunc(t, result(i,:));
    end
end
