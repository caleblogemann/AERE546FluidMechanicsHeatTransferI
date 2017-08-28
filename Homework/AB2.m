function [result] = AB2(RHSFunc, x0, nTimeSteps, tFinal)
    nEquations = length(x0);
    result = zeros(nTimeSteps+1, nEquations);
    result(1,:) = x0;
    deltaT = tFinal/nTimeSteps;

    % take first step with explicit Euler method
    rhsOld = RHSFunc(0, result(1,:));
    result(2,:) = x0 + deltaT*rhsOld;

    for i = 2:nTimeSteps
        t = (i-1)*deltaT;
        rhsNew = RHSFunc(t, result(i,:));
        result(i+1, :) = result(i,:) + deltaT*(1.5*rhsNew + 0.5*rhsOld);
        rhsOld = rhsNew;
    end
end
