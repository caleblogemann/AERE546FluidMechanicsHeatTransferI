function [result] = RK2(RHSFunc, x0, nTimeSteps, tFinal)
    nEquations = length(x0);
    result = zeros(nTimeSteps+1, nEquations);
    result(1,:) = x0;

    deltaT = tFinal/nTimeSteps;

    for i = 1:nTimeSteps
        t = (i-1)*deltaT;
        temp = result(i,:) + 0.5*deltaT*RHSFunc(t, result(i,:));
        result(i+1, :) = result(i,:) + deltaT*RHSFunc(t + 1/2*deltaT, temp);
    end
end
