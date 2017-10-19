function [result] = RK2(RHSFunc, x0, nTimeSteps, deltaT)
    n = length(x0);
    result = zeros(nTimeSteps+1, n);
    result(1,:) = x0;

    for i = 1:nTimeSteps
        t = (i-1)*deltaT;
        temp = result(i,:) + 0.5*deltaT*RHSFunc(t, result(i,:));
        result(i+1, :) = result(i,:) + deltaT*RHSFunc(t + 1/2*deltaT, temp);
        disp(i/nTimeSteps);
    end
end
