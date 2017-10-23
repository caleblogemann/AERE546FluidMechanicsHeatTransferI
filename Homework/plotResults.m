for n = 1:10:nTimeSteps+1
    test = reshapeFunc(adiSol(n,:));
    contour(x, y, test);
    pause(.0000001);
end