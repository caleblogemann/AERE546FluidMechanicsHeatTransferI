function [y, dy] = hermiteInterpolation(x1, y1, dy1, x2, y2, dy2, x)
    a = y1;
    b = dy1;
    temp = (y2 - y1)/(x2 - x1);
    c = (temp - dy1)/(x2 - x1);
    d = ((dy2 - temp)/(x2 - x1) - c)/(x2 - x1);

    yfun= @(x) a + b*(x - x1) + c*(x - x1).*(x - x1) + d*(x - x1).*(x - x1).*(x - x2);
    dyfun = @(x) b + 2*c*(x - x1) + d*(2*(x - x1).*(x - x2) + (x - x1).^2);
    y = yfun(x);
    dy = dyfun(x);
end
