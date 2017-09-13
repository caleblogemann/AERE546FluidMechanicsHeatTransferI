function [y] = tridiag(n, a, b, c, v);
    % solve a tridiagonal system with Gaussian elimination
    % Also known as Thomas Algorithm
    % (a_1 c_1   0   0                  0  )(   y_1   )  (   v_1   )
    % (b_1 a_2 c_2   0                  0  )(   y_2   )  (   v_2   )
    % (  0 b_2 a_3 c_3                  0  )(   y_3   )= (   v_3   )
    % (  0                              0  )(         )= (         )
    % (  .                  a_{n-1} c_{n-1})( y_{n-1} )  ( v_{n-1} )
    % (  0                  b_{n-1}    a_n )(   y_n   )  (   v_n   )

    % create array zero to store solutions
    y = zeros(size(v));

    % eliminate b_i's 
    for(i=1:n-1)
        a(i+1) = a(i+1) + c(i)*(-b(i)/a(i));
        v(i+1) = v(i+1) + v(i)*(-b(i)/a(i));
    end

    % solve for y_n
    y(n) = v(n)/a(n);

    for(i=(n-1):-1:1)
        y(i) = (v(i) - c(i)*y(i+1))/a(i);
    end

end
