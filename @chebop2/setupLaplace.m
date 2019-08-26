function N = setupLaplace( dom )
%SETUPLAPLACE( DOM )  Construct a chebop2 object for Laplace operator.
% A small piece of code that is faster than calling the chebop2
% constructor for forming the Laplace operator on DOM.

    N = chebop2();
    N.domain = dom;
    N.op = @(u) lap(u);
    N.coeffs = [ 0 0 1 ; 0 0 0 ; 1 0 0 ];
    N.xorder = 2;
    N.yorder = 2;
    N.rhs = 0;

end