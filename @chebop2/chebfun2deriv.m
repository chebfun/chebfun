function deriv = chebfun2deriv( op )
%CHEBFUN2DDERIV Summary of this function goes here
%   Detailed explanation goes here

% For now, only work with linear operators. Hence, we can take any function on
% [-1,1] to evaluate the function handle to get the derivative information

u = chebfun2(@(x,y) x.*y);

opu = op(u);

deriv = opu.deriv;

end