function deriv = chebfun2deriv( op )
%CHEBFUN2DDERIV   
% 
% A few operations in Chebfun2 store AD information (see plus, diff, uminus)
% This allows us to pull out the coefficients to a partial differential
% operator. 

% For now, only work with linear operators. Hence, we can take any function on
% [-1,1] to evaluate the function handle to get the derivative information

% Evaluate the operator at a chebfun2 and use the AD information to 
% determine the coefficients of the partial differential operator: 

u = chebfun2(@(x,y) x.*y);
opu = op(u);
deriv = opu.deriv;

end