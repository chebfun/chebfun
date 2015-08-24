function [Y, X] = max2( f )
%MAX2   Global maximum of a SEPARABLEAPPROX.
%   Y = MAX2(F) returns the global maximum of F over its domain. 
%   
%   [Y, X] = MAX2(F) returns the global maximum in Y and its location X.
%
%  This command may be faster if the OPTIMIZATION TOOLBOX is installed.
% 
% See also MIN2, MINANDMAX2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call MINANDMAX2:
[Y, X] = minandmax2(f);   

% Extract out maximum:
Y = Y(2); 
X = X(2,:); 

end
