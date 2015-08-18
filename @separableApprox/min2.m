function [Y, X] = min2( f )
%MIN2   Global minimum of a SEPARABLEAPPROX.
%   Y = MIN2(F) returns the global minimum of F over its domain. 
%   
%   [Y, X] = MIN2(F) returns the global minimum in Y and its location in X.  
%
%  This command may be faster if the OPTIMIZATION TOOLBOX is installed.
%
% See also MAX2, MINANDMAX2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call MINANDMAX2:
[Y, X] = minandmax2(f);

% Extract out minimum:
Y = Y(1); 
X = X(1,:);  

end
