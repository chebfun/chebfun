function [Y, X] = min2( f )
%MAX2   Global maximum of a CHEBFUN2.
%   Y = MIN2(F) returns the global minimum of F over its domain. 
%   
%   [Y, X] = MIN2(F) returns the global minimum in Y and its location in X.  
%
%  This command may be faster if the OPTIMIZATION TOOLBOX is installed.
%
% See also MIN2, MINANDMAX2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call MINANDMAX2:
[Y, X] = minandmax2(f);

% Extract out minimum:
Y = Y(1); 
X = X(1,:);  

end
