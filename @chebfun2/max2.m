function [Y, X]=max2(f)
%MAX2 Global maximum of a chebfun2
%
% Y = MAX2(F) returns the global maximum of F over its domain. 
%   
% [Y X] = MAX2(F) returns the global maximum in Y and its location X.  
%
% For high accuracy results this command requires the Optimization Toolbox.
% 
% See also MIN2, MINANDMAX2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Call minandmax2:
[Y, X] = minandmax2( f );   

% Extract out maximum
Y = Y(2); 
X = X(2,:); 

end