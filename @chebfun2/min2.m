function [Y, X] = min2( f )
%MIN2 Global minimum of a chebfun2. 
%
% Y = MIN2(F) returns the global minium of F.
% 
% [Y, X]=MIN2(F) returns the global minimum of F and its coordinates in 
%    X = (X(1), X(2)). 
%
% For high accuracy results this command requires the Optimization Toolbox.
%
% See also MAX2, MINANDMAX2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Call minandmax2:
[Y, X] = minandmax2( f );
% Extract out minimum:
Y = Y( 1 ); 
X = X( 1, : );  

end