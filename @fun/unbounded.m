function m = unbounded(ends)
%UNBOUNDEDMAP creates a map structure for chebfuns
%   M = UNBOUNDEDMAP(ENDS) returns a structure that defines a nonlinear map
%   from [-1 1] to the unbounded domain [ENDS(1) ENDS(2)].
%   The structure MAP consists of three function handles and one string.
%   M.FOR is a function that maps [-1,1] to ENDS.
%   M.INV is the inverse map.
%   M.DER is the derivative of the map defined in MAP.FOR
%   M.NAME is a string that identifies the map.

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.:

% The domain:
a = ends(1); b = ends(2);

% initialise the map structure:
m = struct('for', [], 'inv', [], 'der', [], 'name', 'unbounded', 'par', [a b]);

% Fixed map parameters:
s = 1;
c = 0;

% Deal with the different cases:
if ( a == -inf && b == inf )
    
    m.for = @(y) 5*s*y./(1 - min(y.^2, 1)) + c;
    m.inv = @(x) 2*x./(5*s + sqrt(25*s^2 + 4*x.^2));
    m.der = @(y) 5*s*(1 + y.^2)./(1 - y.^2).^2;

elseif ( a == -inf )
    
     m.inv = @(x) (15*s + x - b)./(15*s - x + b);
     m.for = @(y) 15*s*(y - 1)./(y + 1) + b;
     m.der = @(y) 15*s*2./(y + 1).^2;

elseif ( b == inf )
    
     m.inv = @(x) (-15*s + x - a)./(15*s + x - a);
     m.for = @(y) 15*s*(y + 1)./(1 - y)+a;
     m.der = @(y) 15*s*2./(y - 1).^2;

else
    
    error('FUN:unbounded:input', 'Error: Check input')
    
end
