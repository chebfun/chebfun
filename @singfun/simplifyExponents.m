function f = simplifyExponents(f)
%SIMPLIFYEXPONENTS  Simplify the exponents of a SINGFUN.
%   F = SIMPLIFYEXPONENTS(F) returns a SINGFUN which has both exponents less
%   than 1 by absorbing the integer part of any boundary exponents larger than
%   or equal to 1 into its smoothPart.
%
% See also EXTRACTBOUNDARYROOTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab the exponents:
exps = get(f, 'exponents');

% Grab the indice for exponents larger or equal to 1:
ind = ( exps >= 1 );

% Both exponents are less than 1:
if ( ~any( ind ) )
    return
end

% Sort out the new exponents and the order of the boundary roots which need to
% be absorbed into the smoothPart:
newExps = exps;
newExps(ind) = exps(ind) - floor(exps(ind));
pow = exps - newExps;

% Compute the factor from function roots: 
mult = singfun.constructSmoothPart(@(x) (x+1).^pow(1).*(1-x).^pow(2), [], []);
f.smoothPart = f.smoothPart.*mult;

% Update the exponents:
f.exponents = newExps;

end
