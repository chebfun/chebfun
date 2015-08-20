function g = fracInt(f, mu)
%FRACINT  Fractional integral of a SINGFUN. 
%   FRACINT(F, MU) gives the order MU fractional integral of a SINGFUN object F.
%
%   Currently this only supports the situation where F is  is smooth at the
%   right boundary (i.e., F.exponents(2) = 0).

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Attempt to simplify the exponents as we must have zero at the right.
f = simplifyExponents(f);
exps = f.exponents;

if ( exps(2) ~= 0 )
    error('CHEBFUN:SINGFUN:fracInt:exponents', ...
        ['FRACINT(F, MU) currently only supports the situation when F is ', ...
         'smooth at the right boundary.']);
end

g = f;
% Compute the smooth part:
g.smoothPart = fracInt(f.smoothPart, mu, exps(1));
% Update the exponents:
g.exponents = f.exponents + [mu, 0];

if ( exps(1) ~= 0 )
    % Simplify the exponents again.
    g = simplifyExponents(g);
end

end
