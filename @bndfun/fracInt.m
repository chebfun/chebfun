function f = fracInt(f, mu)
%FRACINT  Fractional integral of a BNDFUN. 
%   FRACINT(F, MU) gives the order MU fractional integral of a BNDFUN object F.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Scaling for domain:
scl = (diff(f.domain)/2)^mu;

if ( ~isa(f.onefun, 'singfun') )
    % Cast ONEFUN to a SINGFUN since the result must be singular:
    f.onefun = singfun(f.onefun);
end

% Fractional integral of the ONEFUN:
f.onefun = scl * fracInt(f.onefun, mu);

end
