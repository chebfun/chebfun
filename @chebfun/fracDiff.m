function f = fracDiff(f, mu, type)
%FRACDIFF  Fractional derivative of a CHEBFUN. 
%   FRACINT(F, MU) gives the order MU Riemann-Liouville fractional derivative of
%   a CHEBFUN object F.
%
%   FRACINT(F, MU, 'Caputo') instead uses the Caputo definition of the
%   fractional derivative.
%
%   Currently this only supports the situation where F is smooth (i.e., it has
%   no breakpoints or endpoint singularities) and on a finite domain.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to Riemann-Liouville:
if ( nargin < 3  )
    type = 'RL';
end

% Extract the fractional part:
n = ceil(mu);

if ( n == mu )
    f = diff(f, n);
    return
end

% No piecewise support yet:
if ( numel(f.funs) > 1 )
    error('CHEBFUN:CHEBFUN:fracDiff:breakpoints', ...
        'FRACDIFF does not currently support piecewise functions.');
end

if ( strcmpi(type, 'Caputo') )
    % Caputo:
    f = diff(f, n);
    f = fracInt(f, n - mu); 
    
else
    % Riemann-Liouville:
    f = fracInt(f, n - mu);
    f = diff(f, n);
    
end

end
