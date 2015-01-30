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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Default to Riemann-Liouville:
if ( nargin < 3  )
    type = 'RL';
end

% Extract the fractional part:
mu_int = floor(mu);
mu_frac = mu - mu_int;

if ( mu_frac == 0 )
    f = diff(f, mu_int);
    return
end

% No piecewise support yet:
if ( numel(f.funs) > 1 )
    error('CHEBFUN:CHEBFUN:fracInt:breakpoints', ...
        'FRACINT() does not currently support piecewise functions.');
end

if ( strcmpi(type, 'Caputo') )
    % Caputo:
    f = fracInt(f, 1 - mu_frac);
    f = diff(f, mu_int + 1);
    
else
    % Riemann-Liouville:
    f = diff(f, mu_int + 1);
    f = fracInt(f, 1 - mu_frac); 
    
end

end