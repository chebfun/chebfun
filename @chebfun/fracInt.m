function f = fracInt(f, mu)
%FRACINT  Fractional integral of a CHEBFUN. 
%   FRACINT(F, MU) gives the order MU fractional integral of a CHEBFUN object F.
%
%   Currently this only supports the situation where F is smooth (i.e., it has
%   no breakpoints or endpoint singularities) and on a finite domain.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract the fractional part of mu:
mu_int = floor(mu);
mu_frac = mu - mu_int;

% Computer the integer CUMSUMs:
f = cumsum(f, mu_int);

if ( mu_frac == 0 )
    return
end

% No piecewise support yet:
if ( numel(f.funs) > 1 )
    error('CHEBFUN:CHEBFUN:fracInt:breakpoints', ...
        'FRACINT() does not currently support piecewise functions.');
end

% Call BNDFUN/FRACINT():
f.funs{1} = fracInt(f.funs{1}, mu_frac);

f.pointValues = chebfun.getValuesAtBreakpoints(f.funs);

end