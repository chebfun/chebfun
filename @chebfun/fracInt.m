function f = fracInt(f, mu)
%FRACINT  Fractional integral of a CHEBFUN. 
%   FRACINT(F, MU) gives the order MU fractional integral of a CHEBFUN object F.
%
%   Currently this only supports the situation where F is smooth (i.e., it has
%   no breakpoints or endpoint singularities) and on a finite domain.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% No piecewise support yet:
if ( numel(domain(f)) > 2 )
    error('CHEBFUN:CHEBFUN:fracInt:breakpoints', ...
        'FRACINT does not currently support piecewise functions.');
end

% Extract the fractional part of mu:
mu_int = floor(mu);
mu_frac = mu - mu_int;

% Computer the integer CUMSUMs:
if ( mu_int > 0 )
    % Ensure mu is in [0, 1):
    f = cumsum(f, mu_int);
    f = fracInt(f, mu_frac);
    return
elseif ( mu_frac == 0 )
    % Nothing more to do (non-fractional integral).
    return
end

% Quasimatrix / Array-valued support:
% TODO: This should be overhauled once SINGFUN supports array-valuedness.
m = numColumns(f);
if ( m > 1 )
    % Convert to cell-array of scalar CHEBFUNs.
    f = mat2cell(f);    
    % Loop over columns:
    for k = 1:m
        f{k} = fracInt(f{k}, mu_frac);
    end
    % Convert from cell array to quasimatrix.
    f = cat(2-f{1}.isTransposed, f{:});
    return
end

% From here on f is scalar-valued:

% Call BNDFUN/FRACINT():
f.funs{1} = fracInt(f.funs{1}, mu_frac);

% Update the pointValues:
f.pointValues = chebfun.getValuesAtBreakpoints(f.funs);

end
