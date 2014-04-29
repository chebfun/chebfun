function r = roots(f, varargin)
%ROOTS   Roots of an UNBNDFUN in an unbounded domain.
%   ROOTS(F) returns the real roots of the UNBNDFUN F.
%
%   ROOTS(F, OPTIONS) modifies the default ROOTS properties, by passing the
%   OPTIONS to the rootfinding method of the ONEFUN of F.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call ROOTS@FUN:
r = roots@classicfun(f, varargin{:});

% Get rid of spurious roots which are caused by fast decay of function defined
% in an unbounded domain:

% Set a threshold for the 'farfield':
farfield = 1e-1/eps;

ends = get(f, 'domain');

if ( isinf(ends(1)) )
    mask = ( r < -farfield );
    r(mask) = [];
end

if ( isinf(ends(2)) )
    mask = ( r > farfield );
    r(mask) = [];
end

end
