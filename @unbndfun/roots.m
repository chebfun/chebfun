function r = roots(f, varargin)
%ROOTS   Roots of an UNBNDFUN in an unbounded domain.
%   ROOTS(F) returns the real roots of the UNBNDFUN F.
%
%   ROOTS(F, OPTIONS) modifies the default ROOTS properties, by passing the
%   OPTIONS to the rootfinding method of the ONEFUN of F.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call ROOTS@FUN.  Roots of an UNBNDFUN should always be pruned to try and
% remove spurious roots caused by fast decay of f near infinite endpoints.
if ( (nargin > 1) && isstruct(varargin{1}) )
    varargin{1}.filterEndpointRoots = true;
    r = roots@classicfun(f, varargin{:});
else
    r = roots@classicfun(f, varargin{:}, 'filterEndpointRoots', true);
end

% Do further filtering in case something was missed by the pruning filter in
% ROOTS@FUN():

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
