function v = null(N, prefs)
%NULL   Null space of a linear CHEBOP.
%   Z = NULL(N) returns a CHEBMATRIX with orthonormal columns which span the
%   null space of the linear CHEBOP N. That is, N(Z) has negligible elements,
%   SIZE(Z, 2) is the nullity of N, and Z'*Z = I. N may contain linear boundary
%   conditions, but they will be treated as homogeneous.
%
%   NULL(A, PREFS) allows additional preferences to be passed via the 
%   CHEBOPPREF, PREFS.
%
%   Systems of equations are not yet supported.
%
% Example 1:
%   N = chebop(@(u) diff(u), [0, pi]);
%   V = null(N);
%   norm(N(V))
%
% Example 2:
%   N = chebop(@(x, u) 0.2*diff(u, 3) - sin(3*x).*diff(u));
%   N.rbc = 1;
%   V = null(N)
%   norm(N(V))
%   plot(V)
%
% See also CHEBOP/SVDS, CHEBOP/EIGS, NULL.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    prefs = cheboppref();
end

% Linearize and check whether the CHEBOP is linear:
[L, ignored, fail] = linop(N); %#ok<ASGLU>
if ( fail )
    error('CHEBFUN:CHEBOP:null:nonlinear', ...
        ['The operator appears to be nonlinear.\n', ...
         'NULL() supports only linear CHEBOP instances.']);
end

% Determine the discretization:
prefs = determineDiscretization(N, L, prefs);

% Call LINOP/NULL:
v = null(L, prefs);

% Return a CHEBFUN rather than a CHEBMATRIX for scalar problems:
if ( isa(v, 'chebmatrix') && all(size(v, 1) == 1) )
    v = chebfun(v);
end

end
