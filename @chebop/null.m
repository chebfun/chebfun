function v = null(N, prefs, nullity)
%NULL   Null space of a linear CHEBOP.
%   Z = NULL(N) returns a CHEBMATRIX with orthonormal columns which span the
%   null space of the linear CHEBOP N. That is, N(Z) has negligible elements,
%   SIZE(Z, 2) is the nullity of N, and Z'*Z = I. N may contain linear boundary
%   conditions, but they will be treated as homogeneous. The nullity is
%   determined by comparing the differential order of the system and the number
%   of supplied boundary conditions.
%
%   NULL(N, PREFS) allows additional preferences to be passed via the 
%   CHEBOPPREF, PREFS.
%
%   NULL(n, PREFS, K) or NULL(N, K) attempts to find K null vectors. If the
%   nullity of N is determined to be less than K then a warning is thrown. This
%   is useful in situations where the nullity is known in advance and the
%   algorithm struggles to determine it automatically.
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

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    nullity = [];
end
if ( nargin < 2 )
    prefs = cheboppref();
elseif ( ~isa(prefs, 'cheboppref') )
    nullity = prefs;
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
v = null(L, prefs, nullity);

% Return a CHEBFUN rather than a CHEBMATRIX for scalar problems:
if ( isa(v, 'chebmatrix') && all(size(v, 1) == 1) )
    v = chebfun(v);
end

end
