function varargout = polyeigs(N, varargin)
%POLYEIGS   Solve a CHEBOP polynomial eigenvalue problem.
% [X,E] = POLYEIG(A0,A1,..,Ap,K) solves the CHEBOP polynomial eigenvalue
% problem of degree p:
%    (A0 + lambda*A1 + ... + lambda^p*Ap)*x = 0.
% The input is p+1 CHEBOPs, A0, A1, ..., Ap and the output is an inf-by-K
% chebfun quasimatrix, X, whose columns are the K least oscillatory
% eigenfunctions, and a vector of length K, E, whose elements are the
% eigenvalues.
%    for j = 1:K
%       lambda = E(j)
%       u = X(:,j)
%       A0(u) + lambda*A1(u) + ... + lambda^p*Ap(u) %is approximately 0.
%    end
%
% E = POLYEIGS(A0,A1,..,Ap,K) is a vector of length k whose elements are
% the K least oscillatory eigenvalues of the polynomial eigenvalue problem.
%
% POLYEIGS(A0,A1,..,Ap,K,SIGMA) also finds K solutions to the polynomial
% eigenvalue problem. If SIGMA is a scalar, the eigenvalues found are the
% ones closest to SIGMA. Other possibilities are 'LR' and 'SR' for the
% eigenvalues of largest and smallest real part, and 'LM' (or Inf) and 'SM'
% for largest and smallest magnitude. SIGMA must be chosen appropriately
% for the given operator; for example, 'LM' for an unbounded operator will
% fail to converge!
%
% Similarly to CHEBOP/EIGS, this routine uses the built-in POLYEIG on dense
% matrices of increasing size, stopping when the targeted eigenfunctions
% appear to have converged, as determined by the chebfun constructor.
%
% Example:
%   A = chebop(@(x,u) diff(u,2),[-1 1],'dirichlet');
%   B = chebop(@(x,u) -x.*diff(u));
%   C = chebop(@(x,u) u);
%   [V D] = polyeigs(A,B,C,6,0)
%   plot(V)
%
% See also CHEBOP/EIGS, POLYEIG, LINOP/POLYEIGS.

% Copyright 2022 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Did we get preferences passed?
if ( isa(varargin{end}, 'cheboppref') )
    prefs = varargin{end};
    isPrefGiven = 1;
else
    prefs = cheboppref();
    isPrefGiven = 0;
end

% Tell CHEBOP/LINEARIZE() to stop if it detects nonlinearity:
linCheck = true; 

% Linearize, thereby obtaining linearity information, a LINOP, and an input of
% the correct dimensions to pass to N:
[L, ~, isLinear, u0] = linearize(N, N.init, [], linCheck);

% We need the entire operator (including BCs) to be linear:
assert(all(isLinear), 'CHEBFUN:CHEBOP:eigs:nonlinear', ...
    ['The input operator appears to be nonlinear.\n', ...
    'EIGS() supports only linear CHEBOP instances.']);

% Support for generalised problems:
for k = 1:nargin-1
    if ( isa(varargin{k}, 'chebop') )
        % Tell CHEBOP/LINEARIZE() that we don't want it to try to reshape inputs
        % that it believes are parameters to doubles, rather than CHEBFUNs.
        paramReshape = false;
        % Linearise the second CHEBOP:
        [varargin{k}, ~, isLinear] = ...
            linearize(varargin{k}, u0, [], linCheck, paramReshape);
    
        % We need the entire operator (including BCs) to be linear:
        assert(all(isLinear), 'CHEBFUN:CHEBOP:eigs:nonlinear', ...
            ['The second input operator appears to be nonlinear.\n', ...
            'EIGS() supports only linear CHEBOP instances.']);
        
    end
end

% Determine the discretization.
prefs = determineDiscretization(N, length(L.domain), prefs);

% Clear boundary conditions if the dicretization uses periodic functions (since
% if we're using periodic basis functions, the boundary conditions will be
% satisfied by construction).
disc = prefs.discretization();
tech = disc.returnTech();
if ( isPeriodicTech(tech()) )
    [~, L] = clearPeriodicBCs(N, L);
end

% Add the preferences in vargarin to pass them to LINOP/EIGS.
if ( isPrefGiven )
    % If a CHEBOPPREF was passed to the method, it will have been at the last
    % position of varargin, indexed at nargin-1. Overwrite it with the current
    % PREFS, as the discretization might have changed in the periodic case:
    varargin{nargin-1} = prefs;
else
    % Otherwise, add the PREFS to VARARGIN, so that it can be passed to the call
    % to LINOP/EIGS below.
    varargin{nargin} = prefs;
end

% Call LINOP/EIGS.
[varargout{1:nargout}] = polyeigs(L, varargin{:});

% Return a CHEBFUN rather than a CHEBMATRIX for scalar problems:
if ( nargout > 1 && isa(varargout{1}, 'chebmatrix') )
    u = varargout{1};
    if ( all(size(u) == [1 1]) )
        u = u{1};
    end
    varargout{1} = u;
end

end
