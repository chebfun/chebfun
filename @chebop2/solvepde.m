function u = solvepde(N, f, varargin)
%SOLVEPDE   Solve CHEBOP2 partial differential equation.
%    U = SOLVEPDE(N, F) is equivalent to U = N \ F and solves the linear partial
%    differential equation N(U) = F. 
% 
%    U = SOLVEPDE(N, F, M, N) solves the PDE with a discretization size of
%    M-by-N. 
% 
%    U = SOLVEPDE(N, F, M, INF) solves the PDE adaptively in the 2nd variable
%    and a discretization size of M in the 1st. 
% 
%    U = SOLVEPDE(N, F, INF, N) solves the PDE adaptively in the 1st variable
%    and a discretization size of N in the 2nd. 
%
%    U = SOLVEPDE(N, F, INF, INF) is equivalent to U = SOLVEPDE(N, F).
%
% For further details about the PDE solver, see: 
% A. Townsend and S. Olver, The automatic solution of partial differential
% equations using a global spectral method, in preparation, 2014.
% 
% Warning: This PDE solver is an experimental new feature. It has not been
% publicly advertised.  
        
% Copyright 2014 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Try and make a CHEBFUN2 out of right hand side.
if isa(f, 'chebfun')
    f = @(x,y) f(x) + 0*y;
    warning('CHEBFUN:CHEBOP2:solvepde:rhs', 'Univariate righthand side.');
end
if isa(f, 'double')
    f = @(x,y) f + 0*x;
end

% Go get stuff we will need: 
rect = N.domain; 
f = chebfun2(f, rect);
prefs = chebfunpref();
tol = max(prefs.techPrefs.eps, 1e-14); % Be gentle!    
maxDiscretise_x = 2*prefs.cheb2Prefs.maxRank; % Use maxRank to get maxDisc. 
maxDiscretise_y = maxDiscretise_x;
minsample = 9;

% Find out what grid to start on, and which directions to do adaptivity.
if ( nargin == 3 && isa(varargin{1}, 'double') )
    % Nonadaptive solve with an m-by-m discretization.
    m = varargin{1};
    n = m;
    adaptive_x = 0;
    adaptive_y = 0;
    maxDiscretise_x = n + 1; 
    maxDiscretise_y = m + 1; 
elseif ( nargin == 4 && isa(varargin{1},'double') && ~isinf(varargin{1}) ...
        && isa(varargin{2},'double') && ~isinf(varargin{2}) )
    % Nonadaptive solve with an m-by-n discretization.
    m = varargin{1};
    n = varargin{2};
    adaptive_x = 0;
    adaptive_y = 0;
    maxDiscretise_x = n + 1; 
    maxDiscretise_y = m + 1; 
elseif ( nargin == 4 && isa(varargin{1},'double') && ~isinf(varargin{1})...
        && isinf(varargin{2}) )
    % Adaptive solve in the horizontal variable, nonadaptive in the other.
    % (An m-by-inf discretization.)
    m = varargin{1};
    n = minsample;
    adaptive_x = 1;
    adaptive_y = 0;
    maxDiscretise_y = m + 1;
elseif ( nargin == 4 && isa(varargin{2},'double') && ~isinf(varargin{2})...
        && isinf(varargin{1}) )
    % Adaptive solve in the vertical variable, nonadaptive in the other.
    % (An inf-by-n discretization.)
    n = varargin{2};
    m = minsample;
    adaptive_x = 0;
    adaptive_y = 1;
    maxDiscretise_x = n + 1;
elseif ( nargin == 2 || ( nargin == 4 && isinf(varargin{2}) && isinf(varargin{1}) ) )
    % This code does the solver with adaptive calls in both the x- and
    % y-direction.
    n = minsample;
    m = minsample;
    adaptive_x = 1;
    adaptive_y = 1;
else
    error('CHEBFUN:CHEBOP2:solvepde:syntax', 'Unrecognized input syntax.')
end

% Check if the discreizations make sense: 
if ( abs(round(n) - n) > 0 || abs(round(m) - m) > 0)
    error('CHEBFUN:CHEBOP2:solvepde:discSize', ...
        'Discretization size should be an integer.');
end

% Adaptive solver.
Resolved_x = 0; 
Resolved_y = 0; 
Resolved = Resolved_x & Resolved_y;

while ( ( ~Resolved ) && ( m < maxDiscretise_y ) &&...
                         ( n < maxDiscretise_x ) )
                     
    % Solve PDE, return an m x n matrix of coefficients:
    X = chebop2.denseSolve(N, f, n, m);
    
    if ( adaptive_y ) 
        % Resolved in y-direction?
        [Resolved_y, m] = resolveCheck(X, m, tol);
    else
        Resolved_y = 1; 
    end
    
    if ( adaptive_x ) 
        % Resolved in x-direction?
        [Resolved_x, n] = resolveCheck(X.', n, tol);
    else 
        Resolved_x = 1; 
    end
    
    % Update tolerances:
    tol = updateTolerance(tol, m, n);
    
    % Resolved in both directions:
    Resolved = Resolved_x & Resolved_y;
    
    % Check we do not have NANs/INFs:
    if ( any( isnan(X(:)) | isinf(X(:)) ) )
        error('CHEBFUN:CHEBOP2:solvepde:nanInf', 'Nonunique solution to PDE.')
    end
    
end

% Did we stop without resolving:
if ( ( m >= maxDiscretise_y ) || ( n >= maxDiscretise_x ) )
    warning('CHEBFUN:CHEBOP2:solvepde:maxDisc', ...
        'Maximum discretization reached. Solution may not be accurate.')
end

% Truncate at the end:
X = cutTrailingCoefficients(X);

% Form a CHEBFUN2 object to represent the solution:
u = chebfun2(X, 'coeffs', rect); 

end

function [Resolved, newDisc] = resolveCheck(coeffs, oldDisc, tol)
% Basic Resolution check:

tail = max(abs(coeffs(:,end-8:end)));
Resolved = all(tail < 20*oldDisc*tol);

% If unresolved, then increase the grid:
if ( ~Resolved )
    newDisc = 2^(floor(log2(oldDisc))+1) + 1;
else
    newDisc = oldDisc;
end

end

function tol = updateTolerance(tol, m, n)
% Increase tolerance so that weak corner singularities cause less of a problem.

% Increase tolerance on large grids:
if ( max(m, n) > 250 )
    tol = max(tol, 1e-10);
end

end

function X = cutTrailingCoefficients(X)
%[TODO]: Do we want to do this? Bad effects for hyperbolic/parabolic PDEs.
%
%     % Cut trailing coefficients.
%     idy = find(max(abs(X))/max(abs(X(:))) > eps, 1, 'last');
%     idx = find(max(abs(X), [], 2)/max(abs(X(:))) > eps, 1, 'last');
%     if ( isempty(idx) && ~isempty(idy) )
%         X = X(1:idy,:);
%     elseif ( ~isempty(idx) && isempty(idy) )
%         X = X(:,1:idx);
%     else
%         X = X(1:idy,1:idx);
%     end
end
