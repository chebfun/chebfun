function u = mldivide(N, f, varargin )
% Solve partial differential equation using a tensor products.
% The only ACA stuff going on here is in the determining the discretisation
% size required to approximate the solution.

% Try and make a chebfun2 out of right hand side.
if isa(f,'chebfun')
    f = @(x,y) f(x) + 0*y;
    warning('Chebop2:mldivide:rhs','Univariate righthand side.');
end
if isa(f,'double')
    f = @(x,y) f + 0*x;
end

%%% Go get stuff we will need: 
rect = N.domain; 
f = chebfun2(f, rect);
prefs = chebfunpref();
tol = max( prefs.cheb2Prefs.eps, 1e-14 );      % Be gentle!    
maxDiscretise = 2*prefs.cheb2Prefs.maxRank;    % Use maxRank to get max disc. 
minsample = 9;

%%% Find out what grid to start on and which directions to do adaptivity. %%%
if ( nargin == 3 && isa(varargin{1}, 'double') )
    % Nonadaptive solve with an m-by-m discretization
    m = varargin{1};
    n = m;
    adaptive_x = 0;
    adaptive_y = 0;
elseif ( nargin == 4 && isa(varargin{1},'double') && ~isinf(varargin{1}) ...
        && isa(varargin{2},'double') && ~isinf(varargin{2}) )
    % Nonadaptive solve with an m-by-n discretization
    m = varargin{1};
    n = varargin{2};
    adaptive_x = 0;
    adaptive_y = 0;
elseif ( nargin == 4 && isa(varargin{1},'double') && ~isinf(varargin{1})...
        && isinf(varargin{2}) )
    % Adaptive solve in the horizontal variable, nonadaptive in the other.
    % (An m-by-inf discretization.)
    m = varargin{1};
    n = minsample;
    adaptive_x = 1;
    adaptive_y = 0;
elseif ( nargin == 4 && isa(varargin{2},'double') && ~isinf(varargin{2})...
        && isinf(varargin{1}) )
    % Adaptive solve in the vertical variable, nonadaptive in the other.
    % (An inf-by-n discretization.)
    n = varargin{2};
    m = minsample;
    adaptive_x = 0;
    adaptive_y = 1;
elseif ( nargin == 2 || (nargin == 4 && isinf(varargin{2}) && isinf(varargin{1})))
    % This code does the solver with adaptive calls in both the x- and
    % y-direction.
    n = minsample;
    m = minsample;
    adaptive_x = 1;
    adaptive_y = 1;
else
    error('CHEBOP2:MLDIVIDE:SYNTAX','Unrecognized input syntax')
end

% Check the discreizations make sense: 
if ( abs( round(n) - n ) > 0 || abs( round(m) - m ) > 0)
    error('CHEBOP2:MLDIVIDE:DISC','Discretization size should be an integer');
end

%%% Adaptive solver %%%
Resolved_x = 0; Resolved_y = 0; Resolved = Resolved_x & Resolved_y;
while ( ( ~Resolved ) && ( max(m, n) < maxDiscretise ))
    % Solve PDE, return an m x n matrix of coefficients:
    X = chebop2.denseSolve(N, f, m, n);
    
    if ( adaptive_y ) 
        % Resolved in y-direction?
        [Resolved_y, m] = ResolveCheck( X, m, tol );
    end
    
    if ( adaptive_x ) 
        % Resolved in x-direction?
        [Resolved_x, n] = ResolveCheck( X.', n, tol );
    end
    
    % update tolerances:
    tol = updateTolerance(tol, m, n);
    
    % Resolved in both directions:
    Resolved = Resolved_x & Resolved_y;
    
    % Check we do not have nans:
    if ( any( isnan( X ) ) | any( isinf( X ) ) )
        error('CHEBOP2:MLDIVIDE:NANINF', 'Nonunique solution to PDE')
    end
end

% Did we stop without resolving:
if ( max(m,n) >= maxDiscretise )
    warning('CHEBOP2:MAXDISC',...
        'Maximum discretization reached. Solution may not accurate.')
end

% Truncate at the end:
X = cutTrailingCoefficients( X );

% Form a chebfun2 object to represent the solution:
u = chebfun2( rot90( X, 2 ), 'coeffs', rect);

end

function [Resolved, newDisc] = ResolveCheck( coeffs, oldDisc, tol )
% Basic Resolution check:

tail = max(abs(coeffs(:,end-8:end)));
Resolved = all( ( tail < 20*oldDisc*tol ) );

% If unresolved, then increase the grid:
if ( ~Resolved )
    newDisc = 2^(floor(log2(oldDisc))+1)+1;
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


function X = cutTrailingCoefficients( X )
%TODO:
%     % cut trailing coefficients
%     idy = find( max(abs(X))/max(abs(X(:))) > eps, 1, 'last');
%     idx = find( max(abs(X),[],2)/max(abs(X(:))) > eps, 1, 'last');
%     if isempty(idx) && ~isempty(idy)
%         X = X(1:idy,:);
%     elseif ~isempty(idx) && isempty(idy)
%         X = X(:,1:idx);
%     else
%         X = X(1:idy,1:idx);
%     end
end