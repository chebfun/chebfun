function D = diffmat(N, varargin)
%DIFFMAT   Differentiation matrix.
%   D = DIFFMAT(N) returns the NxN differentiation matrix associated with the
%   Chebyshev spectral collocation method at second-kind Chebyshev points. 
%
%   D = DIFFMAT(N, P) returns the differentiation matrix of order P. See
%   COLLOC2.DIFFMAT for further details.
%
%   D = DIFFMAT(N, P, DOM) scales the differentiation matrix D to the domain
%   DOM. DOM should be a 1x2 vector.
%
%   D = DIFFMAT(N, P, DOM, DISC) or DIFF(N, P, DISC) returns the differentiation
%   matrix associated with the CHEBDISCRETIZATION DISC.
%
%   D = DIFFMAT(N, 'periodic') returns the N x N first-order Fourier 
%   differentiation matrix 1st-order Fourier diff mat.
%
%   D = DIFFMAT(N, P, 'periodic') returns the N x N Fourier differentiation 
%   matrix of order P.
%
%   D = DIFFMAT(N, P, 'periodic', DOM) scales the Pth-order Fourier 
%   differetiation matrix to the domain DOM.
%
%   D = DIFFMAT(N, P, 'rect') or D = DIFFMAT([N -P], P) returns (N-P) x N 
%   rectangular differentiation matrix of order P which maps from an N-point 
%   Chebyshev grid of second kind to an (N-P)-point Chebyshev grid of first kind.
%
%   D = DIFFMAT([N -P], P, DOM) returns (N-P) x N rectangular differentiation 
%   matrix of order P which is scaled to the domain DOM.
%
%   D = DIFFMAT([N -P], P, DOM, DISC) returns the Pth-order(N-P) x N rectangular 
%   differentiation matrix on domain DOM associated with the CHEBDISCRETIZATION 
%   DISC. When DISC is specified as 'colloc1' or colloc1() or @colloc1, D maps
%   between Chebyshev grids of first kind of size N and (N-P). when DISC is set
%   'colloc2' or colloc2() or @colloc2, D is same as default and maps from a 
%   second-kind Chebyshev grid.
%
%   D = DIFFMAT([M N]) returns an M x N first-order rectangular differentiation 
%   matrix which maps from an N-point Chebyshev grid of second kind to an 
%   M-point Chebyshev grid of first kind.
%   
%   D = DIFFMAT([M N], P) returns an M x N rectangular differentiation matrix of 
%   order P which maps from an N-point Chebyshev grid of second kind to an 
%   M-point Chebyshev grid of first kind.
%
%   D = DIFFMAT([M N], P, DOM) returns the same D but scaled to the domain DOM.
%
%   D = DIFFMAT([M N], P, DOM, DISC) returns a rectangular differentiation 
%   matrix associated with the CHEBDISCRETIZATION DISC. D maps from the N-point
%   Chebyshev grid of first kind when DISC is specified as 'colloc1' or 
%   colloc1() or @colloc1. It is same as default when DISC is 'colloc2' or 
%   colloc2() or @colloc2. 
%
%   D = DIFFMAT(N, P, DOM, DISC, LBC, RBC) or D = DIFFMAT([N -P], P, DOM, DISC, 
%   LBC, RBC) or D = DIFFMAT([N-P N], P, DOM, DISC, LBC, RBC) returns a square
%   differentiation matrix which is deflated by including information about the
%   boundary conditions specified by the cell structures LBC and RBC. When the
%   original differentiation matrix (without boundary condition information) is 
%   square, the boundary conditions are included by row replacement, whereas 
%   this is done by row appending when the original differentiation matrix is 
%   rectangular. The string specifier 'dirichlet' and 'neumann' indicate
%   Dirichlet and Neumann boundary conditions respectively. The string 'sum' 
%   indicates a side condition of definite integral i.e., sum(U), where U is the
%   solution to the system formed by D.
%
%   Example 1: D = DIFFMAT(N, P, DOM, DISC, {'dirichlet'}, {'neumann'}) replaces
%   the first and the last row of a square differentiation matrix by Dirichlet 
%   and Neumann boundary conditions respectively.
%
%   Example 2: D = DIFFMAT([N -P], P, DOM, DISC, {}, {'neumann' 'sum'}) appends 
%   two rows to an (N-P) x N rectangular differentiation matrix to square it up.
%   The first appended row corresponds to Neumann boundary condition at the 
%   right endpoint while the second appended row corresponds to a side condition 
%   sum(U) = I for a scalar I, where U is the solution to the resulting system.
%   
%   Note that for the case of row appending, the number of the boundary 
%   conditions given by LBC and RBC must total P. Also note that RBC can be 
%   omitted if there are only left boundary conditions.
%
% See also DIFF, COLLOC2.DIFFMAT, CUMSUMMAT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
%   Xu, K. and Hale, N., Explicit construction of rectangular differentiation
%   matrices, submitted 2014.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse the inputs:
[m, n, p, dom, bc, nlbc, nrbc, disc] = parseInputs(N, varargin{:});

%% Different cases:
if ( m == n ) % Square case:
    D = disc.diffmat(n, p);
else
    if ( p == 1 )
        if ( isa(disc, 'colloc1') )
            D = rectdiff1(m, n);
        elseif ( isa(disc, 'colloc2') )
            D = rectdiff2(m, n);
        end
    else
        if ( isa(disc, 'colloc1') )
            D = rectdiff_rec(m, n, p, 1);
        elseif ( isa(disc, 'colloc2') )
            D = rectdiff_rec(m, n, p, 2);
        end
    end
end

%% Rescaling:
scl = (2/(dom(end) - dom(1)))^p;
D = scl*D;

%% Boundary conditions:
if ( ~isempty(bc) )
    nbc = nlbc + nrbc;
    BC = zeros(nbc, n);
    for j = 1:nbc
        if ( j <= nlbc )
            y = -1;
            r = pi;
            idx = 1;
        else
            y = 1;
            r = 0;
            idx = n;
        end
        
        switch bc{j}
            case 'dirichlet'
                if ( isa(disc, 'colloc1') )
                    [x, ignored, v, t] = chebpts(n, 1);
                    BC(j,:) = barymat(y, x, v, r, t);
                else
                    I = eye(n);
                    BC(j,:) = I(idx,:);
                end
                
            case 'neumann'
                
                if ( isa(disc, 'colloc1') )
                    DD = diffmat(n, 1, dom, 'colloc1');
                    [x, ignored, v, t] = chebpts(n, 1);
                    P = barymat(y, x, v, r, t);
                    DD = P*DD;
                    BC(j,:) = DD;
                else
                    DD = diffmat(n, 1, dom);
                    BC(j,:) = DD(idx,:);
                end
                
            case 'sum'
                
                if ( isa(disc, 'colloc1') )
                    [ignored, w] = chebpts(n, dom, 1);
                else
                    [ignored, w] = chebpts(n, dom);
                end
                
                BC(j,:) = w;
                
            otherwise
                error('CHEBFUN:diffmat:wrongBC', ...
                    'Unknown type of boundary conditions.');
        end
    end
end

%% Replace or append the boundary conditions:
if ( ~isempty(bc) && ( m == n ) )
    % Replacement for square case:
    D(1:nlbc, :) = BC(1:nlbc, :);
    D(end-nrbc+1:end,:) = BC(nlbc+1:end,:);
elseif ( ~isempty(bc) )
    D = [BC(1:nlbc,:); D; BC(nlbc+1:end,:)];
end

end


function D = rectdiff1(m, n)
%RECTDIFF1  Explicit constrcution of 1st-order rectangular differentiation  
%matrix mapping from 1st-kind grid.

% mapping-from grid (angles):
T = chebtech1.angles(n).';

% difference between dimensions:
c = n-m;

% mapping-to grid (angles):
TAU = chebtech1.angles(m);

% Denominator:
denom = bsxfun(@(u,v) 2*sin((v+u)/2).*sin((v-u)/2), T, TAU);

% Sign:
sgn = ones(m, n);
if ( mod(c, 2) )
    sgn(1:2:end,2:2:end) = -1;
    sgn(2:2:end,1:2:end) = -1;
else
    sgn(1:2:end,1:2:end) = -1;
    sgn(2:2:end,2:2:end) = -1;
end

D = sgn.*((cos(c*TAU)./sin(TAU))*sin(T)-sin(c*TAU)/n*sin(T)./denom)./denom;

% indices for applying negative-sum trick: 
[~, idx] = min(abs(denom), [], 2);
idx = sub2ind([m n], 1:m, idx.');
D(idx) = 0;
D(idx) = -sum(D, 2);

if ( c == 1 )
    % Flipping trick:
    ii = logical(rot90(tril(ones(m, n)), 2));
    rot90D = rot90(D,2);
    D(ii) = -rot90D(ii);
end

end

function D = rectdiff2(m, n)
%RECTDIFF2  Explicit constrcution of 1st-order rectangular differentiation   
%matrix mapping from a 2nd-kind grid.

nm1 = n - 1;                    % For convenience.
cm1 = nm1 - m;                  % Difference between dimensions:
t = chebpts(n).';               % Second-kind grid.
tau = chebpts(m, 1);            % First-kind grid.
T = chebtech2.angles(n).';          % Second-kind grid (angles).
TAU = chebtech1.angles(m);   % First-kind grid (angles).

% Explicit expression:
denom = 2*bsxfun( @(u,v) sin((v+u)/2) .* sin((v-u)/2), T, TAU );
numer = bsxfun( @times, 1 - tau*t, cos(cm1*TAU)./sin(TAU) );

sgn = (-1)^cm1;

if ( cm1 == 0 )
    D = numer ./ denom.^2 / nm1;
else
    D = repmat(sin(cm1*TAU), 1, n)./denom + numer./denom.^2 / nm1;
    D = sgn*D;
end
D(:,[1,n]) = .5*D(:,[1,n]);     % Scaling for first and last columns.

if ( cm1 == 0 )
    % Flipping trick:
    ii = logical(rot90(tril(ones(m, n)), 2));
    rot90D = rot90(D,2);
    D(ii) = sgn*rot90D(ii);
end

% Sign:
D(1:2:end,1:2:end) = -D(1:2:end,1:2:end);
D(2:2:end,2:2:end) = -D(2:2:end,2:2:end);

% Negative sum trick:
[~, idx] = min(abs(denom), [], 2); 
idx = sub2ind([m n], 1:m, idx.');
D(idx) = 0; D(idx) = -sum(D, 2);

if ( cm1 == 0 )
    % Fix corner values:
    D(1) = -.25/(nm1*sin(pi/(2*m))*sin(pi/(4*m))^2); D(end) = -D(1);
    % Negative sum trick for corner entries:
    D(1,2) = -sum(D(1,[1 3:end])); D(end,end-1) = -D(1,2);
end

end
    
function D = rectdiff_rec(m, n, p, kind)
%Recursive construction for high-order rectangular differentiation matrices

% Sign and scaling:
sgn = ones(1, n);
sgn(1:2:end) = -1;

if ( kind == 1 )
    % mapped-from angles (1st-kind):
    T = chebtech1.angles(n).';
    
    % Compute the first order diff matrix:
    D = rectdiff1(m, n);
    
    % Preparation for higher order (p>1):
    a = [1; zeros(n,1)];
    
    % Signs:
    sgn = (-1)^(n-1)*sgn.*sin(T)/n;
    
else
    % mapped-from grid (2nd-kind):
    T = chebtech2.angles(n).';
    
    % Compute the first order diff matrix:
    D = rectdiff2(m, n);
    
    % Preparation for higher order (p>1):
    a = [1; 0; -1; zeros(n-2,1)];
    
    % Signs:
    sgn = (-1)^(n-1)*sgn/(2*(n-1));
    sgn([1 n]) = sgn([1 n])/2;
    
end

% mapped-to grid (1st-kind):
tau = chebpts(m, 1);
TAU = chebtech1.angles(m);

a = computeDerCoeffs(a);

% Denominator:
denom = bsxfun(@(u,v) 2*sin((v+u)/2).*sin((v-u)/2), T, TAU);

% indices for applying negative-sum trick:
[~, idx] = min(abs(denom), [], 2);
idx = sub2ind([m n], 1:m, idx.');

% Recursion for higher-order matrices:
for l = 2:p
    
    % Compute coefficients of the derivative of T_n:
    a = computeDerCoeffs(a);
    
    % Evaluating at tau by Clenshaw method:
    Tt = chebtech.clenshaw(tau, a);
    D = (Tt*sgn + l*D)./denom;
    
end

% negative-sum trick:
D(idx) = 0;
D(idx) = -sum(D, 2);

end

function cout = computeDerCoeffs(c)
%COMPUTEDERCOEFFS   Recurrence relation for coefficients of derivative.
%   C is the matrix of Chebyshev coefficients of a (possibly array-valued)
%   CHEBTECH object.  COUT is the matrix of coefficients for a CHEBTECH object
%   whose columns are the derivatives of those of the original.
    
    [n, m] = size(c);
    cout = zeros(n-1, m);                     % Initialize vector {c_r}
    w = repmat(2*(n-1:-1:1)', 1, m);
    v = w.*c(1:end-1,:);                      % Temporal vector
    cout(1:2:end,:) = cumsum(v(1:2:end,:));   % Compute c_{n-2}, c_{n-4},...
    cout(2:2:end,:) = cumsum(v(2:2:end,:));   % Compute c_{n-3}, c_{n-5},...
    cout(end,:) = .5*cout(end,:);             % Adjust the value for c_0
end

function [m, n, p, dom, bc, nlbc, nrbc, disc] = parseInputs(N, varargin)
% Parse the inputs to DIFFMAT.

p = 1;
dom = [-1 1];
disc = colloc2();
bc = [];
lbc = {};
nlbc = 0;
rbc = {};
nrbc = 0;
islbc = 1;

if ( isscalar(N) )
    n = N;
    m = n;
elseif ( N(2) < 0 )
    n = N(1);
    m = n + N(2);
else
    m = N(1);
    n = N(2);
end

if ( numel(varargin) == 0 )
    % Trivial case.
    return
end

for j = 1:numel(varargin)
    v = varargin{j};
    if ( isnumeric(v) )
        if ( isscalar(v) )
            p = v;
        else
            dom = v;
        end
    elseif ( ischar(v) && strcmpi(v, 'rect') )
        m = n - p;
    elseif ( isa(v, 'function_handle') || ischar(v) || ...
        isa(v, 'chebDiscretization') )
        disc = v;
    elseif ( iscell(v) )
        if ( islbc )
            lbc = v;
            islbc = 0;
        else
            rbc = v;
        end
    else
        error('CHEBFUN:diffmat:unknown', ...
            'Unknown input of type %s.', class(v));
    end
end

if ( p < 0 )
    error('CHEBFUN:diffmat:wrongInput', ...
            'The order of differentiation matrix must be non-negative.');
end

% Ensure DISC is a discretization:
if ( ischar(disc) )
    if ( strcmpi(disc, 'periodic') )
        disc = fourtech();
        if ( m ~= n )
            error('CHEBFUN:diffmat:wrongInput', ...
                'Rectangular Fourier differentiation matrices are not supported.');
        end
    else
        disc = str2func(disc);
    end
end
if ( isa(disc, 'function_handle') )
    disc = disc();
end

% No breakpoints allowed:
if ( numel(dom) > 2 )
    dom = dom([1 end]);
    warning('CHEBFUN:diffmat:noBreaks', ...
        'DIFFMAT does not support domains with breakpoints.');
end

% Boundary conditions:
bc = [lbc rbc];
if ( ~isempty(bc) )
    if ( isa(disc, 'fourtech')  )
        error('CHEBFUN:diffmat:wrongBC', ...
            ['For periodic functions, there is no need to specify boundary ' ...
             'conditions.']);
    end
    nlbc = numel(lbc);
    nrbc = numel(rbc);
    nbc = nlbc + nrbc;
    if ( nbc ~= p )
        error('CHEBFUN:diffmat:wrongBC', ...
            ['The number of boundary conditions must match differentiation ' ...
             'order p.']);
    end
end

end