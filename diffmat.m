function D = diffmat(N, varargin)
%DIFFMAT   Spectral differentiation matrix.
%   D = DIFFMAT(N) returns the N x N differentiation matrix associated with the
%   Chebyshev spectral collocation method at second-kind Chebyshev points. 
%
%   D = DIFFMAT(N, P) returns the N x N differentiation matrix of order P.
%
%   D = DIFFMAT(N, P, DOM) scales the differentiation matrix D to the domain
%   DOM. DOM should be a 1x2 vector.
%
%   D = DIFFMAT(N, P, DOM, GRID) returns the square differentiation matrix on 
%   grid specified by GRID. The specifier GRID can be 'chebkind1' (first-kind 
%   Chebyshev grid) or 'chebkind2' (second-kind Chebyshev grid) or 'leg' 
%   (Legendre grid).
%
%   D = DIFFMAT(N, P, DOM, GRID, LBC, RBC) returns the square differentiation 
%   matrix on grid GRID with the corresponding row(s) replaced by boundary 
%   conditions specified by string specifiers LBC and RBC for the left and right 
%   boundary respectively. The specifier LBC and RBC can be any of 'dirichlet', 
%   'neumann', and 'sum' with 'dirichlet' and 'neumann' indicating Dirichlet and
%   Neumann boundary conditions respectively. The string 'sum' indicates a side 
%   condition of definite integral i.e., sum(U), where U is the implied solution 
%   to the system formed by D. For a differentiation matrix of order higher than
%   1, if there are multiple boundary conditions at one boundary, the 
%   corresponding specifiers need to be grouped in a cell using curly brackets. 
%   If there is no boundary condition at a boundary, then it needs to be 
%   indicated by either an empty array, i.e. [], or an empty cell, i.e. {}. Note
%   that the number of the boundary conditions given by LBC and RBC must total 
%   P, i.e. the order of differentiation, otherwise an error is thrown. Also 
%   note that RBC can be omitted if there are only left boundary conditions.
%
%   Example 1: D = DIFFMAT(N, 1, DOM, GRID, [], 'dirichlet') replaces the last 
%   row of an N x N differentiation matrix by a Dirichlet boundary condition.
%
%   Example 2: D = DIFFMAT(N, 1, DOM, GRID, 'dirichlet') replaces the first 
%   row of an N x N differentiation matrix by a Dirichlet boundary condition.
%
%   Example 3: D = DIFFMAT(N, 2, DOM, GRID, 'dirichlet', 'neumann') replaces
%   the first and last rows of an N x N differentiation matrix by a Dirichlet 
%   and a Neumann boundary conditions respectively.
%
%   Example 4: D = DIFFMAT(N, 3, DOM, GRID, 'dirichlet', {'dirichlet' 'neumann'}) 
%   replaces the first row of an N x N differentiation matrix by a Dirichlet
%   boundary condition and the last two rows by Dirichlet and Neumann boundary 
%   conditions.
%
%   D = DIFFMAT(N, 'periodic') returns the N x N first-order Fourier 
%   differentiation matrix on the default interval [-1 1]. The tag
%   'periodic' can be replaced by 'trig'.
%
%   D = DIFFMAT(N, P, 'periodic') returns the N x N Fourier differentiation 
%   matrix of order P  on the default interval [-1 1].
%
%   D = DIFFMAT(N, P, 'periodic', DOM) scales the Pth-order Fourier 
%   differetiation matrix to domain DOM.
%
%   The remaining options concern rectangular spectral differentiation matrices,
%   as introduced in [Driscoll & Hale 2014]. Chebfun uses these to solve ODE
%   boundary value problems in the mode GRID1 = 'chebkind2', GRID2 =
%   'chebkind1'. See also [Xu & Hale 2014].
%
%   D = DIFFMAT([M N]) returns the M x N first-order rectangular differentiation 
%   matrix which maps from an N-point Chebyshev grid of the second kind to an 
%   M-point Chebyshev grid of the same kind.
%   
%   D = DIFFMAT([M N], P) returns an M x N rectangular differentiation matrix of 
%   order P which maps from an N-point to an M-point Chebyshev grid, both of
%   second kind.
%
%   D = DIFFMAT([M N], P, DOM) returns the same D but scaled to the domain DOM.
%
%   D = DIFFMAT([M N], P, DOM, GRID) returns an M x N first-order rectangular 
%   differentiation matrix which maps from an N-point grid of type GRID to an 
%   M-point grid of the same type. The specifier GRID can be 'chebkind1', 
%   'chebkind2', or 'leg'.
%
%   D = DIFFMAT([M N], P, DOM, GRID1, GRID2) returns an M x N first-order 
%   rectangular differentiation matrix which maps from an N-point grid of type 
%   GRID1 to an M-point grid of type GRID2. The specifier GRID1 and GRID2 can be
%   any of 'chebkind1', 'chebkind2', and 'leg'.
%
%   D = DIFFMAT(N, P, DOM, GRID1, GRID2, LBC, RBC) returns a square 
%   differentiation matrix which is squared up and deflated by appending 
%   boundary conditions specified by string specifiers LBC and RBC for the left 
%   and the right boundary respectively. 
%
%   Example: D = DIFFMAT([M N], P, DOM, GRID1, GRID2, {}, {'neumann' 'sum'}) 
%   appends two rows to an M x N rectangular differentiation matrix to square it 
%   up. The first appended row corresponds to a Neumann boundary condition at 
%   the right endpoint while the second appended row corresponds to a side 
%   condition sum(U) = I for a scalar I, where U is the solution to the 
%   resulting system.
%
% See also DIFF, CHEBCOLLOC2.DIFFMAT, CUMSUMMAT, DIFFROW, INTMAT, INTROW.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References for rectangular differentiation matrices:
%
%   Driscoll, T. and Hale, N., Rectangular spectral collocation, submitted 2014.
%
%   Xu, K. and Hale, N., Explicit construction of rectangular differentiation
%   matrices, submitted 2014.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse the inputs:
[m, n, p, dom, bc, nlbc, nrbc, mapFrom, mapTo] = parseInputs(N, varargin{:});

%% Different cases:
if ( strcmpi(mapFrom, mapTo) && ( m == n ) ) % Square case:
    if ( strcmpi(mapFrom, 'chebkind1') )
        D = chebcolloc1.diffmat(n, p);
    elseif ( strcmpi(mapFrom, 'chebkind2') )
        D = chebcolloc2.diffmat(n, p);
    elseif ( strcmpi(mapFrom, 'periodic') )
        D = trigcolloc.diffmat(n, p);
    else
        [x, ignored, v] = legpts(n); %#ok<ASGLU>
        [y, ignored, w] = chebpts(n); %#ok<ASGLU>
        z = legpts(n);
        P1 = barymat(y, x, v);
        P2 = barymat(z, y, w);
        D = chebcolloc2.diffmat(n, p);
        D = P2*D*P1;
        
        % Flipping trick for symmetry:
        d = diag(rot90(D));
        d = sign(d).*(abs(d) + abs(flipud(d)))/2;
        D(logical(flipud(eye(N)))) = d;
        DRot = rot90(D, 2);
        idxTo = rot90(~triu(ones(N)));
        D(idxTo) = (-1)^p*DRot(idxTo);
        if ( mod(N, 2) == 1 )
            D((N+1)/2,(N+1)/2) = 0;
        end
    end
    
elseif ( strcmpi(mapTo, 'chebkind1') )
    
    if ( strcmpi(mapFrom, 'chebkind1') )
        if ( p == 1 )
            D = rectdiff1(m, n);
        else
            D = rectdiff_rec(m, n, p, 1);
        end
    elseif ( strcmpi(mapFrom, 'chebkind2') )
        if ( p == 1 )
            D = rectdiff2(m, n);
        else
            D = rectdiff_rec(m, n, p, 2);
        end
    else
        [x, ignored, v] = legpts(n);  %#ok<ASGLU>
        [y, ignored, w] = chebpts(n);  %#ok<ASGLU>
        z = chebpts(m, 1);
        P1 = barymat(y, x, v);
        P2 = barymat(z, y, w);
        D = chebcolloc2.diffmat(n, p);
        D = P2*D*P1;
    end
            
elseif ( strcmpi(mapTo, 'chebkind2') )
    [z, ignored, ignored, s] = chebpts(m);  %#ok<ASGLU>
    if ( strcmpi(mapFrom, 'chebkind1') )
        D = chebcolloc1.diffmat(n, p);
        [x, ignored, v, r] = chebpts(n, 1);  %#ok<ASGLU>
        P = barymat(z, x, v, s, r, 1);
        D = P*D;
    elseif ( strcmpi(mapFrom, 'chebkind2') )
        D = chebcolloc2.diffmat(n, p);
        [x, ignored, v, r] = chebpts(n);  %#ok<ASGLU>
        P = barymat(z, x, v, s, r, 1);
        D = P*D;
    else
        [x, ignored, v] = legpts(n);  %#ok<ASGLU>
        [y, ignored, w] = chebpts(n);  %#ok<ASGLU>
        P1 = barymat(y, x, v);
        P2 = barymat(z, y, w);
        D = chebcolloc2.diffmat(n, p);
        D = P2*D*P1;
    end
    
elseif ( strcmpi(mapTo, 'leg') )
    z = legpts(m);
    if ( strcmpi(mapFrom, 'chebkind1') )
        D = chebcolloc1.diffmat(n, p);
        [x, ignored, v] = chebpts(n, 1);  %#ok<ASGLU>
        P = barymat(z, x, v);
        D = P*D;
    elseif ( strcmpi(mapFrom, 'chebkind2') )
        D = chebcolloc2.diffmat(n, p);
        [x, ignored, v] = chebpts(n);  %#ok<ASGLU>
        P = barymat(z, x, v);
        D = P*D;
    else
        [x, ignored, v] = legpts(n);  %#ok<ASGLU>
        [y, ignored, w] = chebpts(n);  %#ok<ASGLU>
        P1 = barymat(y, x, v);
        P2 = barymat(z, y, w);
        D = chebcolloc2.diffmat(n, p);
        D = P2*D*P1;
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
            z = -1;
            r = pi;
            idx = 1;
        else
            z = 1;
            r = 0;
            idx = n;
        end
        
        switch bc{j}
            case 'dirichlet'
                if ( strcmpi(mapFrom, 'chebkind1') )
                    [x, ignored, v, t] = chebpts(n, 1);  %#ok<ASGLU>
                    BC(j,:) = barymat(z, x, v, r, t);
                elseif ( strcmpi(mapFrom, 'chebkind2') )
                    I = eye(n);
                    BC(j,:) = I(idx,:);
                else
                    [x, ignored, v] = legpts(n);  %#ok<ASGLU>
                    BC(j,:) = barymat(z, x, v);
                end
                
            case 'neumann'
                
                if ( strcmpi(mapFrom, 'chebkind1') )
                    DD = diffmat(n, 1, dom, 'chebkind1');
                    [x, ignored, v, t] = chebpts(n, 1);  %#ok<ASGLU>
                    P = barymat(z, x, v, r, t);
                    DD = P*DD;
                    BC(j,:) = DD;
                elseif ( strcmpi(mapFrom, 'chebkind2') )
                    DD = diffmat(n, 1, dom);
                    BC(j,:) = DD(idx,:);
                else
                    DD = diffmat(n, 1, dom);
                    [x, ignored, v] = legpts(n);  %#ok<ASGLU>
                    y = chebpts(n);
                    P = barymat(y, x, v);
                    DD = DD*P;
                    BC(j,:) = DD(idx, :);
                end
                
            case 'sum'
                
                if ( strcmpi(mapFrom, 'chebkind1') )
                    [ignored, w] = chebpts(n, dom, 1);  %#ok<ASGLU>
                elseif ( strcmpi(mapFrom, 'chebkind2') )
                    [ignored, w] = chebpts(n, dom);  %#ok<ASGLU>
                else
                    [ignored, w] = legpts(n, dom);  %#ok<ASGLU>
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
    a = [zeros(n,1); 1];
    
    % Signs:
    sgn = (-1)^(n-1)*sgn.*sin(T)/n;
    
else
    % mapped-from grid (2nd-kind):
    T = chebtech2.angles(n).';
    
    % Compute the first order diff matrix:
    D = rectdiff2(m, n);
    
    % Preparation for higher order (p>1):
    a = [zeros(n-2,1); -1; 0; 1];
    
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
    cout = zeros(n-1, m);                       % Initialize vector {c_r}
    w = repmat(2*(1:n-1)', 1, m);
    v = w.*c(2:end,:);                          % Temporal vector
    cout(n-1:-2:1,:) = cumsum(v(n-1:-2:1,:));   % Compute c_{n-2}, c_{n-4},...
    cout(n-2:-2:1,:) = cumsum(v(n-2:-2:1,:));   % Compute c_{n-3}, c_{n-5},...
    cout(1,:) = .5*cout(1,:);                   % Adjust the value for c_0
end

function [m, n, p, dom, bc, nlbc, nrbc, mapFrom, mapTo] = parseInputs(N, varargin)
% Parse the inputs to DIFFMAT.

p = 1;
dom = [-1 1];
mapFrom = [];
mapTo = [];
lbc = {};
nlbc = 0;
rbc = {};
nrbc = 0;
isLbcGiven = 0;
isRbcGiven = 0;

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

for j = 1:numel(varargin)
    v = varargin{j};
    if ( isnumeric(v) )
        if ( isempty(v) )
            if ( ~isLbcGiven )
                lbc = {};
                isLbcGiven = 1;
            elseif ( ~isRbcGiven )
                rbc = {};
                isRbcGiven = 1;
            end
        elseif ( isscalar(v) )
            p = v;
        else
            dom = v;
        end
    elseif ( ischar(v) )
        switch v
            case 'rect'
                
                if ( strcmpi(mapFrom, 'periodic') )
                    error('CHEBFUN:diffmat:wrongInput', ...
                        ['Rectangular Fourier differentiation matrices are '...
                        'not supported.']);
                end
                
                if ( isscalar(N) )
                    m = n - p;
                end
                
            case {'periodic', 'trig'}
                mapFrom = 'periodic';
                mapTo = 'periodic';
                if ( m ~= n )
                    error('CHEBFUN:diffmat:wrongInput', ...
                        ['Rectangular Fourier differentiation matrices are '...
                        'not supported.']);
                end
                
            case {'chebkind1', 'chebkind2', 'leg'}
                if ( isempty(mapFrom) )
                    mapFrom = v;
                elseif ( isempty(mapTo) )
                    mapTo = v;
                else
                    error('CHEBFUN:diffmat:unknown', ...
                        'Too many inputs for grid type.');
                end
                
            case {'dirichlet', 'neumann', 'sum', []}
                if ( ~isLbcGiven )
                    lbc = {v};
                    isLbcGiven = 1;
                elseif ( ~isRbcGiven )
                    rbc = {v};
                    isRbcGiven = 1;
                else
                    error('CHEBFUN:diffmat:unknown', ...
                        ['Too many inputs for boundary condition. ' ...
                        'Use curly brackets to group left and right ' ...
                        'boundary conditions, if multiple boundary ' ...
                        'conditions are considered at one boundary.']);
                end
            otherwise
                error('CHEBFUN:diffmat:unknown', ['Unknown input ', v]);
        end

    elseif ( iscell(v) )
        if ( ~isLbcGiven )
            lbc = v;
            isLbcGiven = 1;
        elseif ( ~isRbcGiven )
            rbc = v;
            isRbcGiven = 1;
        else
            error('CHEBFUN:diffmat:unknown', ...
                'Unrecognized boundary condition.');
        end
    else
        error('CHEBFUN:diffmat:unknown', ...
            'Unknown input of type %s.', class(v));
    end
end

% If only one grid is given, then let 
if ( ~isempty(mapFrom) && isempty(mapTo) )
    mapTo = mapFrom;
end

if ( isempty(mapFrom) )
    mapFrom = 'chebkind2';
    mapTo = 'chebkind2';
end   

if ( p < 0 )
    error('CHEBFUN:diffmat:wrongInput', ...
            'The order of differentiation matrix must be non-negative.');
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
    if ( isa(mapTo, 'periodic')  )
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
