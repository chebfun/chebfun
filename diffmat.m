function D = diffmat(N, varargin)
%DIFFMAT   Differentiation matrix.
%   D = DIFFMAT(N) returns the NxN differentiation matrix associated with the
%   Chebyshev spectral collocation method at second-kind Chebyshev points. 
%
%   D = DIFFMAT(N, K) returns the differentiation matrix of order K. See
%   COLLOC2.DIFFMAT for further details.
%
%   D = DIFFMAT(N, K, DOM) scales the differentiation matrix D to the domain
%   DOM. DOM should be a 1x2 vector.
%
%   D = DIFF(N, K, DOM, DISC) or DIFF(N, K, DISC) returns the differentiation
%   matrix associated with the CHEBDISCRETIZATION DISC.
%
% See also DIFF, COLLOC2.DIFFMAT, CUMSUMMAT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Parse the inputs:
p = 1;
dom = [-1 1];
disc = colloc2();

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
    else
        error('CHEBFUN:diffmat:unknown', ...
            'Unknown input of type %s.', class(v));
    end
end

if ( p < 0 )
    error('CHEBFUN:diffmat:wrongInput', ...
            'The order of differentiation matrix must be non-negative');
end

% Ensure DISC is a discretization:
if ( ischar(disc) )
    if ( strcmpi(disc, 'periodic') )
        disc = fourtech();
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

end

function D = rectdiff1(m, n)
%RECTDIFF1  Explicit constrcution rectangular differentiation matrix mapping 
%from 1st-kind grid.

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

end

function D = rectdiff2(m, n)
%RECTDIFF2  Explicit constrcution rectangular differentiation matrix mapping 
%from a 2nd-kind grid.

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

% Flipping trick:
ii = logical(rot90(tril(ones(m, n)), 2));
rot90D = rot90(D,2);
D(ii) = sgn*rot90D(ii);

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
    
    sgn = (-1)^(n-1)*sgn.*sin(T)/n;
    
else
    % mapped-from grid (2nd-kind):
    T = chebtech2.angles(n).';
    
    % Compute the first order diff matrix:
    D = rectdiff2(m, n);
    
    % Preparation for higher order (p>1):
    a = [1; 0; -1; zeros(n-2,1)];
    
    sgn = (-1)^(n-1)*sgn/(2*(n-1));
    sgn([1 n]) = sgn([1 n])/2;
    
end

% mapped-to grid (1st-kind):
tau = chebpts(m, 1);
TAU = chebtech1.angles(m);

a = deri(a);

% Denominator:
denom = bsxfun(@(u,v) 2*sin((v+u)/2).*sin((v-u)/2), T, TAU);

% indices for applying negative-sum trick:
[~, idx] = min(abs(denom), [], 2);
idx = sub2ind([m n], 1:m, idx.');

% Recursion for higher-order matrices:
for l = 2:p
    
    % Compute coefficients of the derivative of T_n:
    a = deri(a);
    
    % Evaluating at tau by Clenshaw method:
    Tt = clen(tau, a);
    D = (Tt*sgn + l*D)./denom;
    
end

% negative-sum trick:
D(idx) = 0;
D(idx) = -sum(D, 2);

end