function out = chebcoeffs(f, varargin)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a CHEBFUN.
%   A = CHEBCOEFFS(F, N) returns the first N Chebyshev coefficients of F, i.e.,
%   the row vector such that F = ... + A(1) T_N(x) + ... + A(N) T_1(x) +
%   A(N+1) T_0(x), where T_M(x) denotes the M-th Chebyshev polynomial.
%
%   If F is a smooth CHEBFUN (i.e., with no breakpoints), then CHEBCOEFFS(F) is
%   equivalent to CHEBCOEFFS(F, LENGTH(F)).
%
%   If F is array-valued with M columns, then A is an MxN matrix.
%
%   C = CHEBCOEFFS(F, N, 'kind', 2) returns the vector of coefficients for the
%   Chebyshev expansion of F in 2nd-kind Chebyshev polynomials F = ... + C(1)
%   U_N(x) + ... + C(N) U_1(x) + C(N+1) U_0(x).
%
% See also LEGCOEFFS, FOURCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial empty case:
if ( isempty(f) )
    out = [];
    return
end

if ( numel(f) > 1 )
    % TODO: Why not?
    error('CHEBFUN:CHEBFUN:chebcoeffs:quasia', ...
        'CHEBCOEFFS does not support quasimatrices.');
end

%% Initialise:
argin = {}; 
N = [];
numFuns = numel(f.funs);
kind = 1;

% Parse inputs:
while ( ~isempty(varargin) )
    if ( any(strcmpi(varargin{1}, {'chebkind', 'kind'})) )
        kind = varargin{2};
        varargin(1:2) = [];
    else
        argin = [argin, varargin{1}];
        varargin(1) = [];
    end
end
if ( numel(argin) > 0 ) 
    N = argin{1};
end

%% Error checking:
if ( isempty(N) && numFuns == 1 )
    N = length(f);
end
if ( isempty(N) )
    error('CHEBFUN:CHEBFUN:chebcoeffs:inputN', ...
        'Input N is required for piecewise CHEBFUN objects.');
end
if ( ~isscalar(N) || isnan(N) )
    error('CHEBFUN:CHEBFUN:chebcoeffs:inputN', 'Input N must be a scalar.');
end
if ( any(isinf(f.domain)) )
% Chebyshev coefficients are not defined on an unbounded domain.
    error('CHEBFUN:CHEBFUN:chebcoeffs:infint', ...
        'Infinite intervals are not supported here.');
end

if ( numFuns ~= 1 )
    f = merge(f);
    numFuns = numel(f.funs);
end

%% Compute the coefficients:
if ( numFuns == 1 )
    
    % CHEBCOEFFS() of a smooth piece:
    out = chebcoeffs(f.funs{1}, N).';    
    
else
    % CHEBCOEFFS() of a piecewise smooth CHEBFUN:

    % Compute coefficients via inner products.
    d = f.domain([1, end]);
    x = chebfun('x', d);
    w = 1./sqrt((x - d(1)).*(d(2) - x));
    numCols = numColumns(f);
    out = zeros(numCols, N);
    f = mat2cell(f);
    for j = 1:numCols
        for k = 1:N
            T = chebpoly(k-1, d);
            I = (f{j}.*T).*w;
            out(j, N-k+1) = 2*sum(I)/pi;
        end
    end
    out(:,N) = out(:,N)/2;
    
end

% Return 2nd-kind coefficients:
if ( (kind == 2) && (numel(out) > 1) )
    out(:,end) = 2*out(:,end);
    % Recurrence relation / conversion matrix:
    out = .5*[out(:,1:2), out(:,3:end) - out(:,1:end-2)];
end

end
