function out = chebcoeffs(f, varargin)
%CHEBCOEFFS   Chebyshev coefficients of a CHEBFUN.
%   A = CHEBCOEFFS(F) returns the Chebyshev coefficients of F assuming 
%   it is a global (not piecewise) chebfun.  This is column vector of
%   coefficients such that  F = A(1) T_0(x) + ... + A(N) T_(N-1)(x),
%   where N is the length of F.
%
%   If F is an array-valued chebfun, then A is a matrix with the
%   same number of columns as F.
%
%   If the domain of F is [a,b] rather than [-1,1], then the 
%   coefficients are those of F transplanted to [-1,1].
%
%   If F is a piecewise chebfun, you can extract the Chebyshev
%   coefficients of the pieces with GET(F, 'COEFFS').
%
%   Alternatively, A = CHEBCOEFFS(F, N) returns the first N Chebyshev
%   coefficients of a piecewise chebfun F even though F is not
%   represented by a global Chebyshev expansion.  Chebfun does this
%   by evaluating appropriate integrals. 
%
%   A = CHEBCOEFFS(F, 'kind', 2) or A = CHEBCOEFFS(F, N, 'kind', 2)
%   return vectors or matrices corresponding to expansions 
%   F = A(1) U_0(x) + ... + A(N) U_(N-1)(x) in Chebyshev polynomials
%   of the second kind.
%
%   Examples:
%    x = chebfun('x');
%    chebcoeffs(exp(x))
%    chebcoeffs(abs(x),10)
%
% See also LEGCOEFFS, TRIGCOEFFS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
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

%% Merge:
% Try to merge out breakpoints if possible, since integrals are expensive:
if ( numFuns ~= 1 )
    f = merge(f);
    numFuns = numel(f.funs);
end

%% Error checking:
if ( isempty(N) && numFuns > 1 )
    error('CHEBFUN:CHEBFUN:chebcoeffs:inputN', ...
        'Input N is required for piecewise CHEBFUN objects.');
end
if ( ~isempty(N) && (~isscalar(N) || isnan(N)) )
    error('CHEBFUN:CHEBFUN:chebcoeffs:inputN', 'Input N must be a scalar.');
end
if ( any(isinf(f.domain)) )
% Chebyshev coefficients are not defined on an unbounded domain.
    error('CHEBFUN:CHEBFUN:chebcoeffs:infint', ...
        'Infinite intervals are not supported here.');
end

%% Compute the coefficients:
if ( numFuns == 1 )
    % CHEBCOEFFS() of a smooth piece:
    out = chebcoeffs(f.funs{1}, N, kind);

else
    % CHEBCOEFFS() of a piecewise smooth CHEBFUN:
    % (Compute coefficients via inner products.)

    % Construct a Chebfun of the appropriate Chebyshev weight:
    d = f.domain;
    x = chebfun([d(1) ; d(end)], [d(1), d(end)]);
    w = sqrt((x - d(1)).*(d(end) - x));
    if ( kind == 1 )
        w = 1./w;
    end
    
    % Chebyshev polynomials up to degree N - 1:
    T = chebpoly(0:(N-1), d, kind);
    
    % Compute the weighted inner products:
    numCols = numColumns(f);
    out = zeros(N, numCols);
    for j = 1:numCols
        fjw = extractColumns(f, j).*w;
        out(:,j) = innerProduct(fjw, T);
    end
     
    if ( kind == 1 )
        % Scale the T_0 term:
        out(1,:) = out(1,:)/2;
    end
    
    % Scale by 2/pi:
    out = (2/pi)*out;
    
end

end
