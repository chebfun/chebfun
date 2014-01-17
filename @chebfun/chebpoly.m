function out = chebpoly(f, varargin)
%CHEBPOLY   Chebyshev polynomial coefficients of a CHEBFUN.
%   A = CHEBPOLY(F) returns the row vector of coefficients such that F_1 = A(1)
%   T_M(x) + ... + A(M) T_1(x) + A(M+1) T_0(x), where T_M(x) denotes the M-th
%   Chebyshev polynomial and F_1 denotes the first FUN of CHEBFUN F.
%
%   A = CHEBPOLY(F, I) returns the coefficients for the I-th FUN.
%
%   A = CHEBPOLY(F, I, N) truncates or pads the vector A so that N coefficients
%   of the I-th FUN of F are returned. However, if I is 0 then the global
%   coefficients of the whole CHEBFUN F are returned (by computing relevant
%   inner products with Chebyshev polynomials).
%
%   If F is array-valued with M columns, then A is an MxN matrix.
%
%   C = CHEBPOLY(F, ..., 'kind', 2) returns the vector of coefficients for the
%   Chebyshev expansion of F in 2nd-kind Chebyshev polynomials F_1 = C(1) U_M(x)
%   + ... + C(M) U_1(x) + C(M+1) U_0(x).
%
%   There is also a CHEBPOLY command in the Chebfun trunk directory, which
%   computes the CHEBFUN corresponding to the Chebyshev polynomial T_n.
%
% See also LEGPOLY.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial empty case:
if ( isempty(f) )
    out = [];
    return
end

%% Initialise:
argin = {}; 
ii = []; 
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
    ii = argin{1}; 
end
if ( numel(argin) > 1 ) 
    N = argin{2}; 
end

%% Error checking:
if ( isempty(ii) )
    if ( numFuns > 1 )
        warning('CHEBFUN:chebpoly:nfuns1',['Chebfun has more than one fun. ' ...
                'Returning Chebyshev coefficients for the first fun only. ' ...
                'Use CHEBPOLY(F, 1) to suppress this warning.']);
    end
    ii = 1; 
    % [TODO]: Perhaps ii = 0 should be the default behaviour?  (If we do this,
    % we need to pick a sane value for N.)
end
if ( ii > numFuns )
    error('CHEBFUN:chebpoly:nfuns2', 'Chebfun only has %s FUNSs.', ...
        num2str(numFuns));
end
if ( (numel(ii) > 1) || (numel(N) > 1) )
    error('CHEBFUN:chebpoly:scalar', 'Inputs I and N must be scalars.');
end
if ( (ii == 0) && isempty(N) )
    error('CHEBFUN:chebpoly:inputs', 'Input N must not be empty if I is zero.');
end
if ( ~isempty(N) && ~isnumeric(N) )
    error('CHEBFUN:chebpoly:inputN', 'Input N must be a scalar.');
end
if ( any(isinf(f.domain)) )
% Chebyshev coefficients are not defined on an unbounded domain.
    error('CHEBFUN:chebpoly:infint', ...
        'Infinite intervals are not supported here.');
end

%% Compute the coefficients:
if ( ii > 0 )
    % CHEBPOLY() of a smooth piece:
    out = chebpoly(f.funs{ii}, N).';
    
elseif ( (ii == 0) && (numFuns == 1))
    
    % CHEBPOLY() of a smooth piece:
    out = chebpoly(f.funs{1}, N).';    
    
else
    % CHEBPOLY() of a piecewise smooth CHEBFUN:
    % [TODO]: Add a test for this code.

    % Compute coefficients via inner products.
    d = f.domain([1, end]);
    x = chebfun('x', d);
    w = 1./sqrt((x-d(1)).*(d(2)-x));
    numCols = min(size(f));
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
