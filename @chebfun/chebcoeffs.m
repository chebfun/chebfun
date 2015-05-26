function out = chebcoeffs(f, varargin)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a CHEBFUN.
%   A = CHEBCOEFFS(F, N) returns the first N Chebyshev coefficients of F is
%   the column vector A such that F = A(1) T_0(x) +  A(2) T_1(x) + ... + 
%   A(N) T_(N-1)(x), where T_M(x) denotes the M-th Chebyshev polynomial of the
%   first kind.
%
%   If F has a 'finite' Chebyshev expansion (i.e., it is a smooth CHEBFUN with
%   no breakpoints or endpoint singularities and is based on a CHEBTECH), then
%   CHEBCOEFFS(F) is equivalent to CHEBCOEFFS(F, LENGTH(F)). This syntax should
%   be used with caution, and passing N is preferred.
%
%   If F is array-valued with M columns, then A is an NxM matrix.
%
%   C = CHEBCOEFFS(F, N, 'kind', 2) returns the vector of coefficients of F
%   such that F = C(1) + C(2) U_1(x) + ... + C(N) U_(N-1)(x), where U_M(x)
%   denotes the M-th Chebyshev polynomial of the second kind.
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
    
    % We compute 2nd-kind coefficients by computing the 1st-kind coefficients and
    % using a recurrence relation.  The recurrence for the coefficient of U_n
    % requires the coefficients of T_n and T_{n + 2}, so we need to compute two
    % extra 1st-kind coefficients if 2nd-kind coefficients have been requested.
    if ( kind == 2 )
        N = N + 2;
    end
    
    % CHEBCOEFFS() of a smooth piece:
    out = chebcoeffs(f.funs{1}, N);   
    
    % Compute 2nd-kind coefficients from 1st-kind ones using the recurrence
    %   T_n(x) = (1/2)*(U_n(x) - U_{n-2}(x)):
    if ( kind == 2 )
        out(1,:) = 2*out(1,:);
        out = 0.5*[out(1:end-2,:) - out(3:end,:) ; out(end-1:end,:)];
        out = out(1:end-2,:);
    end
    
else
    % CHEBCOEFFS() of a piecewise smooth CHEBFUN:

    % Compute coefficients via inner products.
    d = f.domain([1, end]);
    x = chebfun('x', d);
    if ( kind == 1 )
        w = 1./sqrt((x - d(1)).*(d(2) - x));
    elseif ( kind == 2 )
        w = sqrt((x - d(1)).*(d(2) - x));
    end
    
    numCols = numColumns(f);
    out = zeros(N, numCols);
    f = mat2cell(f);
    T = chebpoly(0:(N-1), d, kind);
    T = mat2cell(T);
    for k = 1:N
        Tkw = T{k}.*w;
        for j = 1:numCols
            out(k,j) = innerProduct(f{j}, Tkw);
        end
    end
     
    if ( kind == 1 )
        % Scale the T_0 term:
        out(1,:) = out(1,:)/2;
    end
    
    % Scale by 2/pi:
    out = (2/pi)*out;
    
end

end
