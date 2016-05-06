function out = jaccoeffs(f, n, a, b)
%JACCOEFFS   Jacobi polynomial coefficients of a CHEBFUN.
%   A = JACCOEFFS(F, N, ALPHA, BETA) returns the first N+1 coefficients in the
%   Jacobi series expansion of the CHEBFUN F, so that such that F approximately
%   equals A(1) J_0(x) + ... + A(N+1) J_N(x) where J_N(x) denotes the N-th
%   Jacobi polynomial with parameters ALPHA and BETA. A is a column vector.
%
%   If F is smooth (i.e., numel(f.funs) == 1), then A = JACCOEFFS(F, ALPHA,
%   BETA) will assume that N = length(F).
%
%   JACCOEFFS does not support quasimatrices.
%
% See also CHEBCOEFFS, LEGCOEFFS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(f) > 1 )
    error('CHEBFUN:CHEBFUN:jaccoeffs:quasi', ...
        'JACCOEFFS does not support quasimatrices.');
end

if ( (numel(f.funs) == 1) && (nargin < 4) )
    b = a;
    a = n;
    n = length(f);
elseif ( isempty(n) )
    n = length(f);
end

%%
% Special cases:
if ( a == 0 && b == 0 )
    out = legcoeffs(f, n);
    return
elseif ( a == -.5 && b == -.5 )
    out = chebcoeffs(f, n, 'kind', 1);
    nn = 0:(n-2);
    scl = repmat([1 ; cumprod((nn'+.5)./(nn'+1))], 1, size(out, 2));
    out = out./scl;
    return
elseif ( a == .5 && b == .5 )
    out = chebcoeffs(f, n, 'kind', 2);
    nn = 0:(n-2);
    scl = repmat([1 ; cumprod((nn'+1.5)./(nn'+1))], 1, size(out, 2));
    out = (1:n)'.*out./scl;
    return
end
    
%% 
if ( numel(f.funs) == 1 )
    
    % Compute Jacobi coefficients of the single fun.
    out = jaccoeffs(f.funs{1}, n, a, b);
   
else

    % Jacobi-Vandermonde matrix:
    Enorm = jacpoly(0:n-1, a, b, f.domain);

    % Compute the projection (with correct scaling):
    out = Enorm \ f;
    
end

end
                
