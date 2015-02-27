function [p, s] = trigcf(f, n, M)
%CF   Trigonometric Caratheodory-Fejer approximation
%   P = CF(F, N) computes a degree N trigonometric CF approximant to 
%   CHEBFUN F defined on [a, b], which must consist of just a
%   single FUN.
%
%   [P, S] = CF(F, N) also returns S, the associated CF singular
%   value, an approximation to the minimax error.
%
%   [P, S] = CF(F, N, K) does the same but uses only the K-th
%   partial sum in the Fourier expansion of F.
%
%   TRIGCF = Caratheodory-Fejer approximation is a near-best approximation
%   that is often indistinguishable from the true best approximation (which
%   for tirgonometric polynomials can be computed using TRIGREMEZ(F, N),
%   but often much faster to compute.
%
%   Examples:
%
%   Compute a degree 2 trigonometric polynomial CF approximant to 
%   exp(sin(pi*x)) on [-1, 1]:
%     f = chebfun(@(x) exp(sin(pi*x)), 'trig');
%     [p, s] = trigcf(f, 2);
%
%   References:
%   [1] Javed, M. and Trefethen, L. N.  "Remez and CF approximations of 
%    Periodic Functions". In preparation.
%
% See also CF, TRIGREMEZ, REMEZ.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || nargin < 2 )
    p = f;
    s = 0;
    return;
end

% If M is not provided, set it to empty:
if ( nargin < 3 )
    M = [];
end

numCols = numColumns(f);
if ( numCols > 1 )
    if ( numel(f) == 1 )
        % Convert array-valued CHEBFUNs to quasimatrices first.
        f = cheb2quasi(f);
        [p, s] = trigcf(f, n, M);
    else
        % Deal with quasimatrices.
        s = zeros(1, numCols);
        for k = 1:numCols
            [pk, sk] = trigcfOneColumn(f(:, k), n, M);
            p(k) = pk;            
            s(k) = sk;
        end

        if ( f(1).isTransposed )
            p = p.';            
            s = s.';
        end
    end
else
    [p, s] = trigcfOneColumn(f, n, M);
end

end



%% 
% Main code for Trigonometric CF approximation:
function [p, s] = trigcfOneColumn(f, n, M)

% Check the inputs.
if ( any(isinf(domain(f))) )
    error('CHEBFUN:CHEBFUN:cf:unboundedDomain', ...
        'CF does not work for CHEBFUNs with unbounded domains.');
end

if ( issing(f) )
    error('CHEBFUN:trigcf:singularFunction', ...
        'TRIGCF does not support functions with singularities.');
end

% Deal with complex-valued functions.
if ( ~isreal(f) )
    warning('CHEBFUN:CHEBFUN:trigcf:complex', ...
        'TRIGCF does not work for complex valued functions. Taking real part.');
    f = real(f);
end

% Use the full expansion of f if M is not provided:
if ( isempty(M) || 2*M + 1 > length(f) )
    M = (length(f)-1)/2;
    if ( rem(M, 1) ~= 0 )
        error( 'M is not an integer, even length trigfun!' )
    end
end

% Trivial case: approximation degree exceeds that of the expansion degree.
if ( n >= M )
    p = f;
    % There is no error in this case:
    s = 0;
    return
end


% Map the problem to [-1, 1] and rescale:
dom = f.domain([1, end]);
normf = norm(f);
f = f/normf;
f = newDomain(f, [-1, 1]);

% If the degree of cf approximation is just one less than the expansion, we do
% not solve the eigenvalue problem:
if ( n == M-1 )
    a = trigcoeffs(f, 2*M+1);
    p = chebfun(a(2:end-1), 'coeffs', 'trig');
    s1 = a(1);
    s2 = a(1);
else
    %% Main Algorithm for TRIGCF:
    % Extract the Fourier coefficients to be used in computing the approximation.
    N = (length(f)-1)/2;
    if ( rem(N, 1) ~= 0 )
        error( 'N is not an integer, even length trigfun!' )
    end
    a = trigcoeffs(f, length(f));
    a = a(N+1:(N+1+M));
    
    ck = real(a(n+2:M+1));
    dk = imag(a(n+2:M+1));
    
    % Initialize arrays:
    b1 = zeros(2*n+1, 1);
    b2 = zeros(2*n+1, 1);
    s1 = 0;
    s2 = 0;
    
    % Solve the eigenvalue problem twice:
    if ( norm(ck, inf) > 100*eps )
        [b1, s1] = getCoeffs(ck, M, n);
    end
    
    if ( norm(dk, inf) > 100*eps )
        [b2, s2] = getCoeffs(dk, M, n);
    end
            
    % Construct the CF approximation now:
    a = [conj(a(n+1:-1:2)); a(1:n+1)] - b1 - flipud(b1) - 1i*b2 + 1i*flipud(b2);
    if ( norm(a(n+2:2*n+1) - conj(a(n:-1:1)), inf) > 100*eps || imag(a(n+1)) > 100*eps )
        error( 'why are the coeffs not that of a real function?')
    else
        a(n:-1:1) = conj(a(n+2:2*n+1));
        a(n+1) = real(a(n+1));
    end
    p = chebfun(a, 'coeffs', 'trig');
end

% Re-map p back to original domain:
p = newDomain(p, dom);

% re-normalize p:
p = p*normf;
% How to estimate s here? 1-norm is an upperbound?
% s = norm([s1, s2], 1);
s = normf*2*max(abs([s2, s1]));
end

%%
function [b, s] = getCoeffs(c, N, n)

% Solve the eigenvalue problem:
[V, D] = eig(hankel(c));
d = diag(D);
[s, i] = max(abs(d));
u = V(:,i);
% Compute the coefficients b recursively:
u1 = u(1);
while ( abs(u(1)) < eps )
    d(i) = [];
    V(:, i) = [];
    [s, i] = max(abs(d));
    u = V(:, i);
    u1 = u(1);
end
uu = u(2:(N-n));
b = c.';
for k = n:-1:-n
    b = [-(b(1:(N-n-1))*uu)/u1, b]; %#ok<AGROW>
end
b = b(1:2*n+1).';
end