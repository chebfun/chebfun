function [p, q, r, s] = cf(f, m, n, M)
%CF   Caratheodory-Fejer approximation
%   [P, Q, R_HANDLE] = CF(F, M, N) computes a type (M, N) rational CF
%   approximant to CHEBFUN F defined on [a, b], which must consist of just a
%   single FUN. P and Q are CHEBFUNs representing the numerator and denominator
%   polynomials. R_HANDLE is an anonymous function that evaluates P/Q.
%
%   [P, Q, R_HANDLE, S] = CF(F, M, N) also returns S, the associated CF singular
%   value, an approximation to the minimax error.
%
%   [P, Q, R_HANDLE, S] = CF(F, M, N, K) does the same but uses only the K-th
%   partial sum in Chebyshev expansion of F.
%
%   For polynomial CF approximation, use N = 0 or N = [] or only provide two
%   input arguments. If P and S are required, four output arguments must be
%   specified.
%
%   If F is a quasimatrix then so are the outputs P and Q, R_HANDLE is a cell
%   array of function handles and s is a vector.
%
%   Rational CF approximation can be very ill-conditioned for non-smooth
%   functions. If the program detects this, a warning message is given and the
%   results may be wrong.
%
%   CF = Caratheodory-Fejer approximation is a near-best approximation that is
%   often indistinguishable from the true best approximation (which in Chebfun
%   can be computed with MINIMAX), but often much faster to compute.
%
%   Examples:
%
%   Compute a quadratic polynomial CF approximant to exp(x) on [-1, 1]:
%
%     [p, q, r] = cf(chebfun(@exp), 2);
%
%   Compute a type-(4, 4) rational CF approximant to the same function:
%
%     [p, q, r] = cf(chebfun(@exp), 4, 4);
%
%   References:
%
%   [1] M. H. Gutknecht and L. N. Trefethen, "Real polynomial Chebyshev
%       approximation by the Caratheodory-Fejer method", SIAM J. Numer. Anal. 19
%       (1982), 358-371.
%
%   [2] L. N. Trefethen and M. H. Gutknecht, "The Caratheodory-Fejer method fpr
%       real rational approximation", SIAM J. Numer. Anal. 20 (1983), 420-436.
%
% See also AAA, CHEBPADE, MINIMAX, PADEAPPROX, RATINTERP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Do polynomial approximation if no denominator degree supplied.
if ( nargin < 3 )
    n = 0;
end

% Use full series expansion of F by default.
if ( nargin < 4 )
    M = length(f) - 1;
end

numCols = numColumns(f);
if ( numCols > 1 )
    if ( numel(f) == 1 )
        % Convert array-valued CHEBFUNs to quasimatrices first.
        f = cheb2quasi(f);
        [p, q, r, s] = cf(f, m, n, M);
    else
        % Deal with quasimatrices.
        r = cell(1, numCols);
        s = zeros(1, numCols);
        for k = 1:numCols
            [pk, qk, rk, sk] = cfOneColumn(f(:, k), m, n, M);
            p(k) = pk;
            q(k) = qk;
            r{k} = rk;
            s(k) = sk;
        end

        if ( f(1).isTransposed )
            p = p.';
            q = q.';
            r = r.';
            s = s.';
        end
    end
else
    [p, q, r, s] = cfOneColumn(f, m, n, M);
end

end

function [p, q, r, s] = cfOneColumn(f, m, n, M)

% Check the inputs.
if ( any(isinf(domain(f))) )
    error('CHEBFUN:CHEBFUN:cf:unboundedDomain', ...
        'CF does not work for CHEBFUNs with unbounded domains.');
end

if ( issing(f) )
    error('CHEBFUN:cf:singularFunction', ...
        'CF does not support functions with singularities.');
end

if ( (numel(f.funs) > 1) && (nargin < 4) )
    error('CHEBFUN:CHEBFUN:cf:multipleFuns', ...
        'For CHEBFUNs with multiple FUNs, CF must be called with 4 arguments.');
end

% Form global polynomial approximation if f is only piecewise smooth.
if ( numel(f.funs) > 1 )
    f = chebfun(@(x) feval(f, x), f.domain([1, end]), M + 1);
end

% Trivial case: approximation length exceeds that of the expansion length.
if ( m >= M )
    p = f;
    q = chebfun(1, domain(f));
    r = @(x) feval(p, x);
    s = 0;
    return
end

% Extract the Chebyshev coefficients to be used in computing the approximation.
a = chebcoeffs(f, length(f));
a = a(1:M+1);

% Deal with complex-valued functions.
if ( any(imag(a) ~= 0) )
    warning('CHEBFUN:CHEBFUN:cf:complex', ...
        'CF does not work for complex valued functions. Taking real part.');
    a = real(a);
end

% Compute the CF approximation.
if ( isempty(n) || (n == 0) )
    [p, q, r, s] = polynomialCF(f, a, m, M);
else
    [p, q, r, s] = rationalCF(f, a(end:-1:1).', m, n, M);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial CF approximation.

function [p, q, r, s] = polynomialCF(f, a, m, M)

dom = domain(f);

% Trivial case:  approximation length is the length of the CHEBFUN.
if ( m == M - 1 )
    p = chebfun(a(1:M), dom, 'coeffs');
    q = chebfun(1, dom);
    r = @(x) feval(p, x);
    s = abs(a(M+1));
    return
end

c = a(m+2:M+1);
if ( length(c) > 1024 )
    opts.disp = 0;
    opts.issym = 1;
    opts.isreal = 1;
    opts.v0 = ones(length(c), 1)/length(c);
    [V, D] = eigs(hankel(c), 1, 'lm', opts);
else
    [V, D] = eig(hankel(c));
end

[s, i] = max(abs(diag(D)));

u = V(:,i);
u1 = u(1);
uu = u(2:(M-m));

b = c;
for k = m:-1:-m
    b = [-(b(1:(M-m-1)).'*uu)/u1; b]; %#ok<AGROW>
end
bb = b(m+1:2*m+1) + [0; b(m:-1:1)];
pk = a(1:m+1)-bb;
p = chebfun(pk, dom, 'coeffs');
q = chebfun(1, dom);
r = @(x) feval(p,x);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rational CF approximation.

function [p, q, r, s] = rationalCF(f, a, m, n, M)

dom = domain(f);

tolfft = 1e-14;   % Relative tolerance.
maxnfft = 2^17;   % Maximum FFT length.

% Tolerances for detecting ill-conditioning.
tolCond = 1e13;
tolCondG = 1e3;

% Scaling of T_0 coefficient.
a(end) = 2*a(end);

% Check even / odd symmetries.
if ( max(abs(a(end-1:-2:1)))/vscale(f) < eps )    % f is even.
    if ( ~(mod(m, 2) || mod(n, 2)) )
        m = m + 1;
    elseif ( mod(m, 2) && mod(n, 2) )
        n = n - 1;
        if ( n == 0 )
            [p, q, r, s] = cf(f, m, n, M);
            return;
        end
    end
elseif ( max(abs(a(end:-2:1)))/vscale(f) < eps )  % f is odd.
    if ( mod(m,2) && ~mod(n,2) )
        m = m + 1;
    elseif ( ~mod(m,2) && mod(n,2) )
        n = n - 1;
        if (n == 0)
            [p, q, r, s] = cf(f, m, n, M);
            return;
        end
    end
end

% Obtain eigenvalues and block structure.
[s, u, k, l, rflag] = getBlock(a, m, n, M);
if ( (k > 0) || (l > 0) )
    % f is rational (at least up to machine precision).
    if ( rflag )
        [p, q, r] = chebpade(f, m - k, n - k);
        s = eps;
        %warning('CHEBFUN:CHEBFUN:cf:chebpade', ...
        %  'Function looks close to rational; switching to CHEBPADE.');
        return;
    end

    nnew = n - k;
    [s, u, knew, lnew] = getBlock(a, m + l, nnew, M);
    if ( (knew > 0) || (lnew > 0) )
        n = n + l;
        [s, u, k, l] = getBlock(a, m - k, n, M);
    else
        n = nnew;
    end
end

% Obtain polynomial q from Laurent coefficients using FFT.
N = max(2^nextpow2(length(u)), 256);
ud = polyder(u(end:-1:1));
ud = ud(end:-1:1).';

ac = fft(conj(fft(ud, N)./fft(u, N)))/N;
act = zeros(N, 1);
while ( (norm(1 - act(end-n:end-1)./ac(end-n:end-1), inf) > tolfft) && ...
        (N < maxnfft) )
    act = ac; N = 2*N;
    ac = fft(conj(fft(ud, N)./fft(u, N)))/N;
end
ac = real(ac);

b = ones(1, n + 1);
for j = 1:n,
    b(j+1) = -(b(1:j)*ac(end-j:end-1))/j;
end

z = roots(b);
if ( any(abs(z) > 1) )
    warning('CHEBFUN:CHEBFUN:cf:illConditioned', ...
      'Ill-conditioning detected. Results may be inaccurate.');
end

z = z(abs(z) < 1);
rho = 1/max(abs(z));
z = .5*(z + 1./z);

% Compute q from the roots for stability reasons.
qt = newDomain(chebfun(@(x) real(prod(x - z)/prod(-z)), 'vectorize'), dom);
q = chebfun;
q = defineInterval(q, dom, qt);

% Compute Chebyshev coefficients of approximation Rt from Laurent coefficients
% of Blaschke product using FFT.
v = u(end:-1:1);
N = max(2^nextpow2(length(u)), 256);

ac = fft(exp(2*pi*1i*M*(0:(N-1))'/N).*conj(fft(u, N)./fft(v, N)))/N;
act = zeros(N, 1);

while ( (norm(1 - act(1:(m+1))./ac(1:(m+1)), inf) > tolfft) && ...
        (norm(1 - act(end-m+1:end)./ac(end-m+1:end), inf) > tolfft) && ...
        (N < maxnfft) )
    act = ac; N = 2*N;
    ac = fft(exp(2*pi*1i*M*(0:N-1)'/N).*conj(fft(u, N)./fft(v, N)))/N;
end

ac = s*real(ac);
ct = a(end:-1:(end-m)) - ac(1:(m+1))' - [ac(1) ac(end:-1:(end-m+1))'];
s = abs(s);

% Compute numerator polynomial from Chebyshev expansion of 1./q and Rt.  We
% know the exact ellipse of analyticity for 1./q, so use this knowledge to
% obtain its Chebyshev coefficients (see line below).
qRecip = chebfun(@(x) 1./feval(q, x), dom, ceil(log(4/eps/(rho - 1))/log(rho)));
gam = flipud(chebcoeffs(qRecip, length(qRecip)));
gam = gam.';
gam = [zeros(1, 2*m + 1 - length(gam)) gam];
gam = gam(end:-1:end-2*m);
gam(1) = 2*gam(1);
gam = toeplitz(gam);

% The following steps reduce the Toeplitz system of size 2*m + 1 to a system of
% size m, and then solve it.  If q has zeros close to the domain, then G is
% ill-conditioned, and accuracy is lost.
A = gam(1:m,1:m);
B = gam(1:m,m+1);
C = gam(1:m,end:-1:m+2);
G = A + C - 2*(B*B')/gam(1,1);

if ( (cond(G)/s > tolCond) && (cond(G) > tolCondG) )
  warning('CHEBFUN:CHEBFUN:cf:illConditioned', ...
    'Ill-conditioning detected. Results may be inaccurate.');
end

bc = G\(-2*(B*ct(1)/gam(1,1) -ct(m+1:-1:2)'));
bc0 = (ct(1) - B'*bc)/gam(1,1);
bc = [bc0, bc(end:-1:1)'];
p = chebfun(bc.', dom, 'coeffs');
r = @(x) feval(p, x)./feval(q, x);

end

function [s, u, k, l, rFlag] = getBlock(a, m, n, M)
% Each Hankel matrix corresponds to one diagonal m - n = const in the CF-table;
% when a diagonal intersects a square block, the eigenvalues on the
% intersection are all equal.  k and l tell you how many entries on the
% intersection appear before and after the eigenvalues under consideration.  u
% is the corresponding eigenvector

tol = 1e-14;

if ( n > M + m + 1 )
    c = zeros(1, n - m - M - 1);
    nn = M + m + 1;
else
    c = [];
    nn = n;
end

c = [c, a(M + 1 - abs(m - nn + 1:M))];
if ( length(c) > 1024 )
    % Use eigs() if the matrix is large.
    opts.disp = 0;
    opts.issym = 1;
    opts.isreal = 1;
    opts.v0 = ones(length(c), 1)/length(c);
    [V, D] = eigs(hankel(c), min(n + 10, length(c)), 'lm', opts);
else
    % For modest-sized matrices, just use eig().
    [V, D] = eig(hankel(c));
end

[S, i] = sort(abs(diag(D)), 'descend');
s = D(i(n+1), i(n+1));
u = V(:,i(n+1));

tmp = abs(S - abs(s)) < tol;
k = 0;
l = 0;
while ( (k < n) && tmp(n-k) )
    k = k + 1;
end

while ( ((n + l + 2) < length(tmp)) && tmp(n + l + 2) )
    l = l + 1;
end

% Flag indicating if the function is actually rational.
rFlag = (n + l + 2) == length(tmp);

end
