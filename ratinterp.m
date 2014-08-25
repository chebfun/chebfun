function [p, q, r, mu, nu, poles, residues] = ratinterp(varargin)
%RATINTERP  Robust rational interpolation or least-squares approximation.
%   [P, Q, R_HANDLE] = RATINTERP(F, M, N) computes the (M, N) rational
%   interpolant of F on the M + N + 1 Chebyshev points of the second kind. F
%   can be a CHEBFUN, a function handle or a column vector of M + N + 1 data
%   points.  If F is a CHEBFUN, the rational interpolant is constructed on the
%   domain of F. Otherwise, the domain [-1, 1] is used. P and Q are CHEBFUNs
%   such that P(x)./Q(x) = F(x). R_HANDLE is an anonymous function evaluating
%   the rational interpolant directly.
%
%   [P, Q, R_HANDLE] = RATINTERP(F, M, N, NN) computes a (M, N) rational linear
%   least-squares approximant of F over the NN Chebyshev points of the second
%   kind. If NN = M + N + 1 or NN = [], a rational interpolant is computed.
%
%   [P, Q, R_HANDLE] = RATINTERP(F, M, N, NN, XI) computes a (M, N) rational
%   interpolant or approximant of F over the NN nodes XI. XI can also be one of
%   the strings 'type1', 'type2', 'unitroots' or 'equidistant', in which case
%   NN of the respective nodes are created on the interval [-1, 1].
%
%   [P, Q, R_HANDLE, MU, NU] = RATINTERP(F, M, N, NN, XI, TOL) computes a
%   robustified (M, N) rational interpolant or approximant of F over the NN + 1
%   nodes XI, in which components contributing less than the relative tolerance
%   TOL to the solution are discarded. If no value of TOL is specified, a
%   tolerance of 1e-14 is assumed; set TOL to zero to disable robustness. MU
%   and NU are the resulting numerator and denominator degrees. Note that if
%   the degree is decreased, a rational least-squares approximation is computed
%   over the NN points. The coefficients are computed relative to the
%   orthogonal basis derived from the nodes XI.
%
%   [P, Q, R_HANDLE, MU, NU, POLES, RES] = RATINTERP(F, M, N, NN, XI, TOL)
%   returns the poles POLES of the rational interpolant on the real axis as
%   well as the residues RES at those points. If any of the nodes XI lie in the
%   complex plane, the complex poles are returned as well.
%
%   [P, Q, R_HANDLE] = RATINTERP(D, F, M, N) computes the (M, N) rational
%   interpolant of F on the M + N + 1 Chebyshev points of the second kind on
%   the domain D, given as a two-element row vector.
%
%   Examples:
%
%   Compute a type-(10, 10) robustified rational interpolant to 1/(x -
%   0.2) on [-1, 1] in second-kind Chebyshev nodes:
%
%     [p, q, r] = ratinterp([-1 1], @(x) 1./(x - 0.2), 10, 10, [], 'type2');
%
%   Same thing but with robustness disabled:
%
%     [p, q, r] = ratinterp([-1 1], @(x) 1./(x - 0.2), 10, 10, [], 'type2', 0);
%
%   References:
%
%   [1] P. Gonnet, R. Pachon, and L. N. Trefethen, "ROBUST RATIONAL
%       INTERPOLATION AND LEAST-SQUARES", Electronic Transactions on
%       Numerical Analysis (ETNA), 38:146-167, 2011.
%
%   [2] R. Pachon, P. Gonnet and J. van Deun, "FAST AND STABLE RATIONAL
%       INTERPOLATION IN ROOTS OF UNITY AND CHEBYSHEV POINTS", Submitted to
%       SIAM Journal on Numerical Analysis, 2011.
%
% See also INTERP1, CHEBPADE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO:  Deal with array-valued CHEBFUNs / quasimatrices.

% Parse the inputs.
[dom, f, m, n, NN, xi, xi_type, tol] = parseInputs(varargin{:});

% Set up some values which we will use often.
N = NN - 1;
N1 = N + 1;
ts = tol*norm(f, inf);

% Check for symmetries.
[fEven, fOdd] = checkSymmetries(f, xi, xi_type, N, N1, ts);

% Form matrices for the linear system for the coefficients.
[Z, R] = assembleMatrices(f, n, xi, xi_type, N1);

% Compute coefficients of the numerator and denominator.
[b, n] = computeDenominatorCoeffs(Z, m, n, fEven, fOdd, N1, ts);

a = computeNumeratorCoeffs(f, m, n, xi_type, Z, b, fEven, fOdd, N, N1);
[a, b] = trimCoeffs(a, b, tol, ts);

% Get the exact numerator and denominator degrees.
mu = length(a) - 1;
nu = length(b) - 1;

% Build the numerator and denominator polynomials and create the output
% function handle for evaluating the rational approximation.
[p, q, r] = constructRatApprox(xi_type, R, a, b, mu, nu, dom);

% Compute poles and residues if requested.
if ( nargout > 5 )
    if ( nargout > 6 ) % Compute residues.
        % Compute partial fraction expansion of r.
        [residues, poles] = residue(p, q);
        [poles, ind] = sort(poles);
        residues = residues(ind);

        % Residues are the coefficients of 1/(x - poles(j))
        for j = 1:(length(poles) - 1)
            if (poles(j+1) == poles(j))
                residues(j+1) = residues(j);
            end
        end
    else               % Just compute the poles.
          poles = roots(q, 'all');
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [dom, f, m, n, NN, xi, xi_type, tol] = parseInputs(varargin)

% Make sure we have the correct number of arguments.
if ( nargin < 3 )
    error('CHEBFUN:ratinterp:tooFewArgs', 'Not enough input arguments.');
elseif ( nargin > 7 )
    error('CHEBFUN:ratinterp:tooManyArgs', 'Too many input arguments.');
end

% Deal with the domain input.
if ( nargin == 7 )
    dom = varargin{1};
    varargin(1) = [];

    if ( ~isempty(dom) && ...
        (~isfloat(dom) || ~isequal(size(dom), [1 2])) )
        error('CHEBFUN:ratinterp:badDom1', ...
            'Domain should be a 1 x 2 row vector of endpoints.');
    end

    if ( diff(dom) <= 0 )
        error('CHEBFUN:ratinterp:badDom2', 'Invalid domain.');
    end
elseif ( isempty(varargin{1}) || ...
        (isfloat(varargin{1}) && isequal(size(varargin{1}), [1 2])) )
    dom = varargin{1};
    varargin(1) = [];
else
    dom = [];
end

% Extract the rest of the input arguments. After this, all arguments will
% either be empty or have their user-supplied values.
varargin = [varargin repmat({[]}, [1 (6 - length(varargin))])];
[f, m, n, NN, xi, tol] = deal(varargin{:});

% Mandatory arguments are mandatory.
if ( isempty(f) || isempty(m) || isempty(n) )
    error('CHEBFUN:ratinterp:mandatoryArgs', ...
        'Interpolation data and polynomial degrees must be specified.');
end

% Set the domain if the user did not specify it.
if ( isempty(dom) )
    if ( isa(f, 'chebfun') )
        dom = domain(f);
        dom = dom([1 end]);
    else
        dom = [-1 1];
    end
end

% Ensure dom is not a domain object (since we lazily call diff(dom) below).
if ( isa(dom, 'domain') )
    warning('CHEBFUN:ratinterp:domainDeprecated', ...
        ['Using a DOMAIN object as an input to RATINTERP is deprecated.\n' ...
         'Specify domains using a two-element row vector instead.']);
    warning('off', 'CHEBFUN:ratinterp:domainDeprecated');
    dom = double(dom);
end

% Determine the number of interpolation nodes.
if ( isempty(NN) )
    if ( ~isempty(xi) && isfloat(xi) )
        NN = length(xi);
    elseif ( isfloat(f) )
        NN = length(f);
    else
        NN = m + n + 1;
    end
else
    if ( NN < m + n + 1 )
        error('CHEBFUN:ratinterp:nodeCount', ...
            'Number of nodes NN should be at least M + N + 1.');
    end
end

% Construct the interpolation nodes and get their type.
if ( isempty(xi) )
    xi = chebpts(NN, 2);
    xi_type = 'TYPE2';
elseif ( isfloat(xi) )                                      % Arbitrary nodes.
    if ( length(xi) ~= NN )
        error('CHEBFUN:ratinterp:lengthXI', ...
            'Input vector XI does not have the correct length.');
    end
    xi = 2.0 * ( xi - 0.5*sum(dom) ) / diff(dom);  % Scale nodes to [-1 1].
    xi_type = 'ARBITRARY';
elseif ( ischar(xi) )
    xi_type = xi;
    if ( strcmpi(xi, 'TYPE0') || strcmpi(xi, 'UNITROOTS') ) % Roots of unity.
        xi_type = 'TYPE0';
        xi = exp(2i*pi*(0:1:(NN - 1)).'/NN);
    elseif ( strcmpi(xi, 'TYPE1') )                         % 1st Chebyshev.
        xi = chebpts(NN, 1);
    elseif ( strcmpi(xi, 'TYPE2') )                         % 2nd Chebyshev.
        xi = chebpts(NN, 2);
    elseif ( strncmpi(xi, 'EQUI', 4) )                      % Equispaced.
        xi = linspace(-1, 1, NN).';
    else
        error('CHEBFUN:ratinterp:badXIType', ...
            'Unrecognized type of nodes XI.');
    end
else
    error('CHEBFUN:ratinterp:badXI', ...
        'XI must be a vector of nodes or a string.');
end

% If we were given a function handle or CHEBFUN, sample it on the grid.
if ( ~isfloat(f) )
    f = f(0.5*sum(dom) + 0.5*diff(dom)*xi);
elseif ( length(f) ~= NN )
    error( 'CHEBFUN:ratinterp:lengthF', ...
        'Input vector F does not have the correct length.');
end

% Set the default robustness tolerance.
if ( isempty(tol) )
    tol = 1e-14;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for computing the numerator and denominator coefficients.

function [fEven, fOdd] = checkSymmetries(f, xi, xi_type, N, N1, ts)

fEven = false;
fOdd = false;

if ( strncmpi(xi_type, 'TYPE', 4) )
    if ( xi_type(5) == '0' )  % Roots of unity.
        if ( mod(N, 2) == 1 )
            M = floor(N/2);

            fl = f(2:(M+1));
            fr = f((N+2-M):end);

            fEven = norm(fl - fr, inf) < ts;
            fOdd = norm(fl + fr, inf) < ts;
        end
    else                      % 1st/2nd-kind Chebyshev points.
        M = ceil(N/2);

        fl = f(1:M);
        fr = f(end:-1:(N1-M+1));

        fEven = norm(fl - fr, inf) < ts;
        fOdd = norm(fl + fr, inf) < ts;
    end
else                          % Other nodes.
    M = floor(N/2);
    [xs, ind] = sort(xi);
    xl = xs(1:(M+1));
    xr = xs(end:-1:(N1-M));
    if ( norm(xl + xr, inf) < ts )
        xi = xs;
        f = f(ind);
        M = ceil(N/2);

        fl = f(1:M);
        fr = f(end:-1:(N1-M+1));

        fEven = norm(fl - fr, inf) < ts;
        fOdd = norm(fl + fr, inf) < ts;
    end
end

end

function [Z, R] = assembleMatrices(f, n, xi, xi_type, N1)

if ( strncmpi(xi_type , 'type' , 4) )
    if ( xi_type(5) == '0' )         % Roots of unity.
        row = conj(fft(conj(f))) / N1;
        col = fft(f) / N1;
        col(1) = row(1);
        Z = toeplitz(col, row(1:(n+1)));
    elseif ( xi_type(5) == '1' )     % 1st-kind Chebyshev points.
        D = dct1(diag(f'))';
        Z = dct1(D(:,1:(n+1)));
    else                             % 2nd-kind Chebyshev points.
        D = idct2(eye(N1));
        Z = dct2(diag(f) * D(:,1:(n+1)));
    end
    R = [];
else                                 % Other nodes.
    C = ones(N1);
    C(:,2) = xi;
    for k = 3:N1
        C(:,k) = 2 * xi .* C(:,k-1) - C(:,k-2);
    end

    [C, R] = qr(C);
    Z = C' * diag(f) * C(:, 1:(n+1));
end

end

function [b, n] = computeDenominatorCoeffs(Z, m, n, fEven, fOdd, N1, ts)

shift = xor(fEven, mod(m, 2) == 1);

if ( (n > 0) && (~(fOdd || fEven) || (n > 1)) )
    while ( true )
        % Compute the SVD of the lower part of Z.
        if ~(fOdd || fEven)
            [U, S, V] = svd(Z((m+2):N1, 1:(n+1)), 0);
            ns = n;
            b = V(:,end);
        else
            [U, S, V] = svd(Z((m+2+shift):2:N1, 1:2:(n+1)), 0);
            ns = floor(n/2);
            b = zeros(n + 1, 1);
            b(1:2:end) = V(:,end);
        end

        % Get the smallest singular value.
        ssv = S(ns, ns);

        if ( ssv > ts )
            % Stop if converged.
            break;
        else
            % Reduce denominator degree.
            s = diag(S(1:ns,1:ns));
            if ( fEven || fOdd )
                n = n - 2*sum(s - ssv <= ts);
            else
                n = n - sum(s - ssv <= ts);
            end

            % Terminate if denominator is trivial.
            if ( n == 0 )
                b = 1;
                break;
            elseif ( n == 1 )
                if ( fEven )
                    b = [1 ; 0];
                    break;
                elseif ( fOdd )
                    b = [0 ; 1];
                    break;
                end
            end
        end
    end
elseif ( n > 0 )
    if ( fEven )
        b = [1 ; 0];
    elseif ( fOdd )
        b = [0 ; 1];
    end
else
    b = 1;
end

end

function a = computeNumeratorCoeffs(f, m, n, xi_type, Z, b, fEven, fOdd, N, N1)

if ( strncmpi(xi_type, 'type', 4) )
    if ( xi_type(5) == '0' )      % Roots of unity.
        a = fft(ifft(b, N1, 1) .* f);
        a = a(1:(m+1));
    elseif ( xi_type(5) == '1' )  % 1st-kind Chebyshev points.
        a = dct1(idct1([b ; zeros(N - n, 1)]) .* f);
        a = a(1:(m+1));
    elseif ( xi_type(5) == '2' )  % 2nd-kind Chebyshev points.
        a = dct2(idct2([b ; zeros(N - n, 1)]) .* f);
        a = a(1:(m+1));
    end
else
    a = Z(1:(m+1), 1:(n+1))*b;
end

if ( fEven )
    a(2:2:end) = 0;
elseif ( fOdd )
    a(1:2:end) = 0;
end

end

function [at, bt] = trimCoeffs(a, b, tol, ts)

at = a;
bt = b;

if ( tol > 0 )
    % Nonnegligible coefficients.
    nna = abs(at) > ts;
    nnb = abs(bt) > tol;

    % Discard trailing zeros.
    at = at(1:find(nna, 1, 'last'));
    bt = bt(1:find(nnb, 1, 'last'));

    % Remove small leading coefficients.
    while ( ~isempty(at) && ~isempty(bt) && (abs(at(1)) < ts) && ...
            (abs(bt(1)) < ts) )
        at = at(2:end);
        bt = bt(2:end);
    end
end

% Zero function special case.
if ( isempty(at) )
    at = 0;
    bt = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for assembling the rational interpolant.

function [p, q, r] = constructRatApprox(xi_type, R, a, b, mu, nu, dom)

if ( strncmpi(xi_type, 'type', 4) )
    if ( xi_type(5) == '0' )      % Roots of unity.
        [p, q, r] = constructRatApproxROU(a, b, mu, nu, dom);
    elseif ( xi_type(5) == '1' )  % 1st-kind Chebyshev points.
        [p, q, r] = constructRatApproxCheb1(a, b, mu, nu, dom);
    else                          % 2nd-kind Chebyshev points.
        [p, q, r] = constructRatApproxCheb2(a, b, mu, nu, dom);
    end
else                              % Arbitrary points.
        [p, q, r] = constructRatApproxArb(R, a, b, mu, nu, dom);
end

end

function [p, q, r] = constructRatApproxROU(a, b, mu, nu, dom)

md = 0.5 * sum(dom);
ihd = 2.0 / diff(dom);

% For speed, compute p(x) and q(x) using Horner's scheme and form CHEBFUNs.

% Build the numerator polynomial.
px = a(end) * ones(mu + 1, 1);
x = chebpts(mu + 1);
for k = mu:-1:1
    px = a(k) + x .* px;
end
p = chebfun(px, dom, 'tech', @chebtech2);

% Build the denominator polynomial and form the function handle.
if ( nu > 0 )
    qx = b(end) * ones(nu + 1, 1);
    x = chebpts(nu + 1);
    for k = nu:-1:1
        qx = b(k) + x .* qx;
    end
    q = chebfun(qx, dom, 'tech', @chebtech2);

    r = @(x) polyval(a((mu+1):-1:1) , ihd*(x - md)) ...
        ./ polyval(b((nu+1):-1:1), ihd*(x - md));
else
    q = chebfun(b, dom, 'tech', @chebtech2);
    r = @(x) polyval(a(mu+1:-1:1) , ihd*(x - md)) / b;
end

end

function [p, q, r] = constructRatApproxCheb1(a, b, mu, nu, dom)

md = 0.5 * sum(dom);
ihd = 2.0 / diff(dom);

% Build the numerator polynomial.
px = idct1(a);
p = chebfun(px, dom, 'tech', @chebtech1);

% Build the denominator polynomial and form the function handle.
if ( nu > 0 )
    qx = idct1(b);

    wp = sin((2*(0:mu) + 1)*pi/(2*(mu + 1)));
    wp(2:2:end) = -wp(2:2:end);
    wp = wp * 2^(mu - nu)/(mu + 1)*(nu + 1);

    wq = sin((2*(0:nu) + 1)*pi/(2*(nu + 1)));
    wq(2:2:end) = -wq(2:2:end);

    q = chebfun(qx, dom, 'tech', @chebtech1);

    r = @(x) ratbary(ihd*(x - md), px, qx, ...
        chebpts(mu + 1, 1), chebpts(nu + 1, 1), wp, wq);
else
    q = chebfun(b, dom, 'tech', @chebtech1);
    r = @(x) p(x)/b;
end

end

function [p, q, r] = constructRatApproxCheb2(a, b, mu, nu, dom)

md = 0.5 * sum(dom);
ihd = 2.0 / diff(dom);

% Build the numerator polynomial.
p = chebfun(a, dom, 'coeffs');

% Build the denominator polynomial and form the function handle.
if ( nu > 0 )
    q = chebfun(b, dom, 'coeffs');

    px = idct2(a);
    qx = idct2(b);

    wp = ones(1,(mu+1));
    wp(2:2:end) = -1;
    wp(1) = 0.5;
    wp(end) = 0.5*wp(end);
    wp = wp * (-2)^(mu - nu) / mu * nu;

    wq = ones(1,(nu+1));
    wq(2:2:end) = -1;
    wq(1) = 0.5;
    wq(end) = 0.5*wq(end);

    r = @(x) ratbary(ihd*(x - md), px, qx, ...
        chebpts(mu + 1, 2), chebpts(nu + 1, 2), wp, wq);
else
    q = chebfun(b, dom, 'tech', @chebtech2);
    r = @(x) p(x)/b;
end

end

function [p, q, r] = constructRatApproxArb(R, a, b, mu, nu, dom)

% Compute basis C at the mu + 1 and nu + 1 Chebyshev points and form a CHEBFUN.

% Build the numerator polynomial.
Cf = ones(mu + 1);
Cf(:,2) = chebpts(mu + 1);
for k = 3:(mu + 1)
    Cf(:,k) = 2 * Cf(:,2) .* Cf(:,k-1) - Cf(:,k-2);
end
Cf = Cf / R(1:(mu+1), 1:(mu+1));
p = chebfun(Cf(:,1:(mu+1)) * a, dom, 'tech', @chebtech2);

% Build the denominator polynomial and form the function handle.
if ( nu > 0 )
    Cf = ones(nu + 1);
    Cf(:,2) = chebpts(nu + 1);
    for k = 3:(nu + 1)
        Cf(:,k) = 2 * Cf(:,2) .* Cf(:,k-1) - Cf(:,k-2);
    end
    Cf = Cf / R(1:(nu+1), 1:(nu+1));
    q = chebfun(Cf(:,1:nu+1) * b, dom, 'tech', @chebtech2);
    r = @(x) p(x) ./ q(x);
else
    q = chebfun(b, dom, 'tech', @chebtech2);
    r = @(x) p(x)/b;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rational barycentric formula.

function y = ratbary(x, px, qx, xp, xq, wp, wq)
%RATBARY   Evaluate rational function using first-kind barycentric formula.
%   Y = RATBARY(X, PX, QX, XP, XQ, WP, WQ) evaluates the rational function P/Q
%   at the points in the array X, where P and Q are polynomials.  P and Q are
%   specified in Lagrange form by their values PX and QX at the nodes XP and
%   XQ, respectively.  WP and WQ are the barycentric weights for the nodes XP
%   and XQ, respectively.
%
%   The function works by using the first-kind barycentric formula to evaluate
%   P and Q individually and then takes the quotient.

y = zeros(size(x));

% Handle matrices of evaluation points one column at a time.
if ( (size(x, 1) > 1) && (size(x, 2) > 1) )
    for k = 1:size(x, 2)
        y(:,k) = ratbary(x(:,k), px, qx, xp, xq, wp, wq);
    end
    return
end

% Degrees of the polynomials involved (plus 1).
np = length(px);
nq = length(qx);

% Multiply by the barycentric weights.
pxw = px.' .* wp;
qxw = qx.' .* wq;

for i = 1:length(x)
    % Compute the sum in the first-kind barycentric formula for p.
    dxpinv = 1.0 ./ (x(i) - xp(:));
    ind = find(~isfinite(dxpinv));
    if ( ~isempty(ind) )
        y(i) = px(ind);
    else
        y(i) = (pxw * dxpinv);
    end

    % Compute the sum in the first-kind barycentric formula for q.
    dxqinv = 1.0 ./ (x(i) - xq(:));
    ind = find(~isfinite(dxqinv));
    if ( ~isempty(ind) )
        y(i) = y(i)/qx(ind);
    else
        y(i) = y(i) / (qxw * dxqinv);
    end
end

% Evaluate node polynomial for p.
llp = repmat(x(:), 1, np) - repmat(xp', length(x), 1);
lp = prod(llp, 2);
if ( ~isfinite(lp) )
    lp = exp(sum(log(llp), 2));
end
lp(lp == 0) = 1;

% Evaluate node polynomial for q.
llq = repmat(x(:), 1, nq) - repmat(xq', length(x), 1);
lq = prod(llq, 2 );
if ( ~isfinite(lq) )
    lq = exp(sum(log(llq), 2));
end
lq(lq == 0) = 1;

% Multiply by ratio of node polynomials and reshape if needed.
y = reshape(y(:) .* lp ./ lq , size(x));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete cosine transforms.

% DCT for Chebyshev points of the first kind.
function y = dct1(x)

n = size(x, 1);
w = (2/n)*(exp(-1i*(0:(n-1))*pi/(2*n))).';
w(1) = w(1)/sqrt(2);

if ( (mod(n, 2) == 1) || ~isreal(x) )
    y = fft([x ; x(n:-1:1,:)]);
    y = diag(w)*y(1:n,:);
else
    y = fft([x(1:2:n,:) ; x(n:-2:2,:)]);
    y = diag(2*w)*y;
end

if ( isreal(x) )
    y = real(y);
end

end

% iDCT for Chebyshev points of the first kind.
function x = idct1(y)

n = size(y, 1);
w = (n/2)*(exp(1i*(0:(n-1))*pi/(2*n))).';

if ( (mod(n, 2) == 1) || ~isreal(y) )
    w(1) = w(1)*sqrt(2);
    x = ifft([diag(w)*y
              zeros(1, size(y, 2))
              -1i*diag(w(2:n))*y(n:-1:2,:)]);
    x = x(1:n,:);
else
    w(1) = w(1)/sqrt(2);
    x([1:2:n n:-2:2],:) = ifft(diag(w)*y);
end

if ( isreal(y) )
    x = real(x);
end

end

% DCT for Chebyshev points of the second kind.
function c = dct2(v)

n = size(v, 1);
c = [v(end:-1:2,:) ; v(1:end-1,:)];

if ( isreal(v) )
    c = fft(c)/(2*n - 2);
    c = real(c);
elseif ( isreal(1i*v) )
    c = fft(imag(c))/(2*n - 2);
    c = 1i*real(c);
else
    c = fft(c)/(2*n - 2);
end

c = c(n:-1:1,:);
if ( n > 2 )
    c(2:end-1,:) = 2*c(2:end-1,:);
end
c = c(end:-1:1,:);

end

% iDCT for Chebyshev points of the second kind.
function v = idct2(c)

n = size(c, 1);
ii = 2:(n - 1);
c = c(end:-1:1,:);
c(ii,:) = c(ii,:)/2;
v = [c(end:-1:1,:) ; c(ii,:)];

if ( isreal(c) )
    v = real(ifft(v));
elseif ( isreal(1i*c) )
    v = 1i*real(ifft(imag(v)));
else
    v = ifft(v);
end

v = (n - 1)*[2*v(1,:) ; (v(ii,:) + v(2*n-ii,:)) ; 2*v(n,:)];
v = v(end:-1:1,:);

end
