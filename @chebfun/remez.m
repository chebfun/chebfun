function varargout = remez(f, varargin)
%REMEZ   Best polynomial or rational approximation for real valued chebfuns.
%   P = REMEZ(F, M) computes the best polynomial approximation of degree M to
%   the real CHEBFUN F in the infinity norm using the Remez algorithm.
%
%   [P, Q] = REMEZ(F, M, N) computes the best rational approximation P/Q of type
%   (M, N) to the real CHEBFUN F using the Remez algorithm.
%
%   [P, Q, R_HANDLE] = REMEZ(F, M, N) does the same but additionally returns a
%   function handle R_HANDLE for evaluating the rational function P/Q.
%
%   [...] = REMEZ(..., 'tol', TOL) uses the value TOL as the termination
%   tolerance on the increase of the levelled error.
%
%   [...] = REMEZ(..., 'display', 'iter') displays output at each iteration.
%
%   [...] = REMEZ(..., 'maxiter', MAXITER) sets the maximum number of allowable
%   iterations to MAXITER.
%
%   [...] = REMEZ(..., 'plotfcns', 'error') plots the error after each iteration
%   while the algorithm executes.
%
%   [P, ERR] = REMEZ(...) and [P, Q, R_HANDLE, ERR] = REMEZ(...) also returns
%   the maximum error ERR.
%
%   [P, ERR, STATUS] = REMEZ(...) and [P, Q, R_HANDLE, ERR, STATUS] = REMEZ(...)
%   also return a structure array STATUS with the following fields:
%      STATUS.DELTA  - Obtained tolerance.
%      STATUS.ITER   - Number of iterations performed.
%      STATUS.DIFFX  - Maximum correction in last trial reference.
%      STATUS.XK     - Last trial reference on which the error equioscillates.
%
%   This code is quite reliable for polynomial approximations but rather
%   fragile for rational approximations.  Better results can often be obtained
%   with CF(), especially if f is smooth.
%
% References:
%
%   [1] Pachon, R. and Trefethen, L. N.  "Barycentric-Remez algorithms for best
%   polynomial approximation in the chebfun system", BIT Numerical Mathematics,
%   49:721-742, 2009.
%
%   [2] Pachon, R.  "Algorithms for Polynomial and Rational Approximation".
%   D. Phil. Thesis, University of Oxford, 2010 (Chapter 6).
%
% See also CF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

dom = f.domain([1, end]);
normf = norm(f);

if ( ~isreal(f) )
    error('CHEBFUN:CHEBFUN:remez:real', ...
        'REMEZ only supports real valued functions.');
end

if ( numColumns(f) > 1 )
    error('CHEBFUN:CHEBFUN:remez:quasi', ...
        'REMEZ does not currently support quasimatrices.');
end

if ( issing(f) )
    error('CHEBFUN:CHEBFUN:remez:singularFunction', ...
        'REMEZ does not currently support functions with singularities.');
end

% Parse the inputs.
[m, n, N, rationalMode, opts] = parseInputs(f, varargin{:});

% With zero denominator degree, the denominator polynomial is trivial.
if ( n == 0 )
    qk = 1;
    q = chebfun(1, dom);
    qmin = q;
end

% Initial values for some parameters.
iter = 0;                 % Iteration count.
delta = max(normf, eps);  % Value for stopping criterion.
deltamin = inf;           % Minimum error encountered.
diffx = 1;                % Maximum correction to trial reference.

% Compute an initial reference set to start the algorithm.
xk = getInitialReference(f, m, n, N);
xo = xk;

% Print header for text output display if requested.
if ( opts.displayIter )
    disp('It.     Max(|Error|)       |ErrorRef|      Delta ErrorRef      Delta Ref')
end

% Run the main algorithm.
while ( (delta/normf > opts.tol) && (iter < opts.maxIter) && (diffx > 0) )
    fk = feval(f, xk);     % Evaluate on the exchange set.
    w = baryWeights(xk);   % Barycentric weights for exchange set.

    % Compute trial function and levelled reference error.
    if ( n == 0 )
        [p, h] = computeTrialFunctionPolynomial(fk, xk, w, m, N, dom);
    else
        [p, q, h] = computeTrialFunctionRational(fk, xk, w, m, n, N, dom);
    end

    % Perturb exactly-zero values of the levelled error.
    if ( h == 0 )
        h = 1e-19;
    end

    % Update the exchange set using the Remez algorithm with full exchange.
    [xk, err, err_handle] = exchange(xk, h, 2, f, p, q, N + 2);

    % If overshoot, recompute with one-point exchange.
    if ( err/normf > 1e5 )
        [xk, err, err_handle] = exchange(xo, h, 1, f, p, q, N + 2);
    end

    % Update max. correction to trial reference and stopping criterion.
    diffx = max(abs(xo - xk));
    delta = err - abs(h);

    % Store approximation with minimum norm.
    if ( delta < deltamin )
        pmin = p;
        if ( n > 0 )
            qmin = q;
        end

        errmin = err;
        xkmin = xk;
        deltamin = delta;
    end

    % Display diagnostic information as requested.
    if ( opts.plotIter )
        doPlotIter(xo, xk, err_handle, dom);
    end

    if ( opts.displayIter )
        doDisplayIter(iter, err, h, delta, normf, diffx);
    end

    xo = xk;
    iter = iter + 1;
end

% Take best results of all the iterations we ran.
p = pmin;
err = errmin;
xk = xkmin;
delta = deltamin;

% Warn the user if we failed to converge.
if ( delta/normf > opts.tol )
    warning('CHEBFUN:CHEBFUN:remez:convergence', ...
        ['Remez algorithm did not converge after ', num2str(iter), ...
         ' iterations to the tolerance ', num2str(opts.tol), '.']);
end

% Form the outputs.
status.delta = delta/normf;
status.iter = iter;
status.diffx = diffx;
status.xk = xk;

p = simplify(p);
if ( rationalMode )
    q = simplify(qmin);
    varargout = {p, q, @(x) feval(p, x)./feval(q, x), err, status};
else
    varargout = {p, err, status};
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [m, n, N, rationalMode, opts] = parseInputs(f, varargin)

% Detect polynomial / rational approximation type and parse degrees.
if ( ~mod(nargin, 2) ) % Even number of inputs --> polynomial case.
    m = varargin{1};
    n = 0;
    rationalMode = false;
    varargin = varargin(2:end);
else                   % Odd number of inputs --> rational case.
    [m, n] = adjustDegreesForSymmetries(f, varargin{1}, varargin{2});
    rationalMode = true;
    varargin = varargin(3:end);
end

N = m + n;

% Parse name-value option pairs.
opts.tol = 1e-16*(N^2 + 10); % Relative tolerance for deciding convergence.
opts.maxIter = 20;           % Maximum number of allowable iterations.
opts.displayIter = false;    % Print output after each iteration.
opts.plotIter = false;       % Plot approximation at each iteration.

for k = 1:2:length(varargin)
    if ( strcmpi('tol', varargin{k}) )
        opts.tol = varargin{k+1};
    elseif ( strcmpi('maxiter', varargin{k}) )
        opts.maxIter = varargin{k+1};
    elseif ( strcmpi('display', varargin{k}) )
        opts.displayIter = true;
    elseif ( strcmpi('plotfcns', varargin{k}) )
        opts.plotIter = true;
    else
        error('CHEBFUN:CHEBFUN:remez:badInput', ...
            'Unrecognized sequence of input parameters.')
    end
end

end

function [m, n] = adjustDegreesForSymmetries(f, m, n)
%ADJUSTDEGREESFORSYMMETRIES   Adjust rational approximation degrees to account
%   for function symmetries.
%
%   [M, N] = ADJUSTDEGREESFORSYMMETRIES(F, M, N) returns new degrees M and N to
%   correct the defect of the rational approximation if the target function is
%   even or odd.  In either case, the Walsh table is covered with blocks of
%   size 2x2, e.g.  for even function the best rational approximant is the same
%   for types [m/n], [m+1/n], [m/n+1] and [m+1/n+1], with m and n even. This
%   strategy is similar to the one proposed by van Deun and Trefethen for CF
%   approximation in Chebfun (see @chebfun/cf.m).

% Sample piecewise-smooth CHEBFUNs.
if ( (numel(f.funs) > 1) || (length(f) > 128) )
  f = chebfun(f, f.domain([1, end]), 128);
end

% Compute the Chebyshev coefficients.
c = chebcoeffs(f, length(f));
c(end) = 2*c(end);

% Check for symmetries and reduce degrees accordingly.
if ( max(abs(c(end-1:-2:1)))/vscale(f) < eps )   % f is even.
    if ( mod(m, 2) == 1 )
        m = m - 1;
    end
    if ( mod(n, 2) == 1 )
        n = n - 1;
    end
elseif ( max(abs(c(end:-2:1)))/vscale(f) < eps ) % f is odd.
    if ( mod(m, 2) == 0 )
        m = m - 1;
    end
    if ( mod(n, 2) == 1 )
        n = n - 1;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implementing the core part of the algorithm.

function xk = getInitialReference(f, m, n, N)

% If doing rational Remez, get initial reference from trial function generated
% by CF or Chebyshev-Pade.
flag = 0;
if ( n > 0 )
    if ( numel(f.funs) == 1 )
        %[p, q] = chebpade(f, m, n);
        [p, q] = cf(f, m, n);
    else
        %[p, q] = chebpade(f, m, n, 5*N);
        [p, q] = cf(f, m, n, 5*N);
    end
    [xk, err, e, flag] = exchange([], 0, 2, f, p, q, N + 2);
end

% In the polynomial case or if the above procedure failed to produce a reference
% with enough equioscillation points, just use the Chebyshev points.
if ( flag == 0 )
    xk = chebpts(N + 2, f.domain([1, end]));
end

xo = xk;

end

function [p, h] = computeTrialFunctionPolynomial(fk, xk, w, m, N, dom)

% Vector of alternating signs.
sigma = ones(N + 2, 1);
sigma(2:2:end) = -1;

h = (w'*fk) / (w'*sigma);                          % Levelled reference error.
pk = (fk - h*sigma);                               % Vals. of r*q in reference.

% Trial polynomial.
p = chebfun(@(x) bary(x, pk, xk, w), dom, m + 1);

end

function [p, q, h] = computeTrialFunctionRational(fk, xk, w, m, n, N, dom)

% Vector of alternating signs.
sigma = ones(N + 2, 1);
sigma(2:2:end) = -1;

% Orthogonal matrix with respect to <,>_{xk}.
[C, ignored] = qr(fliplr(vander(xk)));

% Left and right rational interpolation matrices.
ZL = C(:,m+2:N+2).'*diag(fk)*C(:,1:n+1);
ZR = C(:,m+2:N+2).'*diag(sigma)*C(:,1:n+1);

% Solve generalized eigenvalue problem.
[v, d] = eig(ZL, ZR);

% Compute all possible qk and and look for ones with unchanged sign.
qk_all = C(:,1:n+1)*v;
pos =  find(abs(sum(sign(qk_all))) == N + 2);  % Sign changes of each qk.

if ( isempty(pos) || (length(pos) > 1) )
    error('CHEBFUN:CHEBFUN:remez:badGuess', ...
        'Trial interpolant too far from optimal');
end

qk = qk_all(:,pos);       % Keep qk with unchanged sign.
h = d(pos, pos);          % Levelled reference error.
pk = (fk - h*sigma).*qk;  % Vals. of r*q in reference.

% Trial numerator and denominator.
p = chebfun(@(x) bary(x, pk, xk, w), dom, m + 1);
q = chebfun(@(x) bary(x, qk, xk, w), dom, n + 1);

end

function [xk, norme, err_handle, flag] = exchange(xk, h, method, f, p, q, Npts)
%EXCHANGE   Modify an equioscillation reference using the Remez algorithm.
%   EXCHANGE(XK, H, METHOD, F, P, Q) performs one step of the Remez algorithm
%   for the best rational approximation of the CHEBFUN F of the target function
%   according to the first method (METHOD = 1), i.e. exchanges only one point,
%   or the second method (METHOD = 2), i.e. exchanges all the reference points.
%   XK is a column vector with the reference, H is the levelled error, P is the
%   numerator, and Q is the denominator of the trial rational function P/Q.
%
%   [XK, NORME, E_HANDLE, FLAG] = EXCHANGE(...) returns the modified reference
%   XK, the supremum norm of the error NORME (included as an output argument,
%   since it is readily computed in EXCHANGE and is used later in REMEZ), a
%   function handle E_HANDLE for the error, and a FLAG indicating whether there
%   were at least N+2 alternating extrema of the error to form the next
%   reference (FLAG = 1) or not (FLAG = 0).
%
%   [XK, ...] = EXCHANGE([], 0, METHOD, F, P, Q, N + 2) returns a grid of N + 2
%   points XK where the error F - P/Q alternates in sign (but not necessarily
%   equioscillates). This feature of EXCHANGE is useful to start REMEZ from an
%   initial trial function rather than an initial trial reference.

% Compute extrema of the error.
e_num = (q.^2).*diff(f) - q.*diff(p) + p.*diff(q);
rts = roots(e_num, 'nobreaks');
rr = [f.domain(1) ; rts; f.domain(end)];

% Function handle output for evaluating the error.
err_handle = @(x) feval(f, x) - feval(p, x)./feval(q, x);

% Select exchange method.
if ( method == 1 )                           % One-point exchange.
    [ignored, pos] = max(abs(feval(err_handle, rr)));
    pos = pos(1);
else                                           % Full exchange.
    pos = find(abs(err_handle(rr)) >= abs(h)); % Values above levelled error
end

% Add extrema nearest to those which are candidates for exchange to the
% existing exchange set.
[r, m] = sort([rr(pos) ; xk]);
v = ones(Npts, 1);
v(2:2:end) = -1;
er = [feval(err_handle, rr(pos)) ; v*h];
er = er(m);

% Delete repeated points.
repeated = diff(r) == 0;
r(repeated) = [];
er(repeated) = [];

% Determine points and values to be kept for the reference set.
s = r(1);    % Points to be kept.
es = er(1);  % Values to be kept.
for i = 2:length(r)
    if ( (sign(er(i)) == sign(es(end))) && (abs(er(i)) > abs(es(end))) )
        % Given adjacent points with the same sign, keep one with largest value.
        s(end) = r(i);
        es(end) = er(i);
    elseif ( sign(er(i)) ~= sign(es(end)) )
        % Keep points which alternate in sign.
        s = [s ; r(i)];    %#ok<AGROW>
        es = [es ; er(i)]; %#ok<AGROW>
    end
end

% Of the points we kept, choose n + 2 consecutive ones that include the maximum
% of the error.
[norme, index] = max(abs(es));
d = max(index - Npts + 1, 1);
if ( Npts <= length(s) )
    xk = s(d:d+Npts-1);
    flag = 1;
else
    xk = s;
    flag = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for displaying diagnostic information.

% Function called when opts.plotIter is set.
function doPlotIter(xo, xk, err_handle, dom)

xxk = linspace(dom(1), dom(end), 300);
plot(xo, 0*xo, 'or', 'MarkerSize', 12)   % Old reference.
holdState = ishold;
hold on
plot(xk, 0*xk, '*k', 'MarkerSize', 12)   % New reference.
plot(xxk, err_handle(xxk))               % Error function.
if ( ~holdState )                        % Return to previous hold state.
    hold off
end
xlim(dom)
legend('Current Ref.', 'Next Ref.', 'Error')
drawnow

end

% Function called when opts.displayIter is set.
function doDisplayIter(iter, err, h, delta, normf, diffx)

disp([num2str(iter), '        ', num2str(err, '%5.4e'), '        ', ...
    num2str(abs(h), '%5.4e'), '        ', ...
    num2str(delta/normf, '%5.4e'), '        ', num2str(diffx, '%5.4e')])

end
