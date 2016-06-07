function varargout = trigremez(f, varargin)
%TRIGREMEZ   Best trigonometric polynomial approximation for real-valued chebfuns.
%   P = TRIGREMEZ(F, M) computes the best trigonometric polynomial approximation 
%   of degree M to the real CHEBFUN F in the infinity norm using the Remez 
%   algorithm. F may be Chebyshev interpolant to an underlying continuous function 
%   or a trigonometric interpolant to an underlying continous function on the 
%   unit circle.
%
%   [...] = TRIGREMEZ(..., 'tol', TOL) uses the value TOL as the termination
%   tolerance on the increase of the levelled error.
%
%   [...] = TRIGREMEZ(..., 'display', 'iter') displays output at each iteration.
%
%   [...] = TRIGREMEZ(..., 'maxiter', MAXITER) sets the maximum number of allowable
%   iterations to MAXITER.
%
%   [...] = TRIGREMEZ(..., 'plotfcns', 'error') plots the error after each iteration
%   while the algorithm executes.
%
%   [P, ERR] = TRIGREMEZ(...) returns the maximum error ERR.
%
%   [P, ERR, STATUS] = TRIGREMEZ(...)  also return a structure array STATUS 
%   with the following fields:
%      STATUS.DELTA  - Obtained tolerance.
%      STATUS.ITER   - Number of iterations performed.
%      STATUS.DIFFX  - Maximum correction in last trial reference.
%      STATUS.XK     - Last trial reference on which the error equioscillates.
%
%
% References:
%
%   [1] Javed, M. and Trefethen, L. N.  "Remez and CF approximations of 
%    Periodic Functions". In preparation.
%
% See also REMEZ, CF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    varargout = {f};
    return;
end

normf = norm(f);

if ( ~isreal(f) )
    error('CHEBFUN:CHEBFUN:trigremez:real', ...
        'TRIGREMEZ only supports real-valued functions.');
end

if ( numColumns(f) > 1 )
    error('CHEBFUN:CHEBFUN:trigremez:quasi', ...
        'TRIGREMEZ does not currently support quasimatrices.');
end

if ( issing(f) )
    error('CHEBFUN:CHEBFUN:trigremez:singularFunction', ...
        'TRIGREMEZ does not currently support functions with singularities.');
end

if ( any(isinf(abs(domain(f)))) )
    error('CHEBFUN:CHEBFUN:trigremez:unboundedDomain', ...
        'Function must be defined on a bounded domain.');
end

if ( isdelta(f) )
    error('CHEBFUN:CHEBFUN:trigremez:deltaFunctions', ...
        'Function must not have any delta functions.');
end

% Parse the inputs.
[m, N, opts] = parseInputs(f, varargin{:});

% Initial values for some parameters.
iter = 0;       % Iteration count.
delta = normf;  % Value for stopping criterion.
deltamin = inf; % Minimum error encountered.
diffx = 1;      % Maximum correction to trial reference.

% Map everything to [-pi, pi]:
dom = f.domain([1, end]);
a = dom(1);
b = dom(end);
f = newDomain(f, [-pi, pi]);

% Compute an initial reference set to start the algorithm:
xk = trigpts(N, [-pi, pi]);
xo = xk;

% Print header for text output display if requested.
if ( opts.displayIter )
    disp('It.     Max(|Error|)       |ErrorRef|      Delta ErrorRef      Delta Ref')
end

% Run the main algorithm.
while ( (delta/normf > opts.tol) && (iter < opts.maxIter) && (diffx > 0) )
    fk = feval(f, xk);     % Evaluate on the exchange set.
    w = trigBaryWeights(xk);
    
    % Compute trial function and levelled reference error.
    [p, h] = computeTrialFunctionPolynomial(fk, xk, w, m, N, [-pi, pi]);
    
    % Perturb exactly-zero values of the levelled error.
    if ( h == 0 )
        h = 1e-19;
    end

    % Update the exchange set using the Remez algorithm with full exchange.
    [xk, err, err_handle] = exchange(xk, h, 2, f, p, N);

    % If overshoot, recompute with one-point exchange.
    if ( err/normf > 1e5 )
        [xk, err, err_handle] = exchange(xo, h, 1, f, p, N);
    end

    % Update max. correction to trial reference and stopping criterion.
    diffx = max(abs(xo - xk));
    delta = err - abs(h);

    % Store approximation with minimum norm.
    if ( delta < deltamin )
        pmin = p;
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

% Map the points back on the original domain:
forwardMap = @(y) b*(y + pi)/(2*pi) + a*(pi - y)/(2*pi); 
xk = forwardMap(xkmin);

delta = deltamin;

% Warn the user if we failed to converge.
if ( delta/normf > opts.tol )
    warning('CHEBFUN:CHEBFUN:trigremez:convergence', ...
        ['Remez algorithm did not converge after ', num2str(iter), ...
         ' iterations to the tolerance ', num2str(opts.tol), '.']);
end

% Form the outputs.
status.delta = delta/normf;
status.iter = iter;
status.diffx = diffx;
status.xk = xk;

% Map the approximation back to the original domain:
p = newDomain(p, [a, b]);

% return:
varargout = {p, err, status};

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [m, N, opts] = parseInputs(f, varargin)

m = varargin{1};

if ( m < 0 || m ~= round(m) )
    error('CHEBFUN:CHEBFUN:trigremez:parseInputs', ...
        'Degree of approximation must be a nonnegative integer.');
end

varargin = varargin(2:end);

% Number of points for equioscillation:
N = 2*m+2;

% Parse name-value option pairs.
baseTol = 1e-12;
opts.tol = baseTol*(N^2 + 10); % Relative tolerance for deciding convergence.
opts.maxIter = 100;            % Maximum number of allowable iterations.
opts.displayIter = false;      % Print output after each iteration.
opts.plotIter = false;         % Plot approximation at each iteration.

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
        error('CHEBFUN:CHEBFUN:trigremez:badInput', ...
            'Unrecognized sequence of input parameters.')
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implementing the core part of the algorithm.

function [p, h] = computeTrialFunctionPolynomial(fk, xk, w, m, N, dom)

% Vector of alternating signs.
sigma = ones(N, 1);
sigma(2:2:end) = -1;

h = (w'*fk) / (w'*sigma);             % Levelled reference error.
pk = (fk - h*sigma);                  % Vals. of r*q in reference.

% Trial polynomial by interpolation:
p = chebfun(@(x) trigBary(x, pk, xk, dom), dom, 2*m+1, 'trig');

end


% TODO: There is a lot of overlap between the present case of trigonometric
% remez and the usual polynomial remez. The common parts of the code should
% factored out in future.

function [xk, norme, err_handle, flag] = exchange(xk, h, method, f, p, Npts)
%EXCHANGE   Modify an equioscillation reference using the Remez algorithm.
%   EXCHANGE(XK, H, METHOD, F, P) performs one step of the Remez algorithm
%   for the best polynomial approximation of the CHEBFUN F of the target function
%   according to the first method (METHOD = 1), i.e. exchanges only one point,
%   or the second method (METHOD = 2), i.e. exchanges all the reference points.
%   XK is a column vector with the reference, H is the levelled error, P is the
%   trigonometric polynomial.
%
%   [XK, NORME, E_HANDLE, FLAG] = EXCHANGE(...) returns the modified reference
%   XK, the supremum norm of the error NORME (included as an output argument,
%   since it is readily computed in EXCHANGE and is used later in TRIGREMEZ), a
%   function handle E_HANDLE for the error, and a FLAG indicating whether there
%   were at least Npts alternating extrema of the error to form the next
%   reference (FLAG = 1) or not (FLAG = 0).
%
%   [XK, ...] = EXCHANGE([], 0, METHOD, F, P, Q, N) returns a grid of N
%   points XK where the error F - P alternates in sign (but not necessarily
%   equioscillates). This feature of EXCHANGE is useful to start TRIGREMEZ from an
%   initial trial function rather than an initial trial reference.

% Compute extrema of the error.
e_num = diff(f - p);
rts = roots(e_num, 'nobreaks');
% Do not include the other end of the domain, since the domain
% is assumed to be periodic:
rr = [f.domain(1) ; rts];

% Function handle output for evaluating the error.
err_handle = @(x) feval(f, x) - feval(p, x);

% Select exchange method.
if ( method == 1 )                             % One-point exchange.
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

% Of the points we kept, choose Npts consecutive ones that include the maximum
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
