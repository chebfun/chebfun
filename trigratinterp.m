function varargout = trigratinterp(fk, m, n, varargin)
%RATINTERP  Robust trigonometric rational interpolation or least-squares.
%   [P, Q, R_HANDLE] = TRIGRATINTERP(F, M, N) computes the (M, N) 
%   trigonometric rational interpolant of F on 2(M+N)+1 
%   equidistant points. F can be a CHEBFUN, a function handle or a column 
%   vector of 2(M+N)+1 data points.  If F is a CHEBFUN, the rational 
%   interpolant is constructed on the domain of F. Otherwise, the domain 
%   [-1, 1] is used. P and Q are periodic CHEBFUNs such that 
%   P(x)./Q(x) = F(x). R_HANDLE is a function handle evaluating the 
%   rational interpolant directly.
%
%   [P, Q, R_HANDLE] = TRIGRATINTERP(F, M, N, NN) computes a (M, N) 
%   trigonometric rational linear least-squares approximant of F over NN 
%   equidistant points. If NN = 2(M+N)+1 or NN = [], a rational interpolant 
%   is computed.
%
%   [P, Q, R_HANDLE] = TRIGRATINTERP(F, M, N, NN, XI) computes a (M, N) 
%   trigonometric rational interpolant or approximant of F over the vector of 
%   nodes XI. NN is ignored in this case. However, XI can also be the 
%   strings 'equi' or 'equidistant', in which case NN equidistant nodes are 
%   created on the interval [-1, 1) using TRIGPTS.
%
%   [P, Q, R_HANDLE, MU, NU] = TRIGRATINTERP(F, M, N, [], [], TOL)
%   computes a robustified (M, N) trigonometric rational interpolant or 
%   approximant of F, in which components contributing less than 
%   the relative tolerance TOL to the solution are discarded. If no value 
%   of TOL is specified, a tolerance of 1e-14 is assumed; set TOL to zero 
%   to disable robustness. MU and NU are the resulting numerator and 
%   denominator degrees. Note that if the degree is decreased, a rational 
%   least-squares approximation is computed over the 2(M+N)+1 points.
%
%   [P, Q, R_HANDLE, MU, NU, POLES, RES] = TRIGRATINTERP(F, M, N, NN, XI)
%   returns the poles POLES of the rational interpolant on the real axis as
%   well as the residues RES at those points.
%
%   Examples:
%
%   Compute a type (5, 5) robustified rational interpolant to 1/(sin(pi*x) -
%   0.2) on [-1, 1] in equispaced nodes:
%
%     [p, q, r] = trigratinterp(@(x) 1./(sin(pi*x) - 0.2), 5, 5);
%
%   Same thing but with robustness disabled:
%
%     [p, q, r] = trigratinterp(@(x) 1./(sin(pi*x) - 0.2), 5, 5, [], [], 0);
%
%   References:
%
%   [1] Mohsin Javed, "ALGORITHMS FOR TRIGONOMETRIC POLYNOMIAL AND RATIONAL 
%   APPROXIMATION" DPhil thesis. 
%
%
% See also RATINTERP, INTERP1, CHEBPADE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO:  Deal with array-valued CHEBFUNs / quasimatrices.

% Parse the inputs.
[dom, fk, m, n, ~, th, th_type, robustness_flag, ...
    interpolation_flag, tol] = ...
    parseInputs(fk, m, n, varargin{:});

% Set up some values which we will use often.
ts = tol*norm(fk, inf);

% Check for symmetries.
[fEven, fOdd] = checkSymmetries(fk, th, th_type, ts);

% Form matrices for the linear system for the coefficients.
[ac, bc, ~] = trig_rat_interp(fk, m, n, th, fEven, fOdd, robustness_flag, interpolation_flag, ts);

% Get the exact numerator and denominator degrees.
mu = (length(ac) - 1)/2;
nu = (length(bc) - 1)/2;

% Build the numerator and denominator polynomials and create the output
% function handle for evaluating the rational approximation.
[p, q, r] = constructTrigRatApprox(ac, bc, dom, ts);

% Normalize q to 1 if a type (m, 0) approximation is computed
if ( n == 0 || nu == 0 )
    p = p./q;
    q = q./q;
end

% Compute poles and residues if requested.
poles = [];
residues = [];
if ( nargout > 5 )
    if ( nargout > 6 ) % Compute residues.
        % Compute partial fraction expansion of r.
        [residues, poles] = residue(p, q);
        [poles, ind] = sort(poles);
        residues = residues(ind);
        
        % Residues are the coefficients of 1/(x - poles(j)):
        for j = 1:(length(poles) - 1)
            if ( poles(j+1) == poles(j) )
                residues(j+1) = residues(j);
            end
        end
    else               % Just compute the poles.
        poles = roots(q, 'all');
    end
end

outArgs = {p, q, r, mu, nu, poles, residues};
% Return the output based on nargout:
if ( nargout <= 1  )
    varargout{1} = r;
elseif ( nargout <= 7  )
    [varargout{1:nargout}] = outArgs{1:nargout};
else
    error('CHEBFUN:CHEBFUN:ratinterp:nargout', ...
        'Incorrect number of output arguments.');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [dom, fk, m, n, NN, xi, xi_type, robustness_flag, interpolation_flag, tol] = parseInputs(f, m, n, varargin)

% Make sure we have the correct number of arguments.
if ( nargin < 3 )
    error('CHEBFUN:trigratinterp:tooFewArgs', 'Not enough input arguments.');
end

if ( m < 0 || n < 0 || ( m ~= round(m) ) || ( n ~= round(n) ) )
    error('CHEBFUN:trigratinterp:integers', 'M and N must be non-negative integers.');
end


%% Deal with various inputs

% Default domains:
if ( isfloat(f) || isa(f, 'function_handle') )
    dom = [-1, 1];
elseif ( isa(f, 'chebfun') )
    dom = f.domain([1, end]);
else
    error('CHEBFUN:trigratinterp:parseInput', 'Can not recognize the input function F.');
end

% Check if domain is passed:
if ( ~isempty(varargin) )
    domStrIndex = find(strcmpi('dom', varargin));
    if ( ~isempty(domStrIndex) )
        dom = varargin{domStrIndex + 1};
        if ( isa(f, 'chebfun') && any(dom ~= f.domain([1, end])) )
            error('CHEBFUN:trigratinterp:npoints', 'F has different domain from the one passed');
        end
        varargin([domStrIndex, domStrIndex + 1]) = [];
    end
end


% Now process the rest of the varargin
% varargin can only have NN, xi, tol now

% First deal with NN
if ( ~isempty(varargin) )
    NN = varargin{1};
    if ( isfloat(NN) && ( length(NN) == 1 ) )
        if ( rem(NN, 2) == 0 )
            warning('CHEBFUN:trigratinterp:npoints', 'Number of points should be odd.');
        end
        if ( NN < 2*(m+n) + 1 )
            error('CHEBFUN:trigratinterp:npoints', 'NN must be >= 2*(M+N)+1.');
        end
    elseif ( isempty(NN) )
        NN = 2*(m+n)+1;
    else
        error('CHEBFUN:trigratinterp:npoints', 'Can not recognize NN.');
    end
    varargin(1) = [];
else
    % Number of points not provided, use default:
    NN = 2*(m+n)+1;
end


% Now deal with xi or xi_type
% Check if a string specifying the types of points is passed:
if ( ~isempty(varargin) )
    xi = varargin{1};
    if ( strncmpi('equi', xi, 4) || isempty(xi) )
        xi_type = 'equi';
        xi = trigpts(NN, dom);
    elseif ( isfloat(xi) && ( min(size(xi)) == 1 ) )
        % Overwrite the length of points NN:
        NN = length(xi);
        % points are arbitrary:
        xi_type = 'arbi';
    else
        error('CHEBFUN:trigratinterp:xi', 'Can not parse points XI.');
    end
    varargin(1) = [];
else
    xi_type = 'equi';
    xi = trigpts(NN, dom);
end

if ( ~isempty(varargin) )
    % tol should be 0 for no robustness:
    tol = varargin{1};
    if ( all(tol < 0) || all(size(tol) ~= [1, 1]) )
        error('CHEBFUN:trigratinterp:tol', 'tol must be a positive number.');
    end
else
    % default robustness tolerance:
    tol = 1e-14;
end

% Make sure domain is valid:
if ( ~isfloat(dom) || ~isequal(size(dom), [1 2]) )
    error('CHEBFUN:trigratinterp:badDom1', ...
        'Domain should be a 1 x 2 row vector of endpoints.');
end

if ( diff(dom) <= 0 )
    error('CHEBFUN:trigratinterp:badDom2', 'Invalid domain.');
end

% Sort the nodes and make sure they are within the domain:
xi = sort(xi);
if ( ( min(xi) < dom(1) ) || ( max(xi) > dom(2) ) )
    error('CHEBFUN:trigratinterp:domXI', ...
        'Input vector XI must be within the domain.');
end

if ( rem(length(xi), 2) == 0 )
    warning('CHEBFUN:trigratinterp:evenNumberOfPoints', ...
        'Input vector XI does not have odd number of points.');
end

if ( ~isreal(xi) )
    error('CHEBFUN:trigratinterp:realXI', ...
        'Input vector XI must be real.');
end

if ( isfloat(f) )
    fk = f;
else
    % If we were given a function handle or CHEBFUN, sample it on the grid.
    fk = f(xi);
end

if ( length(fk) ~= NN )
    error( 'CHEBFUN:trigratinterp:lengthF', ...
        'Input vector F does not have the correct length.');
end

% Scale nodes to [-1 1].
xi = 2.0 * (xi - 0.5*sum(dom)) / diff(dom);
if ( ( min(xi) == -1 ) && ( max(xi) == 1 ) )
    error('CHEBFUN:trigratinterp:duplicate', ...
        'Periodic interval cannot have both -1 and 1 as points of interpolation.');
end

% Check if robustness is disabled:
if ( tol == 0 )
    % Turn off robustness if tolerance is zero:
    robustness_flag = false;
else
    robustness_flag = true;
end

% Check if we are interpolating or doing least-squares:
interpolation_flag = false;
if ( NN == 2*(m+n)+1 )
    interpolation_flag = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fEven, fOdd] = checkSymmetries(f, xi, xi_type, tol)
% Check symmetries for trigonometric functions.
fEven = false;
fOdd = false;
N = length(xi);
if ( strncmpi(xi_type, 'equi', 4) )
    % The point at the left end i.e. f(-1) is left out while
    % checking for symmetry.
    if ( mod(N, 2) == 1 )
        M = floor(N/2);
        fl = f(2:M+1);
        fr = f(end:-1:M+2);
    else
        M = N/2;
        fl = f(2:M+1);
        fr = f(end:-1:M+1);
    end
    fEven = norm(fl - fr, inf) < tol;
    fOdd = norm(fl + fr, inf) < tol;
else
    % Other nodes.
    [xi, idx] = sort(xi);
    f = f(idx);
    if ( xi(1) == -1 )
        if ( mod(N, 2) == 1 )
            M = floor(N/2);
            xl = xi(2:M+1);
            xr = xi(end:-1:M+2);
            fl = f(2:M+1);
            fr = f(end:-1:M+2);
        else
            M = N/2;
            xl = xi(2:M+1);
            xr = xi(end:-1:M+1);
            fl = f(2:M+1);
            fr = f(end:-1:M+1);
        end
    elseif ( xi(end) == 1 )
        if ( mod(N, 2) == 1 )
            M = floor(N/2);
            xl = xi(1:M);
            xr = xi(end-1:-1:M+1);
            fl = f(1:M);
            fr = f(end-1:-1:M+1);
        else
            M = N/2;
            xl = xi(1:M);
            xr = xi(end-1:-1:M);
            fl = f(1:M);
            fr = f(end-1:-1:M);
        end
    else
        if ( mod(N, 2) == 1)
            M = ceil(N/2);
            xl = xi(1:M);
            xr = xi(end:-1:M);
            fl = f(1:M);
            fr = f(end:-1:M);
        else
            M = N/2;
            xl = xi(1:M);
            xr = xi(end:-1:M+1);
            fl = f(1:M);
            fr = f(end:-1:M+1);
        end
    end
    
    if ( norm(xl + xr, inf) < tol )
        % The nodes are symmetric, check values:
        fEven = norm(fl - fr, inf) < tol;
        fOdd = norm(fl + fr, inf) < tol;
    end
end

end

function [P, Q] = construct_matrices(th, m, n, T)
% Input: th has the points of interpolation, m, n
%        are degrees of p and q and T is the period
% Output P: matrix of dim length(th) x 2*m+1
% Output Q: matrix of dim length(th) x 2*n+1

% Initialize output matrices
P = zeros(length(th), 2*m+1);
Q = zeros(length(th), 2*n+1);

% Construct P:
P(:, 1) = ones(length(th), 1);
for j = 1:m
    P(:, 2*j)   = sin(2*j*pi/T*th);
    P(:, 2*j+1) = cos(2*j*pi/T*th);
end

% Construct Q:
Q(:, 1) = ones(length(th), 1);
for j = 1:n
    Q(:, 2*j)   = sin(2*j*pi/T*th);
    Q(:, 2*j+1) = cos(2*j*pi/T*th);
end

end

function [ac, bc] = getCoeffs(V, m, n, fEven, fOdd)
% Input: U, S, V is the SVD decomposition of the matrix
% of the (linearized) interpolation problem. m, n are degrees
% of p and q and fEven and fOdd are flags indicating symmetries
% Output: The coefficients of p and q are returned in
% the complex exponential basis

if ( fEven || fOdd )
    if ( fEven )
        a = V(1:m+1, end);
        b = V(m+2:end, end);
        tmp = zeros(2*m+1, 1);
        tmp(1:2:end) = a;
        a = tmp;
        tmp = zeros(2*n+1, 1);
        tmp(1:2:end) = b;
        b = tmp;
    elseif( fOdd )
        a = V(1:m, end);
        b = V(m+1:end, end);
        tmp = zeros(2*m+1, 1);
        tmp(2:2:end) = a;
        a = tmp;
        tmp = zeros(2*n+1, 1);
        tmp(1:2:end) = b;
        b = tmp;
    end
else
    a = V(1:2*m+1, end);
    b = V(2*m+2:end, end);
end


ac = sincosine_to_exponential(a);
bc = sincosine_to_exponential(b);

end

function [ac, bc, s] = trig_rat_interp(fk, m, n, th, fEven, fOdd, robustness_flag, interpolation_flag, tol)
stdDomain = [-1, 1];
fEven = 0;
fOdd = 0;
T = stdDomain(end) - stdDomain(1);
while (1)
    % Set up the matrices:
    [P, Q] = construct_matrices(th, m, n, T);
    N = m + n;
    if ( fEven && interpolation_flag )
        th = th(1:N+1);
        fk = fk(1:N+1);
        P = P(1:N+1, 1:2:end); % only cosines
        Q = Q(1:N+1, 1:2:end); % only cosines
    elseif( fOdd && interpolation_flag )
        th = th(1:N);
        fk = fk(1:N);
        P = P(1:N, 2:2:end); % only sines
        Q = Q(1:N, 1:2:end); % only cosines
    end
    
    % Diagonal matrix of function values:
    D = diag(fk);
    
    % Scale the matrix whose null space we want to find:
    sysMat = [P, -D*Q];
    sysMat = diag(1./max(abs(fk), 1)) * sysMat;

    % Apply SVD to find the null space:
    [~, S, V] = svd(sysMat);
    
    % Get the coefficients of p and q in complex exponential basis:
    [ac, bc] = getCoeffs(V, m, n, fEven, fOdd);
    
    % Chop small coefficients at ends:
    ac = chopCoeffs(ac, tol);
    bc = chopCoeffs(bc, tol);
    
    % Singular values of the problem:
    s = diag(S);
    
    % Degree reduction based on small singular values.
    % The code exits the loop if there are less than
    % 2 small singular values.
    
    m = min([(length(ac)-1)/2, m]);
    
    n_big_sing_vals = find(abs(s) > tol, 1, 'last');
    n_small_sing_vals = length(s) - n_big_sing_vals;
    if ( fEven || fOdd )
        if ( ( n_small_sing_vals < 2 ) || ~robustness_flag )
            break
        else
            reduction = n_small_sing_vals;
        end
    elseif ( ( n_small_sing_vals < 2 ) || ~robustness_flag )
        break
    else
        reduction = floor(n_small_sing_vals/2);
    end
    n_new = n - reduction;
    if ( n_new >= 0 )
        n = n_new;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for assembling the rational interpolant.

function [p, q, r] = constructTrigRatApprox(ac, bc, dom, tol)

% Construct chebfuns on the standard domain
stdDomain = [-1, 1];
p = chebfun(ac, 'coeffs', 'trig');
q = chebfun(bc, 'coeffs', 'trig');


% Map to new domain if needed
if ( any(dom ~= stdDomain) )
    p = newDomain(p, dom);
    q = newDomain(q, dom);
end

% Return the function handle:
% [TODO]: this should use trigbary instead
r = @(x) p(x)./q(x);

end

function b = chopCoeffs(a, tol)
% Input: Assume double sided coeffs in vector a with
%       wavenumbers -k to k where k = (length(a)-1)/2
% Output: b is a vector such that at least one of
%       b_k and b{-k} are greather than tol in abs value.
%

% Handle the empty case:
if ( isempty(a) )
    b = a;
    return
end

% The trivial case of a constant:
n = length(a);
if ( n <= 1 )
    b = a;
    return
end

if ( mod(n,2) == 0 )
    error('CHEBFUN:trigratinterp:chopCoeffs', 'coefficients must be odd in length');
end

% Chop the tails:
midIdx = (n + 1)/2;
a1 = abs(a(midIdx:end));
a2 = abs(a(midIdx:-1:1));
aa = (a1 + a2)/2;
idx = find(abs(aa) > tol, 1, 'last');
b = a(midIdx-idx+1:midIdx+idx-1);

end

function c = sincosine_to_exponential(a)
% Input vector a has a_1 + a_2 sin + a_3 cos + a_4 sin + ...
% Output vector has ... + c_-1exp() + c_0 + c_1 exp() + ...
tmp = (a(3:2:end) - 1i*a(2:2:end))/2;
c = [flipud(conj(tmp)); a(1); tmp];
end
