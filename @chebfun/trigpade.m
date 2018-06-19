function varargout = trigpade(f, m, n, varargin)
%TRIGPADE   Trigonometric (Fourier) Pade approximation.
%   [P, Q, R_HANDLE] = TRIGPADE(F, M, N) computes trigonometric polynomials 
%   P and Q of degree M and N, respectively, such that the trigonometric 
%   rational function P/Q is the type (M,N) Fourier-Pade approximation of 
%   the periodic CHEBFUN F. That is, the Fourier series of P/Q coincides 
%   with that of the CHEBFUN F up to the maximum possible order for the 
%   polynomial degrees permitted. R_HANDLE is a function handle for 
%   evaluating the trigonometric rational function P/Q.
% 
%   [P, Q, R_HANDLE, TN_P, TD_P, TN_M, TD_M] = TRIGPADE(F, M, N) also
%   returns the four trigonometric polynomials TN_P, TD_P, TN_M and TD_M
%   such that P/Q = TN_P/TD_P + TN_M/TD_M.
%
%   In both of the above cases, if only one output argument is specified
%   then R_HANDLE is returned, while P and Q are returned if two or more 
%   output arguments are specified. 
%
%   Examples:
%
%   Compute a type (1,3) trigonometric Pade approximation of 
%   exp(sin(pi*x)) on [-1, 1]:
% 
%     f = chebfun(@(t) exp(sin(pi*t)), 'trig')      
%     [p, q, r] = trigpade(f, 1, 3);
%     plotcoeffs(f-p./q, '.')
%   
%   Compute a type (4,2) trigonometric Pade approximation of
%   exp(-160*(x-.5).^2) on [-1, 1]:
% 
%     f = chebfun(@(t) exp(-160*(t-.5).^2), 'trig')
%     [p, q, r] = trigpade(f, 4, 2);
%     plotcoeffs(f-p./q, '.')

%   References:
%   [1] G. A. Baker and P. R. Graves-Morris. Pade Approximants,
%   2nd ed., Cambridge University Press, 1996.
%
%   [2] M. Javed. Algorithms for trigonometric polynomial and rational 
%   approximation. DPhil thesis, Oxford, 2017.
%
% See also PADEAPPROX and CHEBPADE

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( isempty(f) )
    varargout{1} = f;
    return
end

if ( ~isPeriodicTech(f) )
    error('CHEBFUN:CHEBFUN:trigpade:trig', ...
        'Input chebfun F must have a periodic representation');
end

dom = f.domain([1, end]);

% Extract the Fourier coefficients of F:
coeffs  = trigcoeffs(f);
N = (length(coeffs)-1)/2;
if ( N ~= round(N) )
    error('c must have odd length')
end

% set tolerance:
tol = 100*eps*norm(coeffs, inf);

% Handle the trivial case:
if ( n == 0 )
    if ( m >= N )
        p = f;
    else
        p = chebfun(coeffs(N+1-m:N+1+m), dom, 'trig');
    end
    q = chebfun(1, dom, 'trig');
    r_h = @(x) feval(p, x);
    tn_p = p/2;
    tn_m = p/2;
    td_p = q;
    td_m = q;
else
    [p, q, r_h, tn_p, td_p, tn_m, td_m] = trig_pade(coeffs, m, n, N, dom, tol);
end

% Form the outputs
outArgs = {p, q, r_h, tn_p, td_p, tn_m, td_m};
if ( nargout <= 1 )
    varargout{1} = r_h;
elseif ( nargout <= 7 )
    [varargout{1:nargout}] = outArgs{1:nargout};
else
    error('CHEBFUN:CHEBFUN:trigpade:nargout', ...
        'Incorrect number of output arguments.'); 
end
    

end

%%
function [p, q, r_h, tn_p, td_p, tn_m, td_m] = trig_pade(coeffs, m, n, N, dom, tol)
% Make sure that we have sufficiently long
% Laurent series to avoid indexing issues.
d = 2*max([m, n]) - N;
if ( d > 0 )
    % Pad with zeros if needed:
    coeffs = [zeros(d, 1); coeffs; zeros(d, 1)];    
end

[ap, bp, am, bm] = laurent_pade(coeffs, m, n);

%% Construct the four trigonometric polynomials of Fourier Pade approximaiton

% _d is for denominator, _n is for numerator
% Denominator and numerator for the positive part
tn_p = chebfun(ap, dom, 'coeffs', 'trig' );
td_p = chebfun(bp, dom, 'coeffs', 'trig' );


% denominator and numerator for the negative part
tn_m = chebfun(am, dom, 'coeffs', 'trig' );
td_m = chebfun(bm, dom, 'coeffs', 'trig' );


L = (max([length(am),length(ap), length(bm), length(bp)])-1)/2;
am = [zeros(L-(length(am)-1)/2, 1); am; zeros(L-(length(am)-1)/2, 1)];
bm = [zeros(L-(length(bm)-1)/2, 1); bm; zeros(L-(length(bm)-1)/2, 1)];
ap = [zeros(L-(length(ap)-1)/2, 1); ap; zeros(L-(length(ap)-1)/2, 1)];
bp = [zeros(L-(length(bp)-1)/2, 1); bp; zeros(L-(length(bp)-1)/2, 1)];

if ( ~all([length(am),length(ap), length(bm), length(bp)] == 2*L+1) ) 
    error('all vectors of coefficients should be of same length')
end

pk = conv(ap, bm) + conv(am, bp);
qk = conv(bm, bp);
% discard small coeffs:
pk = chop_coeffs(pk, tol);
qk = chop_coeffs(qk, tol);

p = chebfun(pk, dom, 'coeffs', 'trig');
q = chebfun(qk, dom, 'coeffs', 'trig');
p = simplify(p);
q = simplify(q);
% Construct the full approximation:
%p = tn_p.*td_m + tn_m.*td_p; 
%q = td_p.*td_m;
r_h = @(t) feval(p, t)./feval(q, t);

% Discard the imaginary rounding errors:
if ( norm(coeffs - conj(flipud(coeffs)), inf) < tol ) 
    % The input function is real, check the imaginary part
    % of the approximation
    if ( norm(imag(p./q)) > tol )        
        warning('CHEBFUN:CHEBFUN:trigpade:imag', 'imaginary part not negligible.');
    else
        p = real(p);
        q = real(q);
        r_h = @(t) feval(p, t)./feval(q, t);
    end
end

end

%%
function [ap, bp, am, bm] = laurent_pade(c, m, n, tol)
%   Compute the Laurent-Pade approximation from the 
%   given Laurent coefficients in c. 
%   Inputs: It is assumed that c has the form 
%   [c_{-N}, ..., c_{-1}, c_0, c_1, ..., c_{N}].'
%   In particular, c has odd length 2N+1
%   m and n are non-negative integers and tol is 
%   a tolerance.
%
%   Outputs are the coefficients of Laurent fractions
%   The coefficients for a_p/b_p form p+/q+
%   a_m/b_m form p-/q-
%   The output vectors are ready to be used by 
%   the chebfun trig constructor, for example
%   ap is of the form [..., 0, 0, a_0, a_1, a_2, ...]
%   while am has the form 
%   [ ..., a_{-2}, a_{-1}, a_0, 0, 0, ...]. 
%   Zero padding is also done on each side of the
%   output vectors.
       

% Set default tolerance
if ( nargin < 4 )
    tol = 100*eps*norm(c, inf);
end

% make sure c is a column vector
c = c(:);

% Middle of the Laurent series, location of c_{-1}
N = (length(c)-1)/2;

% Handle the special case:
if ( n == 0 )
    bp = 1;
    bm = 1;
    ap = c(N+1:N+1+m);
    am = c(N+1:-1:N+1-m);
    ap(1) = ap(1)/2;
    am(1) = am(1)/2;
    ap = [zeros(length(ap)-1, 1); ap];
    am = [zeros(length(am)-1, 1); am];
    return
end

% Compute one sided Laurent approximation
[ap, bp] = laurent_approx(c, m, n, N, tol);
c_rev = flipud(c);

if ( norm(c-conj(c_rev), inf) < 10*tol )
    % Function is real, don't solve the 
    % same problem again
    am = conj(ap);
    bm = conj(bp);
else
    % Compute the second Laurent approximation
    [am, bm] = laurent_approx(c_rev, m, n, N, tol);
end

% Pad coefficients with zeros
ap = [zeros(length(ap)-1, 1); ap];
bp = [zeros(length(bp)-1, 1); bp];


% Pad coefficients with zeros
am = [zeros(length(am)-1, 1); am];
am = flipud(am);
bm = [zeros(length(bm)-1, 1); bm];
bm = flipud(bm);

end

function [a, b] = laurent_approx(c, m, n, N, tol)

% Input: It is assumed that the Laurent
% coefficients are stored in an array 
% [c_{-N}, ..., c_{-1}, c_0, c_1, ..., c{N}]
% and the length of the array is 2N + 1.
% To convert a Laurent index to Matlab index
% an offset of N+1 is needed.
% c_0 is indexed at N + 1
% c_k is indexed at k + N + 1
% c_{-N} is indexed at 1 and so on
% 
% m, n are non-negative integers
% and tol is a tolerance for zero detection
%
% Output: a and b are vectors with coefficients 
% for a Pade approximation of the power series with
% coefficients c_0, c_1, ...


% Each column has n rows
col = c((m+1)+(N+1):(m+n)+(N+1));

% Each row has n+1 columns
row = c((m+1)+(N+1):-1:(m+1-n)+(N+1));

% C is of size n x (n+1)
C = toeplitz(col, row);

% It has at least one null vector, which 
% forms the coefficients of the denominator
b = null(C);
if ( size(b, 2) > 1)
    b = b(:, 1);
end

% Normalize b, the derivation following Baker and Graves-Morris
% requires this normalization 
if ( abs(b(1)) < tol )
    error('CHEBFUN:TRIGPADE:laurent_approx', 'denominator zero at the origin detected')
else
    b = b/b(1);
end

% Compute the coefficients of the numerator
M = max([m, n]);
col = c(0+(N+1):M+(N+1));
col(1) = col(1)/2;
% C is an (M+1) X (M+1) square matrix
% Note: it is crucial to say toeplitz(col, col) to 
% avoid hermitian formation of the following matrix
C = toeplitz(col, col);
C = tril(C);
% Make sure bb has length M+1
bb = [b; zeros(M+1-length(b), 1)];
a = C*bb;

end


function a = chop_coeffs(c, tol)
% Discard small coeffs from either 
% end of c using tol, preserving 
% the maximum Fourier mode if asymmetric
mid = (length(c)+1)/2;
nm = mid - find(abs(c)>tol, 1, 'first');
np = find(abs(c)>tol, 1, 'last') - mid;
nn = max([nm, np]);
a = c(mid-nn:mid+nn);
end
