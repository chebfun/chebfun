function varargout = trigpade(f, m, n, varargin)
%TRIGPADE   Trigonometric (Fourier) Pade approximation.
%   [P, Q, R_HANDLE] = TRIGPADE(F, M, N) computes trigonometric polynomials 
%   P and Q of degree M and N, respectively, such that the trigonometric 
%   rational function P/Q is the type (M,N) Fourier-Pade approximation of 
%   the periodic CHEBFUN F. That is, the Fourier series of P/Q coincides 
%   with that for the CHEBFUN F up to the maximum possible order for the 
%   polynomial degrees permitted. R_HANDLE is a function handle for 
%   evaluating the trigonometric rational function P/Q.
% 
%   [P, Q, R_HANDLE, TN_P, TD_P, TN_M, TD_M] = TRIGPADE(F, M, N) also
%   returns the four trigonometric polynomials TN_P, TD_P, TN_M and TD_M
%   such that P/Q = TN_P./TD_P + TN_M./TD_M.
%
%   In both of the above cases, if only one output argument is specified
%   then R_HANDLE is returned, while P and Q are returned if two or more 
%   output arguments are specified. 
%
%   References:
%
%   [1] G. A. Baker and P. R. Graves-Morris. Pade Approximants. 
%   Cambridge University Press, second edition, 1996.
%
%   [2] M. Javed. Algorithms for trigonometric polynomial and rational 
%   approximation. DPhil thesis, Oxford.
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

% Make sure that we have sufficiently long
% Laurent series to avoid indexing issues.
d = 2*max([m, n]) - N;
if ( d > 0 )
    % Pad with zeros if this is not the case:
    coeffs = [zeros(d, 1); coeffs; zeros(d, 1)];    
end

[a_p, b_p, a_m, b_m] = laurent_pade(coeffs, m, n);

%% Construct the four trigonometric polynomials of Fourier Pade approximaiton

% Pad coefficients with zeros:
aa_p = [zeros(length(a_p)-1, 1); a_p];
bb_p = [zeros(length(b_p)-1, 1); b_p];

% _d is for denominator, _n is for numerator
% Denonminator and numerator for the +ve part
tn_p = chebfun(aa_p, dom, 'coeffs', 'trig' );
td_p = chebfun(bb_p, dom, 'coeffs', 'trig' );

% Pad coefficients with zeros
aa_m = [zeros(length(a_m)-1, 1); a_m];
aa_m = flipud(aa_m);
bb_m = [zeros(length(b_m)-1, 1); b_m];
bb_m = flipud(bb_m);

% denonminator and numerator for the -ve part
tn_m = chebfun(aa_m, dom, 'coeffs', 'trig' );
td_m = chebfun(bb_m, dom, 'coeffs', 'trig' );

% Construct the full approximation:
p = tn_p.*td_m + tn_m.*td_p; 
q = td_p.*td_m;
r_h = @(t) p(t)./q(t);

% Discard the imaginary rounding errors:
tol = 1e-13;
if ( norm(coeffs - conj(flipud(coeffs)), inf) < tol ) 
    % The input function is real, check the imaginary part
    % of the approximation
    if ( norm(imag(p./q)) > tol )        
        warning('CHEBFUN:CHEBFUN:trigpade:imag', 'imaginary part not negligible.');
    else
        p = real(p);
        q = real(q);
        r_h = @(t) p(t)./q(t);
    end
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
function [a_p, b_p, a_m, b_m] = laurent_pade(c, m, n)
%   Compute the Laurent-Pade approximation from the 
%   given Laurent coefficients in c. 
%   Inputs: It is assumed that c has the form 
%   [c{-N}, ..., c_{-1}, c0, c_1, ..., c_{N}].'
%   In particular, c has odd length 2N+1
%   m and n are non-negative integers
%
%   Outputs are the coefficients of Laurent fractions
%   The coefficients for a_p/b_p form p+/q+
%   a_m/b_m form p-/q-
%   The output vectors correspond to polynomials
%   a(1) + a(2)z + ... + a(m) z^m, etc. The ordering
%   is sucht that the lowest order coeff is at the 
%   top of each output vector
       
% make sure c is a column vector
c = c(:);

% Middle of the laurent series, location of c_{-1}
N = (length(c)-1)/2;

% Handle the special case:
if ( n == 0 )
    b_p = 1;
    b_m = 1;
    a_p = c(N+1:N+1+m);
    a_m = c(N+1:-1:N+1-m);
    a_p(1) = a_p(1)/2;
    a_m(1) = a_m(1)/2;
    return
end

% Compute one sided Laurent approximation
[a_m, b_m] = laurent_approx(c, m, n, N);
c_rev = flipud(c);

if ( norm(c-conj(c_rev), inf) < 1e-13 )
    % Function is real, don't solve the 
    % same problem again:
    a_p = conj(a_m);
    b_p = conj(b_m);
else
    % Compute the second Laurent approximation
    [a_p, b_p] = laurent_approx(c_rev, m, n, N);
end

end

function [a, b] = laurent_approx(c, m, n, N)

% Input: It is assumed that the laurent
% coefficients are stored in an array 
% [c_{-N}, ..., c_{-1}, c_0, c_1, ..., c{N}]
% and the length of the array is 2N + 1.
% To convert a laurent index to matlab index
% an offset of N+1 is needed.
% c_0 is indexed at N + 1
% c_k is indexed at k + N + 1
% c_{-N} is indexed 1
% c_{N} is indexed at 2N+1
% 
% m, n are non-negative integers
%
% Output: a and b are vectors with coefficents 
% for a Pade approximation of the negative 
% coefficients within c

% Each column has n rows
col = c((m+1)+(N+1):(m+n)+(N+1));

% Each row has n+1 columns
row = c((m+1)+(N+1):-1:(m+1-n)+(N+1));

% C is of size n x (n+1):
C = toeplitz(col, row);

% It has at least one null vector, which 
% forms the coefficents of the denominator
b = null(C);
if ( size(b, 2) > 1)
    b = b(:, 1);
end

% Compute the coefficients of the numerator
M = max([m, n]);
col = c(0+(N+1):M+(N+1));
col(1) = col(1)/2;
% C is an (M+1) X (M+1) square matrix
C = toeplitz(col);
C = tril(C);
% Make sure bb has length M+1
bb = [b; zeros(M+1-length(b), 1)];
a = C*bb;

end