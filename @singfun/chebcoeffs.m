function out = chebcoeffs(f, N)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a SINGFUN.
%   A = CHEBCOEFFS(F, N) returns the row vector of first N Chebyshev coefficients
%   of F such that F = ... + A(1) T_{N-1}(x) + ... + A(N-1) T_1(x) + A(N) T_0(x)
%   where T_k(x) denotes the k-th Chebyshev polynomial.
%
% See also LEGPOLY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: This requires a test.

if ( any(f.exponents <= -.5 ) )
    error('CHEBFUN:SINGFUN:chebcoeffs:notIntegrable', ...
        'F does not have a well-defined Chebyshev expansion.');
end

if ( (nargin == 1) || isempty(N) )
    error('CHEBFUN:SINGFUN:chebcoeffs:input', ...
        'F does not have a finite Chebyshev series. Please input N.');
end

if ( isa(f.smoothPart, 'chebtech') )
    % Compute the required inner products in coefficient space efficiently.
    
    n = length(f);
    % Compute the Chebyshev moments:
    w = computeWeights(f.exponents - [.5 .5], n + N - 1).';
    % The Chebyshev coefficients of f:
    b = .5*flipud([zeros(N-1, 1) ; get(f, 'coeffs')]);
    % Multiplication matrix for coefficients of f: (fast, using FFTs)
    out = flipud(fastToeplitzMult(b, w) + fastHankelMult(b, w) - w(1)*b);
    % Trim to required length and scale T_0 term:
    out = [out(end-N+1:end-1) ; out(end)/2].';
    
else
    % Compute the required inner products by calling SUM().
   
    for k = N-1:-1:0
        % Make the kth Chebyshev polynomial:
        Tk = f.smoothPart.make(@(x) cos(k*acos(x)));
        % Compute the product of f with T_k:
        g = f.*Tk;
        % Include the Chebyshev weight:
        g.exponents = g.exponents - [.5, .5];
        % Integrate using SUM():
        out(N-k) = sum(g);
    end
    % Scale the final coefficient:
    out(end) = out(end)/2;
    
end

end

function w = computeWeights(exps, n)
% Compute the weights for the Chebyshev-Jacobi moments.

% TODO: De-duplicate this (also in SUM) by making a static method.

% Grab the exponents:
a = exps(1);

if ( diff(exps) == 0 )
    % If the exponents at the endpoints are same, then compute the
    % appropriate modified moments for Gegenbauer weights.
    
    r = a + .5;
    m0 = gamma(r + .5)*sqrt(pi)/gamma(r + 1);
    k = 1:floor((n-1)/2);
    % Even modified moments for M_2k = \int_{-1}^1 (1-x)^a(1+x)^a T_2k(x) dx
    % and notice that the odd moments vanish due to parity.
    m = m0*[1, cumprod((k - r - 1)./(k + r))];
    % Form the modified moments vector:
    w(1:2:n) = m;
    w(2:2:n) = 0;
    
else
    
    % The general case:
    b = exps(2);
    
    % Common coefficients for the modified moments:
    c1 = a + 1;
    c2 = b + 1;
    c3 = a + b + 1;
    c4 = c1 + c2;
    c5 = a - b;
    c0 = (2^c3)*beta(c1, c2);
    
    % Compute the hypergeometric function related to the modified moments:
    w = zeros(1,n);
    w(1) = 1;
    if ( n > 1 )
        w(2) = c5/c4;
        % Sister Celine's three-term recurrence:
        for j = 3:n
            w(j) = (2*c5*w(j-1) + (j - 2 - c4)*w(j-2)) / (c3 + j - 1);
        end
    end
    % Compute the modified moments:
    w = c0*w;
    
end

end

function y = fastToeplitzMult(a, b, x)
% FASTTOEPLITZMULT(a,b,x) computes y = toeplitz(a,b)*x in O(n logn) operations.
% FASTTOEPLITZMULT(a,x) computes y = toeplitz(a)*x in O(n logn)

% TODO: Move this to @CHEBTECH?
a(1) = 2*a(1);
if ( nargin == 2 )
    x = b;
    b = a;
end
n = length(a);
p = ifft(fft([a ; 0 ; flipud(b(2:end))]) .* fft([x ; zeros(n,1)]));
y = p(1:n);
end

function y = fastHankelMult(a, b, x)
% FASTHANKELMULT(a,b,x) computes hankel(a,b) * x in O(nlogn) operations.
% FASTHANKELMULT(a,x) computes hankel(a) * x in O(nlogn) operations.

% TODO: Move this to @CHEBTECH?

if ( nargin == 2 )
    n = length(b);
    x = b;
    b = [a(end) ; zeros(n-1,1)];
end
n = length(x);
p = ifft(fft([b ; 0 ; (a(1:end-1))]) .* fft([x(end:-1:1) ; zeros(n, 1)]));
y = p(1:n);
end
