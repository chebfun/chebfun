function [r, a, b, mu, nu, poles, residues] = padeapprox(f, m, n, tol)
%PADEAPPROX   Pade approximation to a function or Taylor series.
%   [R, A, B, MU, NU, POLES, RESIDUES] = PADEAPPROX(F, M, N, TOL) constructs a
%   Pade approximant to F using the robust algorithm from [1] based on the SVD.
%   F must be a function handle or a vector of coefficients f_0, ..., f_{m + n},
%   and if F is a function handle, the function must be analytic in a
%   neighborhood of the unit disc, since the coefficients are computed via FFT.
%   M and N are the desired numerator and denominator degrees, respectively, and
%   must be nonnegative. The optional TOL argument specifies the relative
%   tolerance; if omitted, it defaults to 1e-14. Set TOL to 0 to turn off
%   robustness. The output is a function handle R of for an exact type (MU, NU)
%   Pade approximant to F with coefficient vectors A and B and, optionally, the
%   POLES and RESIDUES.
%
%   This code is included in the Chebfun distribution for the convenience of
%   readers of _Approximation Theory and Approximation Practice_, but it is not
%   actually a Chebfun code. A Chebfun analogue is CHEBPADE.
%
%   References:
%
%   [1] P. Gonnet, S. Guettel, and L. N. Trefethen, "ROBUST PADE APPROXIMATION 
%       VIA SVD", SIAM Rev., 55:101-117, 2013.
%
% See also CHEBPADE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to relative tolerance of 1e-14.
if ( nargin < 4 )
    tol = 1e-14;
end

% Compute coefficients if necessary.
if ( ~isnumeric(f) )
    % Sample at many roots of unity and use FFT to get coeffs.
    N = 2048;
    z = exp(2i*pi*(0:N-1)'/N);
    f = fft(f(z))/N;

    % Discard near-zero coefficients.
    tc = 1e-15*norm(f);
    f(abs(f) < tc) = 0;

    % Remove imaginary rounding errors.  (Make real functions real.)
    if ( norm(imag(f), inf) < tc )
        f = real(f);
    end
end

% Make sure c is long enough but not longer than necessary.
c = [f(:) ; zeros(m + n + 1 - length(f), 1)];
c = c(1:m+n+1);

% Compute absolute tolerance.
ts = tol*norm(c);

% Compute the Pade approximation.
if ( norm(c(1:m+1), inf) <= tol*norm(c, inf) ) % Special case r = 0.
    a = 0;
    b = 1;
    mu = -inf;
    nu = 0;
else                                           % General case.
    % First row/column of Toeplitz matrix.
    row = [c(1) zeros(1,n)];
    col = c;

    % Do diagonal hopping across block.
    while ( true )
        % Special case n == 0.
        if (n == 0)
            a = c(1:m+1);
            b = 1;
            break
        end

        % Form Toeplitz matrix.
        Z = toeplitz(col(1:m+n+1), row(1:n+1));

        % Compute numerical rank.
        C = Z(m+2:m+n+1,:);
        rho = sum(svd(C) > ts);

        if ( rho == n)
            break
        end

        % Decrease mn, n if rank-deficient.
        m = m - (n - rho);
        n = rho;
    end

    % Hopping finished. Now compute b and a.
    if ( n > 0 )
        [U, S, V] = svd(C,0);

        % Null vector gives b.
        b = V(:,n+1);

        % Do final computation via reweighted QR for better zero preservation.
        D = diag(abs(b) + sqrt(eps));
        [Q, R] = qr((C*D).');

        % Compensate for reweighting.
        b = D*Q(:,n+1);
        b = b/norm(b);

        % Multiplying gives a.
        a = Z(1:m+1,1:n+1)*b;

        % Count leading zeros of b.
        lam = find(abs(b) > tol, 1, 'first') - 1;

        % Discard leading zeros of b and a.
        b = b(lam+1:end);
        a = a(lam+1:end);

        % Discard trailing zeros of b.
        b = b(1:find(abs(b) > tol, 1, 'last'));
    end

    % Discard trailing zero coefficients in a.
    a = a(1:find(abs(a) > ts, 1, 'last'));

    % Normalize.
    a = a/b(1);
    b = b/b(1);

    % Exact numerator, denominator degrees.
    mu = length(a) - 1;
    nu = length(b) - 1;
end

% Function handle for r.
r = @(z) polyval(a(end:-1:1), z)./polyval(b(end:-1:1), z);

if ( nargout > 5 )
    % Compute poles if requested.
    poles = roots(b(end:-1:1));

    % Estimate residues.
    t = max(tol, 1e-7);  % Perturbation for residue estimate.
    residues = t*(r(poles + t) - r(poles - t))/2;
end

end
