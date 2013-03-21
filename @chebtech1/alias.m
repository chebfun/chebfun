function coeffs = alias(coeffs, m)
%ALIAS   Alias Chebyshev coefficients on the 1st kind Chebyshev grid.
%   ALIAS(C, M) aliases the Chebyshev coefficients stored in the column
%   vector C to have length M. If M > LENGTH(C), the coefficients are padded
%   with zeros. If C is a matrix of coefficients, each of the columns is aliased
%   to length M.

% [TODO]: Simplify these comments and include an appropriate reference
% if we can find one.

% Note that aliasing for the 1st kind Chebyshev grid is different from its
% counterpart for 2nd kind Chebyshev grid, though they both
% live in the coefficient space with respect to Chebyshev polynomials of 
% 1st kind. For more information about aliasing and de-aliasing on the 
% Chebyshev points of 1st kind, please read on.

% Aliasing of Chebyshev polynomials of the first kind T_m on Chebyshev 
% points of the first kind x_k = cos(pi*(2k+1)/(2(n+1))
%
% The following statements are adapted from Theorem 4.1 in ATAP.
%
% (1) For any n >= 1 and 0 <= m <= n, the following Chebyshev polynomials
% take the same values on the (n+1)-point Chebyshev grid:
%
% T_m, -T_{2n+2-m}, -T_{2n+2+m}, T_{4n+4-m}, T_{4n+4+m}, -T_{6n+2-m}, 
% -T_{6n+6+m}, ...
%
% (2) Equivalently, for any k >= 0, T_k takes the same value on the grid
% as (-1)^p*T_m with
%
% m = abs( mod((k+n), 2(n+1)) - n ),                           (I)
%
% a number in the range 0 <= m <= n. Here, p = floor( (n+1+m)/(2(n+1)) ).
%
% Proof:
%
% (1) The first assertion is clear if we notice that 
% T_m(x_k) = cos( m*pi*(2k+1)/(2(n+1))
%          = (-1)^p * cos( (2*p*(n+1) (+/-) m) * ( pi*(2k+1)/(2(n+1)) ) )
%          = (-1)^p T_{2p(n+1) (+/-) m} (x_k)
%
% (2) Suppose first that 0 <= mod(k, 2(n+1)) <= n+1. Then 
% n <= mod(k+n, 2(n+1)) <= 2n+1, so Eqn (I) reduces to m = mod(k, 2(n+1)), 
% with 0 <= m <= n+1, and we have just shown that this implies that 
% (-1)^p*T_k and T_m take the same values on the grid. On the other hand,
% suppose that n+2 <= mod(k, 2(n+1)) <= 2n+1. Then 
% 0 <= mod(k+n, 2(n+1)) <= n-1, so the absolute value becomes a negation 
% and Eqn (I) reduces to m = 2(n+1)-mod(k, 2(n+1)), with 1 <= m <= n. Again
% we have just shown that this implies that (-1)^p*T_k and T_m take the 
% same values on the grid.
%
% Remark: Note that on a (n+1)-point Chebyshev grid of 1st kind, T_{j(n+1)}
% for j = 1,3,5,... are blind spots of aliasing, since their values on the
% 1st kind Chebyshev grid are not coincident with the values of any T_m for
% 0 <= m <= n.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = size(coeffs, 1);
n2 = size(coeffs, 2);

% Pad with zeros:
if ( m > n )
    coeffs = [ zeros(m-n, size(coeffs, 2)); coeffs ];
    return
end

% It's more natural to work with the coefficients in the other order:
coeffs = coeffs(end:-1:1,:);

% Alias coefficients (see discussion above):
if ( m == 1 )
    % Reduce to a single point:
    e = ones(1, ceil(n/2)); 
    e(2:2:end) = -1;
    coeffs = e*coeffs(1:2:end,:);
elseif ( m > n/2 )
    % If m > n/2, only single coefficients are aliased, and we can vectorise.
    j = ((m + 1):n).';
    k = abs(mod(j + m - 2, 2*m) - m + 1) + 1;
    p = floor((j-1+m)/(2*m));
    t = (-1).^p;
    coeffs(k,:) = coeffs(k,:) + repmat(t,1,n2).*coeffs(j,:);
else
    % Otherwise we must do everything in a tight loop. (Which is slower!)
    for j = (m + 1):n
        k = abs(mod(j + m - 2, 2*m) - m + 1) + 1;
        p = floor((j-1+m)/(2*m));
        t = (-1)^p;
        coeffs(k,:) = coeffs(k,:) + t*coeffs(j,:);
    end
end

% Flip the coefficients back again:
coeffs = coeffs(m:-1:1,:);

end
