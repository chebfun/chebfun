function [f, R, E] = qr(f, flag)
%QR     QR factorisation of a multivalued FUNCHEB2.
%   [Q, R] = QR(F) returns a QR factorisation of F such that F = Q*R, where the
%   FUNCHEB2 Q is orthogonal (wrt the continuous L^2 norm on [-1,1]) and of the
%   same size as F and R is an m x m upper-triangular matrix when F has m
%   columns.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Grab the size of f:
m = size(f, 2);

% If F has only one column we simply scale it.
if ( m == 1 )
    R = sqrt(sum(f.*f));
    f = f./R;
    E = 1;
    return
end

% Simplify so that we don't do any extra work: (QR is O(m*n^2)? :/ )
f = simplify(f);
n = size(f, 1);

% Project the values onto a Legendre grid: (where integrals of polynomials
% p_n*q_n will be computed exactly and on an n-point grid)
xc = funcheb2.chebpts(n);
[xl, wl, vl] = legpts(n);
P = barymat(xl, xc);
W = spdiags(sqrt(wl.'), 0, n, n);

% Compute the weighted QR factorisation:
if ( nargout == 3 )
    [Q, R, E] = qr(W * P * f.values, 0);
    % For consistency with MATLAB's QR behavior:
    if ( nargin == 1 || ~( strcmpi(flag, 'vector') || flag == 0 ) )
        % Return E in matrix form:
        I = eye(m);
        E = I(:,E);
    end
else
    [Q, R] = qr(W * P * f.values, 0);
end

% Revert to the Chebyshev grid (and remove the weight and enforce diag(R) >= 0).
Winv = diag(1./sqrt(wl));   % Undo the weighting used for QR.
Pinv = barymat(xc, xl, vl); % Revert to Chebyshev grid (from Legendre).
S = spdiags(sign(diag(R)), 0, m, m);    % Enforce diag(R) >= 0.
Q = Pinv*Winv*Q*S;          % Fix Q.
R = S*R;                    % Fix R.

f.values = Q;                           % Adjust values of f.
f.coeffs = funcheb2.chebpoly(Q);        % Compute new coeffs.

end




