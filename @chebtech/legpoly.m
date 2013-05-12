function b = legpoly(f)
%LEGPOLY   Compute Legendre series coefficients of a CHEBTECH object.
%   B = LEGPOLY(F) returns the Legendre series coefficients of CHEBTECH F, so
%   that F = B(N+1)*P_N + ... + B(1)*P_0, where P_k is the kth Legendre
%   polynomial.
%
%   If F is an array-valued CHEBETCH, then a matrix of coefficients is returned
%   so that F(:,k) = B(N+1,k)*P_N + ... + B(1,k)*P_0.
%
% See also CHEBPOLY.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is currently implemented by 'inverting' the Vandermonde-type matrix. In
% fact the matrix can be inverted explicitly by nothing that D*P'W*P = I, where
% P is the Legendre Vandemonde matrix evaluated on a Chebyshev grid of size
% 2*n+1, W is a diagonal matrix of the corresponding Clenshaw-Curtis weights,
% and the matrix D is the diagonal scaling ((0:n)'+1/2);
%
% In theory, the O(N^2) coefficient-based bethod due to Piessen's [1] should be
% faster (i.e., have a smaller implied constant), but the tight FOR loop is
% prohibitively slow in MATLAB.
%
% [1] R. Piessens, "Computation of Legendre series coefficients", Communications
% of the ACM, Volume 17 Issue 1, Jan. 1974, Page 25.
%
% [TODO]: Perhaps it's better to evaluate on a Legendre grid? Then no doubling.
% [TODO]: Develop a fast algorithm for this transformation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(f, 1) - 1;

% Create a Chebyshev grid of double the length:
x = f.chebpts(2*n + 1);
w = f.quadwts(2*n + 1);

% Evaluate the function on this grid (using FFTs):
v = get(prolong(f, 2*n + 1), 'values');

% Initialise the Vandermonde-like recurrence:
V = zeros(length(x), n);
Pm2 = 1; 
Pm1 = x;
V(:,1:2) = [(1 + 0*x), x];
for k = 1:n-1
    % Legendre polynomial recurrence:
    P = ((2*k + 1)*Pm1.*x - k*Pm2) / (k + 1);
    Pm2 = Pm1; 
    Pm1 = P; 
    % Store in a matrix:
    V(:,2+k) = P;
end

% Use weighted discrete orthogonality to invert the system (see doc above):
% b = ((0:n)'+1/2).*(V'*(w'.*v)); % b = (V\v) (works for scalar v)

% For array-valued CHEBTECH objects we need to use BSXFUN():
b = bsxfun(@times, (0:n)' + 0.5, V'*bsxfun(@times, w', v)); % b = (V\v)

end

%% Piessens' method [1] is much slower in MATLAB because of the tight FOR loops.
% % Obtain and flip the coefficiets:
% c = f.coeffs;
% c = c(end:-1:1,:);
% n = size(c, 1) - 1;
% m = size(c, 2);
% 
% % Initialise the output vector:
% b = zeros(n+1, m);
% 
% % The first entry is different:
% % Make the diagonal entries of the conversion matrix:
% d = [1 ; 2/3 ; zeros(n-2, 1)];
% for k = 2:n
%     d(k+1,1) = k/(k+.5)*d(k,1);
% end
% evenIdx = 0:2:n;
% b(1,:) = 2*(evenIdx./(evenIdx.^2 - 1) - 1./(evenIdx - 1)) * c(1:2:end,:);
% 
% % Loop over remaining entries:
% for j = 1:n
%     % Initialise e to the diagonal entry:
%     e = d(j+1);
%     % Loops along entries in the jth row:
%     jj1 = j*(j+1);
%     for k = j:2:n
%         % Add the kth component:
%         b(j+1,:) = b(j+1,:) + e*c(k+1,:);
%         % Update e via recurrence:
%         e = ((k-1)*k - jj1)*(k+2) / (((k+3)*(k+2) - jj1)*k) * e;
%     end
% end
% b = ((0:n)'+1/2).*b;
%
% % Flip the entries for output:
% % b = b(end:-1:1,:)

