function b = legpoly(f)
%LEGPOLY    Compute Legendre series coefficients of a CHEBTECH object.
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
% This is currently implemented using an O(N^2) algorithm due to Piessen's: R.
% Piessens, "Computation of Legendre series coefficients", Communications of the
% ACM, Volume 17 Issue 1, Jan. 1974, Page 25.
%
% [TODO]: Develop a fast algorithm for this transformation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain and flip the coefficiets:
c = f.coeffs;
c = c(end:-1:1,:);
n = size(c, 1) - 1;
m = size(c, 2);

% Initialise the output vector:
b = zeros(n+1, m);

% The first entry is different:
% Make the diagonal entries of the conversion matrix:
d = [1 ; 2/3 ; zeros(n-2, 1)];
for k = 2:n
    d(k+1,1) = k/(k+.5)*d(k,1);
end
evenIdx = 0:2:n;
b(1,:) = 2*(evenIdx./(evenIdx.^2 - 1) - 1./(evenIdx - 1)) * c(1:2:end,:);

% Loop over remaining entries:
for j = 1:n
    % Initialise e to the diagonal entry:
    e = d(j+1);
    % Loops along entries in the jth row:
    for k = j:2:n
        % Add the kth component:
        b(j+1,:) = b(j+1,:) + e*c(k+1,:);
        % Update e via recurrence:
        e = ((k-1)*k - j*(j+1))*(k+2) / (((k+3)*(k+2) - j*(j+1))*k) * e;
    end
end

% Flip the entries for output:
b = b(end:-1:1,:);

end