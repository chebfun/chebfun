function g = spherefun(f, r)
% SPHEREFUN returns the evaluation of F at the given radius 0<=R<=1
% G = SPHEREFUN(f, r) is the SPHEREFUN function g(lambda, theta) = f(r, :, :).
%
% % See also BALLFUN/DISKFUN.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( f )
   g = spherefun();
   return
end

F = f.coeffs;
[m,n,p] = size(f);

% r = 1 by default
if nargin == 1
    r = 1;
end

if m == 1
    G = reshape(F(1,:,:),n,p);
else
    % Chebyshev functions evaluated at r
    T = zeros(1,m);
    T(1) = 1; T(2) = r;
    for i = 3:m
        T(i) = 2*r*T(i-1)-T(i-2);
    end

    % Build the array of coefficient of the spherefun function
    G = zeros(n,p);
    for i = 1:p
        G(:,i) = T*F(:,:,i);
    end
end
% Build the spherefun function; coeffs2spherefun takes the theta*lambda matrix
% of coefficients
g = spherefun.coeffs2spherefun(G.');
end