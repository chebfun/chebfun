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
[~,n,p] = size(f);

% r = 1 by default
if nargin == 1
    r = 1;
end

G = zeros(1, n, p);
% Evaluate f at the points r
for i = 1:p
    G(:, :, i) = clenshaw_vec( r, F(:, :, i) );
end
G = reshape(G, n, p);

% Build the spherefun function; coeffs2spherefun takes the theta*lambda matrix
% of coefficients
g = spherefun.coeffs2spherefun(G.');
end

% Pulled from chebtech/clenshaw: 
function y = clenshaw_vec(x, c)
% Clenshaw scheme for array-valued functions.
x = repmat(x(:), 1, size(c, 2));
bk1 = zeros(size(x, 1), size(c, 2)); 
bk2 = bk1;
e = ones(size(x, 1), 1);
x = 2*x;
n = size(c, 1)-1;
for k = (n+1):-2:3
    bk2 = e*c(k,:) + x.*bk1 - bk2;
    bk1 = e*c(k-1,:) + x.*bk2 - bk1;
end
if ( mod(n, 2) )
    [bk1, bk2] = deal(e*c(2,:) + x.*bk1 - bk2, bk1);
end
y = e*c(1,:) + .5*x.*bk1 - bk2;
end