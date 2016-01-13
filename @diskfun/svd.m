function varargout = svd( f )
%SVD    Singular value decomposition of a DISKFUN.
%
%   SVD(F) returns the singular values of F. The number of singular values
%   returned is equal to the rank of the DISKFUN.
%
%   S = SVD(F) returns the singular values of F. S is a vector of singular
%   values in decreasing order.
%
%   [U, S, V] = SVD(F) returns the SVD of F. V is a quasi-matrix of
%   orthogonal CHEBFUN objects, U is a quasimatrix of CHEBFUN objects that
%   are orthogonal with respect to the r weight on [0,1] (derived from the
%   measure on the disk) and S is a diagonal matrix with the singular
%   values on the diagonal.
%
%   The length and rank of a DISKFUN are slightly different quantities.
%   LENGTH(F) is the number of pivots used by the constructor, and
%   RANK(F) is the number of significant singular values of F. The relation
%   RANK(F) <= LENGTH(F) should always hold.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) )
    varargout = { [ ] };
    return
end

% Get CDR decomposition of f:
[C, D, R] = cdr( f );

% Generate normalized Bessel functions J_0( l(k) r ), where l(k) 
% is the kth root of J_0.  These set of functions satisfy the following:
%
%   int_0^\pi J_0(l(j)r).*J_0(l(k)r).*r dr = delta_{jk}.

E = C; 
rts = GLRbesselroots(0, size(C,2));
r = chebfun(@(r) r );
for kk = 0:size(C,2)-1
    nrm = sqrt(2)./besselj(1, rts(kk+1));
    E(:,kk+1) = besselj( 0, rts(kk+1)*r )*nrm;
end

% The weighted inner-product that will be used in the theta-variable:
weight = chebfun( @(r) r );
myInnerProduct = @(u, v) sum( weight.*conj(u).*v, [0, 1] );

% Do QR in both variables, one with the weighted inner-product and one with
% the standard L2 inner-product.
[QwC, RwC] = abstractQR(C, E, myInnerProduct);
[QwR, RwR] = qr( R );

% Use the QR factorizations of the columns and rows to make up the SVD of
% the DISKFUN object.  Since
%
%        C * D * R = QwC * ( RwC * D * RwR' ) * QwR'
%
% we compute the SVD of ( RwC * D * RwR' ).
[U, S, V] = svd( RwC * D * RwR.' );
U = QwC * U;
V = QwR * V;

% Output just like the svd of a matrix.
if ( nargout > 1 )
    varargout = { U, S, V };
else
    varargout = { diag( S ) };
end

end


%%%%%%%%%%%%%%%%% BESSELROOTS %%%%%%%%%%%%%%%%%%%
% Calculate the bessel roots using the GLR algorithm.  We could use any
% other algorithm here: 
function [roots, ders] = GLRbesselroots(n,numroots)
roots = zeros(numroots,1); ders = zeros(numroots,1);

newton_roots = 3;
if n == 0
    xs = 2.404825557695773;
elseif n > 0
    % See Hethcote 1970
    xs = n + 1.8557*n^(1/3);
else
    n1 = n + 1;
    % See Piessens 1984
    xs = 2*sqrt(n+1)*(1 + n1/4 - 7*n1^2/96  + 49*n1^3/1152 - 8363*n1/276480);
    newton_roots = min(max(2*ceil(abs(log10(n1))),3),numroots);
end

if n ~=0
    % The first root (Newton iteration)
    [roots(1) ders(1)] = alg3_Bes(n,xs);
else
    roots(1) = xs;
end
if numroots == 1, return, end

% The second root (Newton iteration)
[roots(2) ders(2)] = alg3_Bes(n,roots(1)+.9*pi);
if numroots == 2, return, end

% Some more roots (Newton)
for k = 3:newton_roots
    [roots(k) ders(k)] = alg3_Bes(n,roots(k-1)+.99*pi);
end

% The larger roots (use GLR).
[roots ders] = alg1_Bes(roots,ders,n,numroots,newton_roots);

end
% -------------------------------------------------------------------------

function [roots ders] = alg1_Bes(roots,ders,n,numroots,start)
m = 50; kk = 0:m; fact = 1./factorial(kk);
do_rk = 1;

for j = start:numroots-1
    x = roots(j);
    if do_rk
        h = rk2_Bes(pi/2,-pi/2,roots(j),n) - x;
    else
        h = pi;
    end
    if do_rk && abs(h-pi)<1e-3, do_rk = 0; end
    
    u = zeros(1,m); u(1) = 0;   u(2) = ders(j);
    r = x.^2-n^2;   p = x.^2;
    u(3) = (-x*u(2)-r*u(1))/p;
    u(4) = (-3*x*u(3)-(1+r)*u(2)-2*x*u(1))/p;
    for k = 2:m-2
        u(k+3) = (-x*(2*k+1)*u(k+2)-(k*(k-1)+k*1+r)*u(k+1)-2*k*x*u(k)-k*(k-1)*u(k-1))/p;
    end
    G1 = (u.*fact); G2 = (u(2:m+1).*fact(1:m));

    for l = 1:5
        hh = [1;cumprod(h+zeros(m,1))]; 
        dh = (G1*hh)./(G2*hh(1:m));
        h = h - dh;
    end
    
    roots(j+1) = roots(j) + h;
    ders(j+1) = G2*[1;cumprod(h+zeros(m-1,1))];
end
end

% -------------------------------------------------------------------------

function [x1 d1] = alg3_Bes(n,xs)
u = besselj(n,xs,0);
up =  besselj(n-1,xs,0)-n/xs*u;
theta = atan(xs./sqrt(xs.^2-n^2)*up/u);
x1 = real(rk2_Bes(theta,-pi/2,xs,n));

for k = 1:20
    u = besselj(n,x1,0);
    du = besselj(n-1,x1,0)-n/x1*u;
    dx = u/du;
    x1 = x1 - dx;
    if abs(dx) < eps, break; end
end
d1 = besselj(n-1,x1,0)-n/x1*u;
end

% -------------------------------------------------------------------------

function x = rk2_Bes(t,tn,x,n)
m = 10; h = (tn-t)/m;
for j = 1:m
    k1 = -h./(sqrt(1-(n/x)^2) + .25*(x^2+n^2)/(x^2-n^2)/x.*sin(2*t));
    t = t+h;  x = x+k1;
    k2 = -h./(sqrt(1-(n/x)^2) + .25*(x^2+n^2)/(x^2-n^2)/x.*sin(2*t));
    x = x+.5*(k2-k1);
end
end