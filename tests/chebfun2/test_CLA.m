function pass = test_CLA( pref )
% Try out some basic continuous linear algebra.

if ( nargin < 1 )
    pref = chebfunpref; 
end
tol = 200*pref.cheb2Prefs.chebfun2eps;

seedRNG(0)
gam = 10;
% Add some Gaussians:
f = @(x, y) 0;
for n = 1:20
    x0 = 2*rand-1;
    y0 = 2*rand-1;
    df = @(x,y) exp(-gam*((x-x0).^2 + (y-y0).^2));
    f = @(x, y) f(x, y) + df(x, y);
end

% Make a CHEBFUN2:
F = chebfun2(f);

% Check QR of a chebfun2

% Test QR of the columns:
[Q, R] = qr(F.cols);
pass(1) = ( norm(Q*R - F.cols) < tol );

% Test QR of the rows:
[Q, R] = qr(F.rows);
pass(2) = ( norm(Q*R - F.rows) < tol );

% Check CDR
[C, D, R] = cdr(F);
% Make D positive:
for k = 1:size(D)
    if ( D(k, k) < 0 )
        D(k, k) = -D(k, k);
        C(:, k) = -C(:,k);
    end
end

G1 = (C*D)*R.';
pass(3) = ( norm(G1 - F) < tol );
pass(4) = ( norm(F - G1) < tol );

G2 = C*(D*R.');
pass(5) = ( norm(G2 - F) < tol );
pass(6) = ( norm(F - G2) < tol );

G3 = (C*sqrt(D))*(sqrt(D)*R.');
pass(7) = ( norm(G3 - F) < tol );
pass(8) = ( norm(F - G3) < tol );

end
