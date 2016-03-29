function pass = test_laplacian( )

tol = 1e4*chebfunpref().cheb2Prefs.chebfun2eps;

k = 1;
for ell = [1 2 4 5 7 8 9]
    for m = -ell:ell
        f = spherefun.sphharm(ell, m);
        lap = laplacian(f);
        pass(k, 1) = numel(lap.pivotValues) == numel(f.pivotValues);
        err(k) = SampleError(-ell*(ell+1)*f, lap)/(ell*(ell+1));
        pass(k, 2) = SampleError(-ell*(ell+1)*f, lap) < ell*(ell+1)*tol;
        k = k+1;
    end
end
pass = pass(:)';

% Gaussian
lam0 = pi/sqrt(2);
th0 = (sqrt(5) - 2)/2;
sig2 = 1;
r2 = @(lam,th) 2*(1 - (sin(th)*sin(th0)).*cos(lam-lam0)-cos(th)*cos(th0));
f = @(lam,th) exp(-sig2 * r2(lam, th));
exact = @(lam,th) sig2*exp(-sig2*r2(lam,th)).*(-4 + r2(lam,th).*...
    (2 - sig2*(-4 + r2(lam, th))));
lap = laplacian(spherefun(f));
pass(end+1) = SampleError(exact, lap) < 10*tol;

end

function sample_error = SampleError(h, g)
m = 6; 
n = m;
[x, y] = getPoints(m, n);
[L2, T2] = meshgrid(x, y);
F = feval(h, L2, T2);
% hn = pi/length(g.cols);  % Have to adjust for the shift in y points. Ugly!
% approx = feval(g.cols,(y-hn)/pi-.5) * g.blockDiag * feval(g.rows,x/pi)';
% sample_error = norm( F - approx , inf );
% [C,D,R] = cdr(g);
% approx = feval(g.cols,y/pi) * g.blockDiag * feval(g.rows,x/pi)';
approx = fevalm(g, x, y);
sample_error = norm(F(:) - approx(:), inf);
end

function [x, y] = getPoints(m, n)

% x = linspace(-pi, pi, 2*n+1)';  x( end ) = [ ];
% GBW: You can't just remove the pole and keep everything else equally
% spaced between -pi/2 and 3*pi/2.  The issue is that you can't keep
% both point -pi/ and 3*pi/2.
% y = linspace(-pi/2, 3*pi/2, 2*m+1)'; y( m+1 ) = [ ];

x = trigpts(2*n, [-pi pi]);
% GBW: I believe we have to sample at equally spaced points shifted by h/2
% to not sample the poles and keep an even total number of points.
% y = trigpts(2*m,[-pi/2 3*pi/2]);
% y = trigpts(m,[0 pi]);
y = linspace(0, pi, m).';

% y = y+0.5*pi/m; % Shift y by h/2 to avoid the poles
end