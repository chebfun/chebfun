% Test file for @chebfun/chebpoly.m.

function pass = test_chebpoly(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

warnState = warning('off', 'CHEBFUN:CHEBFUN:chebpoly:deprecated');

rev = @(A) flipud(A).';

% Test on a smooth chebfun.
f = chebfun(@(x) cos(x));
n = 5;
c = chebpoly(f, n);
c_exact = [besselj(0, 1) ;
           0 ;
           -2*besselj(2, 1) ;
           0 ;
           2*(-23*besselj(0, 1) + 40*besselj(1, 1))];
err = c - rev(c_exact);
pass(1) = norm(err, inf) < 1e2*vscale(f)*eps;


g = [f f];
c = chebpoly(g, n);
c_exact = [c_exact c_exact];
err = c - rev(c_exact);
pass(2) = norm(err(:), inf) < 1e2*vscale(f)*eps;


c = chebpoly(f, n, 'kind', 2);
c_exact = [2*besselj(1, 1) ;
           0 ;
           -6*besselj(3, 1) ;
           0 ;
           2*(-235*besselj(1, 1) + 900*besselj(2, 1))];
err = c - rev(c_exact);
pass(3) = norm(err, inf) < 1e3*vscale(f)*eps;


c = chebpoly(g, n, 'kind', 2);
c_exact = [c_exact c_exact];
err = c - rev(c_exact);
pass(4) = norm(err(:), inf) < 1e3*vscale(f)*eps;


% Test on a piecewise-smooth chebfun
f = chebfun(@(x) abs(x), [-1 0 1]);
n = 7;
c = chebpoly(f, n);
c_exact = [2/pi ; 0 ; 4/(3*pi) ; 0 ; -4/(15*pi) ; 0 ; 4/(35*pi)];
err = c - rev(c_exact);
pass(5) = norm(err, inf) < 10*vscale(f)*eps;

g = [f f];
c = chebpoly(g, n);
c_exact = [c_exact c_exact];
err = c - rev(c_exact);
pass(6) = norm(err(:), inf) < 10*vscale(f)*eps;

c = chebpoly(f, n, 'kind', 2);
c_exact = [4/(3*pi) ; 0 ; 4/(5*pi) ; 0 ; -4/(21*pi) ; 0 ; 4/(45*pi)];
err = c - rev(c_exact);
pass(7) = norm(err, inf) < 1e2*vscale(f)*eps;

c = chebpoly(g, n, 'kind', 2);
c_exact = [c_exact c_exact];
err = c - rev(c_exact);
pass(8) = norm(err(:), inf) < 1e2*vscale(f)*eps;

warning(warnState);

end
