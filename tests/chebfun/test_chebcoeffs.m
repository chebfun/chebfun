function pass = test_chebcoeffsy(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% Test on a smooth chebfun.
f = chebfun(@(x) cos(x));
n = 5;
c = chebcoeffs(f, n);
c_exact = [besselj(0, 1) ;
           0 ;
           -2*besselj(2, 1) ;
           0 ;
           2*(-23*besselj(0, 1) + 40*besselj(1, 1))];
err = c - c_exact;
pass(1) = norm(err, inf) < 10*vscale(f).*epslevel(f);

% Test on a smooth array-valued chebfun.
f = [f f];
c = chebcoeffs(f, n);
c_exact = [c_exact c_exact];
err = c - c_exact;
pass(2) = norm(err(:), inf) < 10*vscale(f).*epslevel(f);

% Test on piecewise-smooth chebfun
f = chebfun(@(x) abs(x), [-1 0 1]);
n = 7;
c = chebcoeffs(f, n);
c_exact = [2/pi ; 0 ; 4/(3*pi) ; 0 ; -4/(15*pi) ; 0 ; 4/(35*pi)];
err = c - c_exact;
pass(3) = norm(err, inf) < 10*vscale(f).*epslevel(f);

% Test on a piecewise-smooth array-valued chebfun.
f = [f f];
c = chebcoeffs(f, n);
c_exact = [c_exact c_exact];
err = c - c_exact;
pass(4) = norm(err(:), inf) < 10*vscale(f).*epslevel(f);

end
