function pass = test_vals2coeffs( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
X = ones(21,22,23);
V = ballfun.vals2coeffs(X);
V2 = zeros(21,22,23);
V2(1,12,12) = 1;
pass(1) = norm(V(:)-V2(:)) < tol;

% Example 2
S = [100,150,200];
m = S(1); n = S(2); p = S(3);
X = ones(S);
V = ballfun.vals2coeffs(X);
V2 = zeros(S);
V2(1,floor(n/2)+1,floor(p/2)+1) = 1;
pass(2) = norm(V(:)-V2(:)) < tol;

% Example 3
f = @(lam)exp(1i*lam);
lam = pi*trigpts(3);
vals = f(lam).';
cfs = [0,0,1];
pass(3) = norm(ballfun.vals2coeffs(vals)-cfs) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
