function pass = test_sum3( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(r,lam,th)1,[20,20,20]);
I = sum3(f);
exact = 4*pi/3;
pass(1) = abs(I-exact)<tol;

% Example 2
f = ballfun(@(r,lam,th)cos(lam),[20,20,20]);
I = sum3(f);
exact = 0;
pass(2) = abs(I-exact)<tol;

% Example 3
f = ballfun(@(r,lam,th)r.^2.*sin(lam),[20,20,20]);
I = sum3(f);
exact = 0;
pass(3) = abs(I-exact)<tol;

% Example 4
f = ballfun(@(r,lam,th)cos(th),[20,20,20]);
I = sum3(f);
exact = 0;
pass(4) = abs(I-exact)<tol;

% Example 5
f = ballfun(@(r,lam,th)sin(th),[20,20,20]);
I = sum3(f);
exact = pi^2/3;
pass(5) = abs(I-exact)<tol;

% Example 6
f = ballfun(@(r,lam,th)r.*sin(th).*(1+sin(lam)),[20,20,20]);
I = sum3(f);
exact = pi^2/4;
pass(6) = abs(I-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end

end
