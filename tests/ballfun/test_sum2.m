function pass = test_sum2( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(r,lam,th)1);
ff = sum2(f,2,3);
g = chebfun(@(r)4*pi*r.^2);
pass(1) =norm( ff - g ) <tol;

% Example 2
f = ballfun(@(r,lam,th)cos(lam));
ff = sum2(f,3,2);
g = chebfun(@(r)0);
pass(2) =norm( ff - g ) <tol;

% Example 3
f = ballfun(@(r,lam,th)r.^2.*sin(lam));
ff = sum2(f,2,3);
g = chebfun(@(r)0);
pass(3) =norm( ff - g ) <tol;

% Example 4
f = ballfun(@(r,lam,th)r.*sin(th));
ff = sum2(f,1,2);
g = chebfun(@(th)pi*sin(th),[-pi,pi],'trig');
pass(4) = norm( ff - g ) <tol;

% Example 5
f = ballfun(@(r,lam,th)r.^2.*sin(th).*cos(lam));
ff = sum2(f,1,3);
g = chebfun(@(lam)2*cos(lam)/3,[-pi,pi],'trig');
pass(5) = norm( ff - g ) < tol;


if (nargout > 0)
    pass = all(pass(:));
end

end
