function pass = test_PTdecomposition( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e5*pref.techPrefs.chebfuneps;

% Example 1 :
p = ballfun(@(r,lam,th)cos(r.^2.*sin(th).^2.*cos(lam).*sin(lam)), 'spherical');
t = ballfun(@(r,lam,th)sin(r.^2.*sin(th).*cos(th).*sin(lam)), 'spherical');
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(1) = norm(diff(p,2,'spherical')-diff(p2,2,'spherical'))<tol;
pass(2) = norm(diff(p,3,'spherical')-diff(p2,3,'spherical'))<tol;
pass(3) = norm(diff(t,2,'spherical')-diff(t2,2,'spherical'))<tol;
pass(4) = norm(diff(t,3,'spherical')-diff(t2,3,'spherical'))<tol;


% Example 2 :
p = ballfun(@(x,y,z)x.^2+y.*z);
t = ballfun(@(x,y,z)x.*y.*z);
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(5) = norm(diff(p,2,'spherical')-diff(p2,2,'spherical'))<tol;
pass(6) = norm(diff(p,3,'spherical')-diff(p2,3,'spherical'))<tol;
pass(7) = norm(diff(t,2,'spherical')-diff(t2,2,'spherical'))<tol;
pass(8) = norm(diff(t,3,'spherical')-diff(t2,3,'spherical'))<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
