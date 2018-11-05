function pass = test_PTdecomposition( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e5*pref.techPrefs.chebfuneps;

% Example 1 :
p = ballfun(@(r,lam,th)cos(r.^2.*sin(th).^2.*cos(lam).*sin(lam)));
t = ballfun(@(r,lam,th)sin(r.^2.*sin(th).*cos(th).*sin(lam)));
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(1) = norm(diff(p,2)-diff(p2,2))<tol;
pass(2) = norm(diff(p,3)-diff(p2,3))<tol;
pass(3) = norm(diff(t,2)-diff(t2,2))<tol;
pass(4) = norm(diff(t,3)-diff(t2,3))<tol;


% Example 2 :
p = ballfun(@(x,y,z)x.^2+y.*z,'cart');
t = ballfun(@(x,y,z)x.*y.*z,'cart');
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(5) = norm(diff(p,2)-diff(p2,2))<tol;
pass(6) = norm(diff(p,3)-diff(p2,3))<tol;
pass(7) = norm(diff(t,2)-diff(t2,2))<tol;
pass(8) = norm(diff(t,3)-diff(t2,3))<tol;

% Example 3:
S = [50,50,50];
p = ballfun(@(x,y,z)x.^2+y.*z,'cart',S);
t = ballfun(@(x,y,z)x.*y.*z,'cart',S);
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(9) = norm(diff(p,2)-diff(p2,2))<tol;
pass(10) = norm(diff(p,3)-diff(p2,3))<tol;
pass(11) = norm(diff(t,2)-diff(t2,2))<tol;
pass(12) = norm(diff(t,3)-diff(t2,3))<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end