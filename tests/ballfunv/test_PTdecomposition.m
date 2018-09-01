function pass = test_PTdecomposition( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1 :
S = [51,51,51];
p = ballfun(@(r,lam,th)cos(r.^2.*sin(th).^2.*cos(lam).*sin(lam)),S);
t = ballfun(@(r,lam,th)sin(r.^2.*sin(th).*cos(th).*sin(lam)),S);
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(1) = norm(diff(p,[0,1,0])-diff(p2,[0,1,0]))<tol;
pass(2) = norm(diff(p,[0,0,1])-diff(p2,[0,0,1]))<tol;
pass(3) = norm(diff(t,[0,1,0])-diff(t2,[0,1,0]))<tol;
pass(4) = norm(diff(t,[0,0,1])-diff(t2,[0,0,1]))<tol;

% Example 2 :
S = [53,53,53];
p = ballfun(@(r,lam,th)cos(r.^2.*sin(th).^2.*cos(lam).*sin(lam)),S);
t = ballfun(@(r,lam,th)sin(r.^2.*sin(th).*cos(th).*sin(lam)),S);
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(5) = norm(diff(p,[0,1,0])-diff(p2,[0,1,0]))<tol;
pass(6) = norm(diff(p,[0,0,1])-diff(p2,[0,0,1]))<tol;
pass(7) = norm(diff(t,[0,1,0])-diff(t2,[0,1,0]))<tol;
pass(8) = norm(diff(t,[0,0,1])-diff(t2,[0,0,1]))<tol;

% Example 3 :
S = [41,41,41];
p = randnfunball(5,S);
t = randnfunball(5,S);
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(9) = norm(diff(p,[0,1,0])-diff(p2,[0,1,0]))<tol;
pass(10) = norm(diff(p,[0,0,1])-diff(p2,[0,0,1]))<tol;
pass(11) = norm(diff(t,[0,1,0])-diff(t2,[0,1,0]))<tol;
pass(12) = norm(diff(t,[0,0,1])-diff(t2,[0,0,1]))<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
