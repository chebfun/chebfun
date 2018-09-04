function pass = test_imag( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(r,lam,th)exp(1i*lam));
V = imag(ballfunv(f,f,f));
g = ballfun(@(r,lam,th)sin(lam));
exact = ballfunv(g,g,g);
pass(1) = norm(V-exact)<tol;

% Example 2
f = ballfun(@(r,lam,th)exp(1i*th.*lam));
V = imag(ballfunv(f,f,f));
g = ballfun(@(r,lam,th)sin(lam.*th));
exact = ballfunv(g,g,g);
pass(2) = norm(V-exact)<tol;

% Example 3
f1 = ballfun(@(r,lam,th)exp(1i*lam));
f2 = ballfun(@(r,lam,th)exp(1i*th));
f3 = ballfun(@(r,lam,th)exp(1i*th.*lam));
V = imag(ballfunv(f1,f2,f3));
g1 = ballfun(@(r,lam,th)sin(lam));
g2 = ballfun(@(r,lam,th)sin(th));
g3 = ballfun(@(r,lam,th)sin(lam.*th));
exact = ballfunv(g1,g2,g3);
pass(3) = norm(V-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
