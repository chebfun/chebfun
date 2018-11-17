function pass = test_feval( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th),'polar');
F = feval(f,[0.5,0.7], [0,0], [pi/2,pi/2]);
exact = [0.5, 0.7];
pass(1) = norm(F(:)-exact(:)) < tol;

% Example 2
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th),'polar');
F = feval(f,[1,1], [pi/4,pi/3], [pi/2,pi/2]);
exact = [cos(pi/4); cos(pi/3)];
pass(2) = norm(F(:)-exact(:)) < tol;

% Example 3
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th),'polar');
F = feval(f,[1,1], [0,0], [pi/5, pi/7]);
exact = zeros(1,1,2);
exact(1,1,1) = sin(pi/5);
exact(1,1,2) = sin(pi/7);
pass(3) = norm(F(:)-exact(:)) < tol;

% Example 4
S = [22,23,25];
f = ballfun(@(x,y,z)x.*y);
r = chebpts(S(1));
lam = pi*trigpts(S(2));
th = pi*trigpts(S(3));
exact = ballfun.coeffs2vals(coeffs3(f,S(1),S(2),S(3)));
F = fevalm(f,r,lam,th);
pass(4) = norm(F(:)-exact(:)) < tol;

%% EXTRACT_SPHEREFUN EXAMPLES: 
% Example 1
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th),'polar');
g = f(1,:,:);
h = spherefun(@(lam,th)cos(lam).*sin(th));
pass(5) = norm( g - h ) < tol;

% Example 2
f = ballfun(@(r,lam,th)r.*cos(th),'polar');
g = f(.5,:,:);
h = spherefun(@(lam,th)0.5.*cos(th));
pass(6) = norm( g - h ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
