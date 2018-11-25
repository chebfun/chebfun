function pass = test_feval( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

%% Evaluate at spherical coordinates

% Example 1
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th),'spherical');
F = feval(f,[0.5,0.7], [0,0], [pi/2,pi/2], 'spherical');
exact = [0.5, 0.7];
pass(1) = norm(F(:)-exact(:)) < tol;

% Example 2
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th),'spherical');
F = feval(f,[1,1], [pi/4,pi/3], [pi/2,pi/2], 'spherical');
exact = [cos(pi/4); cos(pi/3)];
pass(2) = norm(F(:)-exact(:)) < tol;

% Example 3
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th),'spherical');
F = feval(f,[1,1], [0,0], [pi/5, pi/7], 'spherical');
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

%% Evaluate at cartesian coordinates

% Example 5
f = ballfun(@(x,y,z)x);
F = f(1,0,0);
pass(5) = abs(1-F) < tol;

% Example 6
f = ballfun(@(x,y,z)y);
F = f(0,0.5,0);
pass(6) = abs(0.5-F) < tol;

% Example 7
f = ballfun(@(x,y,z)z);
F = f(0,0,-0.3);
pass(7) = abs(-0.3-F) < tol;

%% EXTRACT_SPHEREFUN EXAMPLES: 

% Example 8
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th),'spherical');
g = f(1,:,:,'spherical');
h = spherefun(@(lam,th)cos(lam).*sin(th));
pass(8) = norm( g - h ) < tol;

% Example 9
f = ballfun(@(r,lam,th)r.*cos(th),'spherical');
g = f(.5,:,:,'spherical');
h = spherefun(@(lam,th)0.5.*cos(th));
pass(9) = norm( g - h ) < tol;

%% COMPLEX EVALUATION
f = ballfun(@(x,y,z) cos(x)*1i);
pass(10) = abs(f(0,0,0)-1i);

if (nargout > 0)
    pass = all(pass(:));
end
end
