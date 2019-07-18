function pass = test_mixed_tech( pref )
% This tests the chebfun2 constructor for mixed technology.

if ( nargin < 1 )
    pref = chebfunpref;
end
tol = 1e2 * pref.cheb2Prefs.chebfun2eps;

pass = [];

%% 'trig' in one or both dimensions

% Construct from function handle
u = @(x,y) sin(2*pi*x).*cos(2*pi*y);

f1 = chebfun2(u, {'trig', []});
f2 = chebfun2(u, {[], 'periodic'});
f3 = chebfun2(u, {'periodic', 'trig'});
f4 = chebfun2(u, 'trig');
f5 = chebfun2(u);

pass(end+1) = ( norm(f1 - f2) < tol );
pass(end+1) = ( norm(f2 - f3) < tol );
pass(end+1) = ( norm(f3 - f4) < tol );
pass(end+1) = ( norm(f4 - f5) < tol );

% Construct from data
m = 31;
n = 32;
u = @(x,y) sin(2*pi*x).*cos(2*pi*y);
[xc, yc] = chebpts2(n, m);
[xt, yt] = meshgrid(trigpts(n), trigpts(m));

f1 = chebfun2(u(xt,yc), {'trig', []});
f2 = chebfun2(u(xc,yt), {[], 'periodic'});
f3 = chebfun2(u(xt,yt), {'periodic', 'trig'});
f4 = chebfun2(u(xt,yt), 'trig');
f5 = chebfun2(u(xc,yc));

pass(end+1) = ( norm(f1 - f2) < tol );
pass(end+1) = ( norm(f2 - f3) < tol );
pass(end+1) = ( norm(f3 - f4) < tol );
pass(end+1) = ( norm(f4 - f5) < tol );

%% 'equi' in one or both dimensions
m = 31;
n = 32;
u = @(x,y) sin(x.*y);
[xc, yc] = chebpts2(n, m);
[xe, ye] = meshgrid(linspace(-1,1,n), linspace(-1,1,m));

f1 = chebfun2(u(xe,yc), {'equi', []});
f2 = chebfun2(u(xc,ye), {[], 'equi'});
f3 = chebfun2(u(xe,ye), {'equi', 'equi'});
f4 = chebfun2(u(xe,ye), 'equi');
f5 = chebfun2(u(xc,yc));

pass(end+1) = ( norm(f1 - f2) < tol );
pass(end+1) = ( norm(f2 - f3) < tol );
pass(end+1) = ( norm(f3 - f4) < tol );
pass(end+1) = ( norm(f4 - f5) < tol );

%% 'coeffs' in one or both dimensions
m = 31;
n = 32;
u = @(x,y) sin(x.*y);
[xc, yc] = chebpts2(n, m);

vals_vals     = u(xc,yc);
vals_coeffs   = chebtech2.vals2coeffs( vals_vals );
coeffs_vals   = chebtech2.vals2coeffs( vals_vals.' ).';
coeffs_coeffs = chebtech2.vals2coeffs( vals_coeffs.' ).';

f1 = chebfun2(u);
f2 = chebfun2(coeffs_vals, {'coeffs', []});
f3 = chebfun2(vals_coeffs, {[], 'coeffs'});
f4 = chebfun2(coeffs_coeffs, {'coeffs', 'coeffs'});
f5 = chebfun2(coeffs_coeffs, 'coeffs');

pass(end+1) = ( norm(f1 - f2) < tol );
pass(end+1) = ( norm(f2 - f3) < tol );
pass(end+1) = ( norm(f3 - f4) < tol );
pass(end+1) = ( norm(f4 - f5) < tol );

%% 'trig' and 'coeffs' in one or both dimensions
m = 31;
n = 32;
u = @(x,y) sin(2*pi*x).*cos(2*pi*y);
[xc, yc] = chebpts2(n, m);
[xt, yt] = meshgrid(trigpts(n), trigpts(m));

trigvals_chebvals     = u(xt,yc);
chebvals_trigvals     = u(xc,yt);
trigvals_trigvals     = u(xt,yt);
trigvals_chebcoeffs   = chebtech2.vals2coeffs( trigvals_chebvals     );
trigcoeffs_trigvals   =  trigtech.vals2coeffs( trigvals_trigvals.'   ).';
trigcoeffs_chebvals   =  trigtech.vals2coeffs( trigvals_chebvals.'   ).';
trigcoeffs_chebcoeffs =  trigtech.vals2coeffs( trigvals_chebcoeffs.' ).';
chebvals_trigcoeffs   =  trigtech.vals2coeffs( chebvals_trigvals     );

f1 = chebfun2(u);
f2 = chebfun2(trigcoeffs_chebvals, {'coeffs', []}, {'trig', ''});
f3 = chebfun2(trigvals_chebcoeffs, {'', 'coeffs'}, {'trig', ''});
f4 = chebfun2(trigcoeffs_chebcoeffs, 'coeffs', {'trig', []});
f5 = chebfun2(trigcoeffs_trigvals, 'trig', {'coeffs', ''});
f6 = chebfun2(chebvals_trigcoeffs, {[], 'trig'}, {[], 'coeffs'});

pass(end+1) = ( norm(f1 - f2) < tol );
pass(end+1) = ( norm(f2 - f3) < tol );
pass(end+1) = ( norm(f3 - f4) < tol );
pass(end+1) = ( norm(f4 - f5) < tol );
pass(end+1) = ( norm(f5 - f6) < tol );

%% Preferences in one or both dimensions
p = pref;
p.tech = @trigtech;

f = chebfun2(@(x,y) sin(2*pi*x).*y, {p, []});
g = chebfun2(@(x,y) sin(2*pi*x).*y);
pass(end+1) = ( norm(f - g) < tol );

f = chebfun2(@(x,y) sin(2*pi*y).*x, {[], p});
g = chebfun2(@(x,y) sin(2*pi*y).*x);
pass(end+1) = ( norm(f - g) < tol );

f = chebfun2(@(x,y) sin(2*pi*x).*sin(2*pi*y), {p, p});
g = chebfun2(@(x,y) sin(2*pi*x).*sin(2*pi*y));
pass(end+1) = ( norm(f - g) < tol );

%% The last argument takes precedence
f = chebfun2(@(x,y) sin(2*pi*x).*y, 'trig', {'trig', ''});
g = chebfun2(@(x,y) sin(2*pi*x).*y, {'trig', []});
pass(end+1) = ( norm(f - g) < tol );

m = 31;
n = 32;
u = @(x,y) sin(x.*y);
[xx, yy] = meshgrid(linspace(-1,1,n), linspace(-1,1,m));
uu = u(xx,yy);
f = chebfun2(uu, {'equi', []}, 'equi');
g = chebfun2(uu, 'equi');
pass(end+1) = ( norm(f - g) < tol );

prefx = pref;
prefx.tech = @trigtech;
f = chebfun2(@(x,y) x, {prefx, []}, pref);
g = chebfun2(@(x,y) x, pref);
pass(end+1) = ( norm(f - g) < tol );

end
