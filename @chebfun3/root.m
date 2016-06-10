function r = root(f, g, h, varargin)
%root   Find just ONE common roots of three 3D functions.
%
%   We try to find the minimum value of the following objective function:
%   objFun = f.^2 + g.^2 + h.^2.
%   There are two steps in the algorithm:
%   1) Find an initial guess from the tensor of values of obj_fun (usually,
%   a root with 3 accurate digits).
%   2) Improve the initial guess using Newton's method to make the root 
%   accurate to 16 digits).
%
%   See also chebfun2/roots.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

pref = chebfunpref; 
tol = pref.cheb3Prefs.chebfun3eps;

%% Step 1: Find an initial guess via a simulated construction of the 
% objective function f.^2 + g.^2 + h.^2.
[mF, nF, pF] = length(f);
[mG, nG, pG] = length(g);
[mH, nH, pH] = length(h);
len=max([mF, nF, pF; 
         mG, nG, pG; 
         mH, nH, pH]);
dom = f.domain;
m = 5*len(1); n = 5*len(2); p = 5*len(3); 

xx = chebpts(m, dom(1:2)); 
yy = chebpts(n, dom(3:4)); 
zz = chebpts(p, dom(5:6));
[xx,yy,zz] = ndgrid(xx, yy, zz);
f = f./f.vscale; 
g = g./g.vscale; 
h = h./h.vscale;
T = chebpolyval3(f,m,n,p).^2 + chebpolyval3(g,m,n,p).^2 + ...
    chebpolyval3(h,m,n,p).^2;

[ignore, ind]=min(abs(T(:)));
[indX, indY, indZ] = ind2sub(size(T), ind);
r = [xx(indX, indY, indZ), yy(indX, indY, indZ), zz(indX, indY, indZ)];

%% Step 2: Newton's method.
[diffF1, diffF2, diffF3] = grad(f);
[diffG1, diffG2, diffG3] = grad(g);
[diffH1, diffH2, diffH3] = grad(h);

Jac = @(x,y,z) [feval(diffF1, x, y, z),  feval(diffF2, x, y, z),  feval(diffF3, x, y, z)
                feval(diffG1, x, y, z),  feval(diffG2, x, y, z),  feval(diffG3, x, y, z)
                feval(diffH1, x, y, z),  feval(diffH2, x, y, z),  feval(diffH3, x, y, z)];
            
update = 1; 
iter = 1;
while ( norm(update) > 10*tol && iter < 10 )
    update = Jac(r(1),r(2),r(3))\ [feval(f,r(1),r(2),r(3)); 
                                   feval(g,r(1),r(2),r(3)); 
                                   feval(h,r(1),r(2),r(3))];
    r = r - update.'; 
    iter = iter + 1;
end

end

%% TODO: Could Step 2 be improved with the following alternative?
% options = optimset('Display', 'off', 'TolFun', eps, 'TolX', eps);
% ff = @(x) feval(obj_fun, x(1), x(2), x(3));
% [r, val_obj_f] = fminsearch(@(x) ff(x), r, options);