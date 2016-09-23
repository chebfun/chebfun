function r = root(f, g, h)
%ROOT   Find just ONE common root of three CHEBFUN3 objects.
%
%   We try to find the minimum value of the following objective function:
%   objFun = f.^2 + g.^2 + h.^2.
%   There are two steps in the algorithm:
%   1) Find an initial guess from the tensor of values of obj_fun (usually,
%   a root with 3 accurate digits).
%   2) Improve the initial guess using Newton's method to make the root 
%   accurate to 16 digits).
%
% See also CHEBFUN2/ROOTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

pref = chebfunpref; 
tol = pref.cheb3Prefs.chebfun3eps;

%% Step 1: Find an initial guess via a simulated construction of the 
% objective function f.^2 + g.^2 + h.^2.
% Always use a fixed grid to sample. A bigger tensor is probably not
% possible to try and the accuracy won't be great, but at the same time 
% this is big enough to get good accuracy at least for easy functions. See 
% e.g., issue #1900.
m = 161; n = m; p = m;

dom = f.domain;
xx = chebpts(m, dom(1:2)); 
yy = chebpts(n, dom(3:4)); 
zz = chebpts(p, dom(5:6));
[xx,yy,zz] = ndgrid(xx, yy, zz);
f = f./f.vscale; 
g = g./g.vscale; 
h = h./h.vscale;
T = chebpolyval3(f, m, n, p).^2 + chebpolyval3(g, m, n, p).^2 + ...
    chebpolyval3(h, m, n, p).^2;

[~, ind] = min(abs(T(:)));
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
while ( ( norm(update) > 10*tol ) && ( iter < 10 ) )
    update = Jac(r(1), r(2), r(3)) \ [ feval(f, r(1), r(2), r(3)); 
                                       feval(g, r(1), r(2), r(3)); 
                                       feval(h, r(1), r(2), r(3)) ];
    r = r - update.'; 
    iter = iter + 1;
end

end