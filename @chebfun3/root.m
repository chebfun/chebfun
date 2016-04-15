function r = root(f, g, h, varargin)
% Find just one common roots of three 3D functions.

pref = chebfunpref; 
%tol = pref.eps;
tol = pref.chebfuneps;

% There are 3 steps in this algorithm:
% 1) Find an initial guess from the tensor of values of obj_fun (usually, a root with 3 accurate digits)
% 2) Turn the problem to an optimization problem: Improve that initial guess by using MATLAB's fminsearch (typically, this make the root accurate to 6-9 digits)
% Note that the value of obj_fun at the output is usually down to machine
% precision. But, this root is accurate only to 6-9 digits as a common root
% of f, g and h.
% 3) Improve the root again using Newton's method (the goal of this step is to make the root accurate to 16 digits)

% TODO: Can we remove Step 2 above?

%obj_fun = f.^2+g.^2+h.^2;
%len=length(obj_fun);
%len=max([length(f);length(g); length(h)]);
[mF, nF, pF] = length(f);
[mG, nG, pG] = length(g);
[mH, nH, pH] = length(h);
len=max([mF, nF, pF; 
         mG, nG, pG; 
         mH, nH, pH]);

%m = 3*len(1); n = 3*len(2); p = 3*len(3); 
m = 5*len(1); n = 5*len(2); p = 5*len(3); 

xx = chebpts(m, f.domain(1:2)); 
yy = chebpts(n, f.domain(3:4)); 
zz = chebpts(p, f.domain(5:6));
[xx,yy,zz] = ndgrid(xx, yy, zz);
%T = chebpolyval3(obj_fun,m,n,p);

%T = chebpolyval3(f,m,n,p).^2 + chebpolyval3(g,m,n,p).^2 + chebpolyval3(h,m,n,p).^2;

f = f./f.vscale; 
g = g./g.vscale; 
h = h./h.vscale;
T = chebpolyval3(f,m,n,p).^2 + chebpolyval3(g,m,n,p).^2 + ...
    chebpolyval3(h,m,n,p).^2;
%T2 = abs(chebpolyval3(f,m,n,p)) + abs(chebpolyval3(g,m,n,p)) + abs(chebpolyval3(h,m,n,p));

% Step 1:
[ignore, initInd]=min(abs(T(:)));
[indX, indY, indZ] = ind2sub(size(T), initInd);
initialloc = [xx(indX, indY, indZ), yy(indX, indY, indZ), ...
    zz(indX, indY, indZ)];
r = initialloc;

% % Step 2:
% options = optimset('Display', 'off', 'TolFun', eps, 'TolX', eps);
% ff = @(x) feval( obj_fun, x(1), x(2), x(3) );
% [r, val_obj_f] = fminsearch(@(x) ff(x), r, options);

% Step 3: Newton's method.
[diffF1, diffF2, diffF3] = grad(f);
[diffG1, diffG2, diffG3] = grad(g);
[diffH1, diffH2, diffH3] = grad(h);

Jac = @(x,y,z) [feval(diffF1, x, y, z),  feval(diffF2, x, y, z),  feval(diffF3, x, y, z)
                feval(diffG1, x, y, z),  feval(diffG2, x, y, z),  feval(diffG3, x, y, z)
                feval(diffH1, x, y, z),  feval(diffH2, x, y, z),  feval(diffH3, x, y, z)];
            
%r = r';
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