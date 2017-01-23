function pass = test_diff(pref)
% Check the diff command in Chebfun3

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e6 * pref.cheb3Prefs.chebfun3eps;
j = 1; 

% Battery: Functions, df/dx, df/dy and df/dz.
ff = {@(x,y,z) x, @(x,y,z) cos(x).*exp(y).*sin(z), @(x,y,z) cos(x.*y.*z), ...
    @(x,y,z) x.^2 + x.*y.^2.*z.^3};

ffx = {@(x,y,z) 1+0*x, @(x,y,z) -sin(x).*exp(y).*sin(z), ...
    @(x,y,z) -y.*z.*sin(x.*y.*z), @(x,y,z) 2*x+y.^2.*z.^3};

ffy = {@(x,y,z) 0*x, @(x,y,z) cos(x).*exp(y).*sin(z), ...
    @(x,y,z) -x.*z.*sin(x.*y.*z), @(x,y,z) 2*x.*y.*z.^3};

ffz = {@(x,y,z) 0*x, @(x,y,z) cos(x).*exp(y).*cos(z), ...
    @(x,y,z) -x.*y.*sin(x.*y.*z), @(x,y,z) 3*x.*y.^2.*z.^2};

% Domains:
dom = [-1 1 -1  1 -1  1;
       -3 1 -1  2  3  5];
     
% loop over the various options. 
for jj = 1:length(ff)
    for r = 1:size(dom, 1)
        f = chebfun3(ff{jj}, dom(r, :));
        fx = chebfun3(ffx{jj}, dom(r, :));
        fy = chebfun3(ffy{jj}, dom(r, :));
        fz = chebfun3(ffz{jj}, dom(r, :));

        dx = diff(f, 1, 1); 
        dy = diff(f, 1, 2);
        dz = diff(f, 1, 3);

        pass(j) = norm(dx - fx) < tol; j=j+1;
        pass(j) = norm(dy - fy) < tol; j=j+1;
        pass(j) = norm(dz - fz) < tol; j=j+1;
        
    end
end

% Check different syntax for diff. 
jj = 3;  
f = chebfun3(ff{jj}, dom(r, :)); 

err = norm(diff(f) - diff(f, 1, 1));
pass(j) = err < tol; 
j = j + 1;

err = norm(diff(f, 1) - diff(f, 1, 1));
pass(j) = err < tol; 
j = j + 1;

err = norm(diff(f, 2) - diff(f, 2, 1));
pass(j) = err < 100*tol; 
j = j + 1;

err = norm(diff(f, 2, 2) - diff(diff(f, 1, 2), 1, 2));
pass(j) = err < 10*tol;

end