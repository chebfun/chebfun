function pass = test_diffy(pref)
% Check the diffy command in Chebfun3

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e6 * pref.cheb3Prefs.chebfun3eps;
j = 1; 

% Battery:  functions, df/dy and d^2f/dy^2
ff = {@(x,y,z) x, @(x,y,z) cos(x).*exp(y).*sin(z), @(x,y,z) cos(x.*y.*z), ...
    @(x,y,z) x.^2 + x.*y.^2.*z.^3};

ffy = {@(x,y,z) 0*x, @(x,y,z) cos(x).*exp(y).*sin(z), ...
    @(x,y,z) -x.*z.*sin(x.*y.*z), @(x,y,z) 2*x.*y.*z.^3};

ffyy = {@(x,y,z) 0*x, @(x,y,z) cos(x).*exp(y).*sin(z), ...
    @(x,y,z) -x.^2.*z.^2.*cos(x.*y.*z), @(x,y,z) 2*x.*z.^3};

% Domains
dom = [-1 1 -1  1 -1  1;... 
       -3 1 -1  2  3  5];
     
% loop over the various options. 
for jj = 1:length(ff)
    for r = 1:size(dom, 1)
        f = chebfun3(ff{jj}, dom(r,:));
        fy = chebfun3(ffy{jj}, dom(r,:));
        fyy = chebfun3(ffyy{jj}, dom(r,:));

        dy = diffy(f, 1); 
        dy2 = diffy(f, 2); 

        pass(j) = norm(dy - fy) < tol; j=j+1;
        pass(j) = norm(dy2 - fyy) < tol; j=j+1;
        
    end
end
end