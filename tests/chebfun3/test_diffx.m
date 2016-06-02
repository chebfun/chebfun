function pass = test_diffx(pref)
% Check the diffx command in Chebfun3
% [reviewed by LNT 31.05.16]

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3* pref.cheb3Prefs.chebfun3eps;
j = 1; 

% Battery:  functions, df/dx and d^2f/dx^2
ff = {@(x,y,z) x, @(x,y,z) cos(x).*exp(y).*sin(z), @(x,y,z) cos(x.*y.*z), ...
    @(x,y,z) x.^2 + x.*y.^2.*z.^3};

ffx = {@(x,y,z) 1+0*x, @(x,y,z) -sin(x).*exp(y).*sin(z), ...
    @(x,y,z) -y.*z.*sin(x.*y.*z), @(x,y,z) 2*x+y.^2.*z.^3};


ffxx = {@(x,y,z) 0*x, @(x,y,z) -cos(x).*exp(y).*sin(z), ...
    @(x,y,z) -y.^2.*z.^2.*cos(x.*y.*z), @(x,y,z) 2};

% Domains
dom = [-1 1 -1  1 -1  1;... 
       -3 1 -1  2  3  5];
     
% loop over the various options. 
for jj = 1:length(ff)
    for r = 1:size(dom, 1)
        f = chebfun3(ff{jj}, dom(r,:));
        fx = chebfun3(ffx{jj}, dom(r,:));
        fxx = chebfun3(ffxx{jj}, dom(r,:));

        dx = diffx(f, 1); 
        dx2 = diffx(f, 2); 

        pass(j) = norm(dx - fx) < tol; j=j+1;
        pass(j) = norm(dx2 - fxx) < tol; j=j+1;
        
    end
end
end