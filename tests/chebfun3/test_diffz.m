function pass = test_diffz(pref)
% Check the diffz command in Chebfun3

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e6 * pref.cheb3Prefs.chebfun3eps;
j = 1; 

% Battery:  functions, df/dz and d^2f/dz^2
ff = {@(x,y,z) x, @(x,y,z) cos(x).*exp(y).*sin(z), @(x,y,z) cos(x.*y.*z), ...
    @(x,y,z) x.^2 + x.*y.^2.*z.^3};

ffz = {@(x,y,z) 0*x, @(x,y,z) cos(x).*exp(y).*cos(z), ...
    @(x,y,z) -x.*y.*sin(x.*y.*z), @(x,y,z) 3*x.*y.^2.*z.^2};


ffzz = {@(x,y,z) 0*x, @(x,y,z) -cos(x).*exp(y).*sin(z), ...
    @(x,y,z) -x.^2.*y.^2.*cos(x.*y.*z), @(x,y,z) 6*x.*y.^2.*z};

% Domains
dom = [-1 1 -1  1 -1  1;... 
       -3 1 -1  2  3  5];
     
% loop over the various options. 
for jj = 1:length(ff)
    for r = 1:size(dom, 1)
        f = chebfun3(ff{jj}, dom(r,:));
        fz = chebfun3(ffz{jj}, dom(r,:));
        fzz = chebfun3(ffzz{jj}, dom(r,:));

        dz = diffz(f, 1);
        dz2 = diffz(f, 2); 

        pass(j) = norm(dz - fz) < tol; j=j+1;
        pass(j) = norm(dz2 - fzz) < tol; j=j+1;
        
    end
end
end