function pass = test_optimization(pref)
% Can we do global optimization over [0 1]^3?

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1; 
dom = [0 1 0 1 0 1];

Battery = {@(x,y,z) cos(pi*x.*y.*z),...
    @(x,y,z) cos(2*pi*x.*y.*z), ...
    @(x,y,z) cos(3*pi*x.*y.*z),...
    @(x,y,z) cos(4*pi*x.*y.*z),...
    @(x,y,z) cos(5*pi*x.*y.*z),...
    @(x,y,z) cos(6*pi*x.*y.*z),...
    @(x,y,z) cos(7*pi*x.*y.*z),...
    @(x,y,z) sin(pi*x.*y.*z),...
    @(x,y,z) cos(0*pi*(x-y-z).^2),...
    @(x,y,z) cos(pi*(x-y-z).^2),...
    @(x,y,z) cos(2*pi*(x-y-z).^2),...
    @(x,y,z) sin(x+y+z)
    };

maxima = [1
    1
    1
    1
    1
    1
    1
    1
    1
    1
    1
    1];

minima = [-1
    -1
    -1
    -1
    -1
    -1
    -1
     0
    +1
    -1
    -1
     0];

for jj = 1:length(Battery)
    f = Battery{jj};
    g = chebfun3(f, dom);
     [Y, ignored] = minandmax3(g); 
    err(jj) = norm(Y(1) - minima(jj)) + norm(Y(2) - maxima(jj));
end

%%
pass = err < tol;

end