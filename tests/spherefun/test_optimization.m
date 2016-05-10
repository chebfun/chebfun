function pass = test_optimization( pref ) 
%% Can we do global optimization?

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1; 

Battery = {@(x,y,z) cos(pi*z),...
    @(x,y,z) sin(2*pi*z),...
    @(x,y,z) cos(pi*x), ...
    @(x,y,z) cos(pi*y),...
    @(x,y,z) cos(pi*x).*cos(pi*y),...
    @(x,y,z) cos(2*pi*x).*cos(2*pi*y),...
    @(x,y,z) exp(-10*((x-1/sqrt(2)).^2 + (z-1/sqrt(2)).^2 + y.^2))
    };

Maxi = [
    1
    1
    1
    1
    1
    1
    1
    ];

Mini = [
    -1
    -1
    -1
    -1
    -1
    -1
    0
    ];

tt = [];
for jj = 1:length(Battery)
    f = Battery{jj};
    g = spherefun(f);
    s = tic; 
    [Y, X] = minandmax2(g); 
    t = toc(s);  
    tt(jj) = t;
    err(jj) = norm(Y(1) - Mini(jj)) + norm(Y(2) - Maxi(jj));
end

%%
pass = err < tol;

end