function pass = test_optimization( pref ) 
%% Can we do global optimization?

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 5e4*pref.cheb2Prefs.chebfun2eps;
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

%% Check double wrapping is hidden:
j = length(pass);

if ( isempty(ver('optim')) )
  % increase tolerance for last few tests if we don't have fmincon
  tol = 1000*tol;
end

f = spherefun(@(x,y,z) exp(x));
[m,loc] = max2(f);
pass(j+1) = abs(loc(1) - 0) < tol;

f = spherefun(@(x,y,z) exp(x));
[m,loc] = min2(f);
pass(j+2) = abs(mod(loc(1),pi) - 0) < tol;
pass(j+3) = abs(m-0.367879441171442) < 1e-6;

end
