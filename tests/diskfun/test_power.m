function pass = test_power( ) 
% Test power in diskfun 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = diskfun(@(x,y) x);
g = diskfun(@(x,y) x.^2);
pass(1) = norm( sample(f, 100, 100).^2 - ...
                sample(g, 100, 100), inf) < tol;

% Another example:
f = diskfun(@(x,y) cos(x.*y));
g = diskfun(@(x,y) cos(x.*y).^cos(x.*y));
pass(2) = ( norm(sample(f.^f, 100, 100) - sample(g, 100, 100) ) < tol );

end