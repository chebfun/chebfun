function pass = test_gaussfilt( ) 
% Check that the smoothing function works.

tol = 1e3*chebfunpref().cheb2Prefs.chebfun2eps;
 
% Smoothing a constant should just give that constant back.
f = spherefun(@(lam,th) 1 + 0*lam);
g = gaussfilt(f);
pass(1) = norm(f-g) < tol;

% The smooth of a constant shouldn't depend on the smooth parameter.
f = spherefun(@(lam,th) 1 + 0*lam);
g = gaussfilt(f,2);
pass(2) = norm(f-g) < tol;

% Using a function with mean zero, check that the norm of the smoothed
% solution decreases as the smoothing parameter increases
f = spherefun.sphharm(13,7);
sig = [1 10 100];
for j=1:numel(sig)
    g = gaussfilt(f,100);
    pass(2+j) = norm(g) < norm(f);
    f = g;
end

% Check that smoothing does not change the mean of the function too much.
f = 2 + spherefun.sphharm(12,5); % Mean of this function is 2.
g = gaussfilt(f,2);
pass(j+1) = abs(mean2(g)-2) < tol;

end