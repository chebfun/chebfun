function pass = test_sphharm( ) 
% Test spherical harmonic command.

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Linear spherical harmonics
f = spherefun.sphharm(2,1);
g = spherefun(@(x,y,z) -sqrt(15/2/pi)*x.*z);
pass(1) = norm(f-g,inf) < tol;

f = spherefun.sphharm(2,-1);
g = spherefun(@(x,y,z) sqrt(15/2/pi)*y.*z);
pass(2) = norm(f-g,inf) < tol;

f = spherefun.sphharm(3,1);
g = spherefun(@(x,y,z) 1/4*sqrt(21/pi)*x.*(1-5*z.^2));
pass(3) = norm(f-g,inf) < tol;

% Loosen the tolerance
tol = 1e4*chebfunpref().cheb2Prefs.chebfun2eps;

jj = 4; 
% Low order Y's:
for m = 1:10 
    for ll = m:m+3
        u = spherefun.sphharm(ll, m);
        v = laplacian(u);
        pass(jj) = ( norm(v + ll*(ll+1)*u, inf) < tol );
        jj = jj + 1; 
    end
end 

end