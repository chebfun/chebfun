function pass = test_sphharm( ) 
% Test spherical harmonic command.

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Linear spherical harmonics
f = spherefun.sphharm(2,1);
g = spherefun(@(x,y,z) -sqrt(15/4/pi)*x.*z);
pass(1) = norm(f-g,inf) < tol;

f = spherefun.sphharm(2,-1);
g = spherefun(@(x,y,z) sqrt(15/4/pi)*y.*z);
pass(2) = norm(f-g,inf) < tol;

f = spherefun.sphharm(3,1);
g = spherefun(@(x,y,z) 1/4*sqrt(21/2/pi)*x.*(1-5*z.^2));
pass(3) = norm(f-g,inf) < tol;

% Test that the spherical harmonics are orthonormal
pass(4) = abs(sum2(spherefun.sphharm(11,9).^2)-1) < tol;
pass(5) = abs(sum2(spherefun.sphharm(20,-8).^2)-1) < tol;

% Loosen the tolerance
tol = 1e4*chebfunpref().cheb2Prefs.chebfun2eps;

jj = 6; 
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