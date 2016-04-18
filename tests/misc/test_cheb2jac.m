function pass = test_cheb2jac(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% Note: We compare a=b=.5 and a=b=.5 as these are special (easy) cases.

seedRNG(1234);
n = 100;
tol = n^2*eps;

%% scalar
C = rand(n,1);

err = cheb2jac(C, .5, .5) - cheb2jac(C, .5+eps, .5);
err = norm(err(:), inf);
pass(1) = err < tol;

err = cheb2jac(C, -.5, -.5) - cheb2jac(C, -.5+eps, -.5);
err = norm(err(:), inf);
pass(2) = err < tol;

%% matrix
C = rand(n);

err = cheb2jac(C, .5, .5) - cheb2jac(C, .5+eps, .5);
err = norm(err(:), inf);
pass(3) = err < tol;

err = cheb2jac(C, -.5, -.5) - cheb2jac(C, -.5+eps, -.5);
err = norm(err(:), inf);
pass(4) = err < tol;

%% matrix
C = rand(513);

err = cheb2jac(C, .5, .5) - cheb2jac(C, .5+eps, .5);
err = norm(err(:), inf);
pass(3) = err < tol;

err = cheb2jac(C, -.5, -.5) - cheb2jac(C, -.5+eps, -.5);
err = norm(err(:), inf);
pass(4) = err < tol;

end
