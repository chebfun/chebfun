function pass = test_poisson( pref )

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e7*pref.techPrefs.chebfuneps;

%% Poisson with Dirichlet BC
% Example 1:
f = ballfun(@(x,y,z)0);       % Forcing function
exact = @(r, lam, th) r.*sin(lam).*sin(th);
bc = @(lam, th) exact(1, th, lam);           % Dirichlet condition around r=1
u = poisson(f, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(1) = norm(u - exact) < tol;

%% Poisson with Neumann BC
% Example 2
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2,'spherical');
f = laplacian(exact);
bc = @(lam, th) 2*sin(th).^2;  % Neumann condition
u = poisson(f, bc, 39, 40, 41, 'neumann');
pass(2) = norm(diff(u,1) - diff(exact,1)) < tol ...
       && norm(diff(u,2) - diff(exact,2)) < tol ...
       && norm(diff(u,3) - diff(exact,3)) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

